//! Random access to a reference FASTA, and the replay that uses it.
//!
//! The replay is what makes the identity numbers tool independent: each record's SEQ is walked
//! against the actual reference bases through its CIGAR and the edit distance recomputed. That is
//! necessary because protal emits no `NM` tag at all, and a tool's own `NM` comes from its own
//! CIGAR, so it cannot reveal a CIGAR that describes no real alignment.
//!
//! Random access needs an index, so `<fasta>.fai` is built on first use — by scanning the FASTA
//! here rather than shelling out to samtools, which benchpro should not require at runtime.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use log::{debug, info};

use super::cigar;
use super::error::{AlignError, AlignResult};
use super::sam::{AlnRecord, ParseCounters};
use super::truth::ReadKey;

/// One `.fai` record: where a sequence's bases live in the FASTA.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FaiEntry {
    /// Sequence length in bases.
    pub length: u64,
    /// Byte offset of the first base.
    pub offset: u64,
    /// Bases per line.
    pub linebases: u64,
    /// Bytes per line, including the terminator.
    pub linewidth: u64,
}

impl FaiEntry {
    /// Byte offset of a given base.
    ///
    /// # Arguments
    ///
    /// * `base` - 0-based base index within the sequence
    ///
    /// # Returns
    ///
    /// The absolute byte offset in the FASTA.
    pub fn byte_of(&self, base: u64) -> u64 {
        self.offset + base / self.linebases * self.linewidth + base % self.linebases
    }
}

/// Path of the index beside a FASTA.
fn fai_path(fasta: &Path) -> PathBuf {
    let mut path = fasta.as_os_str().to_os_string();
    path.push(".fai");
    PathBuf::from(path)
}

/// Returns the path of `<fasta>.fai`, building it if absent or stale.
///
/// Written to a temporary file and renamed, so a killed run cannot leave a truncated index that a
/// later run would trust.
///
/// # Arguments
///
/// * `fasta` - Reference FASTA
///
/// # Returns
///
/// Path of the usable index.
///
/// # Errors
///
/// Returns [`AlignError::Io`] when the FASTA cannot be read or the index cannot be written.
pub fn ensure_fai(fasta: &Path) -> AlignResult<PathBuf> {
    let fai = fai_path(fasta);
    if is_fresh(&fai, fasta) {
        return Ok(fai);
    }

    info!("Indexing {} (no .fai found)", fasta.display());
    let file = File::open(fasta).map_err(|e| AlignError::io(fasta, e))?;
    let mut reader = BufReader::with_capacity(1 << 20, file);
    let mut index = String::new();

    let mut name: Option<String> = None;
    let (mut length, mut offset, mut linebases, mut linewidth) = (0u64, 0u64, 0u64, 0u64);
    let mut position = 0u64;
    let mut line = Vec::new();

    loop {
        line.clear();
        let read = reader
            .read_until(b'\n', &mut line)
            .map_err(|e| AlignError::io(fasta, e))?;
        if read == 0 {
            break;
        }
        if line.first() == Some(&b'>') {
            if let Some(previous) = name.take() {
                index.push_str(&format!(
                    "{previous}\t{length}\t{offset}\t{linebases}\t{linewidth}\n"
                ));
            }
            let header = &line[1..];
            let end = header
                .iter()
                .position(|b| b.is_ascii_whitespace())
                .unwrap_or(header.len());
            name = Some(String::from_utf8_lossy(&header[..end]).into_owned());
            length = 0;
            offset = position + read as u64;
            linebases = 0;
            linewidth = 0;
        } else {
            let bases = line.iter().filter(|b| !b.is_ascii_whitespace()).count() as u64;
            if linebases == 0 {
                linebases = bases;
                linewidth = read as u64;
            }
            length += bases;
        }
        position += read as u64;
    }
    if let Some(previous) = name {
        index.push_str(&format!(
            "{previous}\t{length}\t{offset}\t{linebases}\t{linewidth}\n"
        ));
    }

    // The temp name carries this process's id: two benchpro runs indexing the same reference at
    // once would otherwise share one temp path, and the second rename would find it already
    // consumed by the first. Rename is atomic, so whichever finishes last simply wins.
    let mut tmp = fai.as_os_str().to_os_string();
    tmp.push(format!(".tmp.{}", std::process::id()));
    let tmp = PathBuf::from(tmp);
    std::fs::write(&tmp, index).map_err(|e| AlignError::io(&tmp, e))?;
    std::fs::rename(&tmp, &fai).map_err(|e| AlignError::io(&fai, e))?;
    Ok(fai)
}

/// Is an existing index at least as new as the FASTA it indexes?
fn is_fresh(fai: &Path, fasta: &Path) -> bool {
    let modified = |path: &Path| std::fs::metadata(path).and_then(|m| m.modified()).ok();
    match (modified(fai), modified(fasta)) {
        (Some(index), Some(source)) => index >= source,
        _ => false,
    }
}

/// Random access to a reference FASTA through its `.fai`.
///
/// The index is filtered to the contigs the alignments actually name: a marker reference has
/// 14.5 M sequences and loading all of them would cost gigabytes for a scoring run that touches a
/// small fraction of them.
#[derive(Debug)]
pub struct Reference {
    file: File,
    index: HashMap<Box<str>, FaiEntry>,
}

impl Reference {
    /// Opens a reference, loading only the wanted contigs from its index.
    ///
    /// # Arguments
    ///
    /// * `fasta` - Reference FASTA; its `.fai` is built if missing
    /// * `wanted` - Contigs to keep in the index; `None` keeps every one
    ///
    /// # Returns
    ///
    /// A handle that can fetch reference bases.
    ///
    /// # Errors
    ///
    /// Returns [`AlignError::Io`] when the FASTA or index cannot be read, and
    /// [`AlignError::Parse`] when an index line is malformed.
    pub fn open(fasta: &Path, wanted: Option<&HashSet<Box<str>>>) -> AlignResult<Self> {
        let fai = ensure_fai(fasta)?;
        let handle = File::open(&fai).map_err(|e| AlignError::io(&fai, e))?;
        let mut index = HashMap::new();

        for (i, line) in BufReader::new(handle).lines().enumerate() {
            let line = line.map_err(|e| AlignError::io(&fai, e))?;
            if line.trim().is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 5 {
                return Err(AlignError::parse(
                    &fai,
                    i + 1,
                    "expected 5 tab separated fields (name, length, offset, linebases, linewidth)",
                ));
            }
            if wanted.is_some_and(|set| !set.contains(fields[0])) {
                continue;
            }
            let number = |value: &str, what: &str| -> AlignResult<u64> {
                value.parse::<u64>().map_err(|_| {
                    AlignError::parse(&fai, i + 1, format!("{what} '{value}' is not a number"))
                })
            };
            index.insert(
                fields[0].into(),
                FaiEntry {
                    length: number(fields[1], "length")?,
                    offset: number(fields[2], "offset")?,
                    linebases: number(fields[3], "linebases")?.max(1),
                    linewidth: number(fields[4], "linewidth")?.max(1),
                },
            );
        }

        debug!(
            "{}: {} contig(s) loaded from the index",
            fasta.display(),
            index.len()
        );

        Ok(Self {
            file: File::open(fasta).map_err(|e| AlignError::io(fasta, e))?,
            index,
        })
    }

    /// Whether a contig is present in the loaded index.
    ///
    /// # Arguments
    ///
    /// * `name` - Contig name
    ///
    /// # Returns
    ///
    /// `true` when the reference declares it.
    pub fn contains(&self, name: &str) -> bool {
        self.index.contains_key(name)
    }

    /// Number of contigs in the loaded index.
    ///
    /// # Returns
    ///
    /// The contig count.
    pub fn len(&self) -> usize {
        self.index.len()
    }

    /// Whether the loaded index is empty.
    ///
    /// # Returns
    ///
    /// `true` when no contig was loaded.
    pub fn is_empty(&self) -> bool {
        self.index.is_empty()
    }

    /// The length of a contig.
    ///
    /// # Arguments
    ///
    /// * `name` - Contig name
    ///
    /// # Returns
    ///
    /// Its length in bases, or `None` when it is unknown.
    pub fn length_of(&self, name: &str) -> Option<u64> {
        self.index.get(name).map(|entry| entry.length)
    }

    /// Every loaded contig's length.
    ///
    /// # Returns
    ///
    /// A map from contig name to length, for callers that need geometry rather than bases.
    pub fn lengths(&self) -> HashMap<Box<str>, u64> {
        self.index
            .iter()
            .map(|(name, entry)| (name.clone(), entry.length))
            .collect()
    }

    /// Fetches reference bases `[start, end)`, uppercased and clipped to the contig.
    ///
    /// # Arguments
    ///
    /// * `name` - Contig name
    /// * `start` - 0-based inclusive start
    /// * `end` - 0-based exclusive end
    ///
    /// # Returns
    ///
    /// The bases, or an empty vector when the contig is unknown or the range is empty. A short
    /// return is meaningful: it means the reference does not have the bases the record claims.
    ///
    /// # Errors
    ///
    /// Returns [`AlignError::Io`] when the FASTA cannot be read.
    pub fn fetch(&mut self, name: &str, start: u64, end: u64) -> AlignResult<Vec<u8>> {
        let Some(entry) = self.index.get(name).copied() else {
            return Ok(Vec::new());
        };
        let end = end.min(entry.length);
        if start >= end {
            return Ok(Vec::new());
        }

        let from = entry.byte_of(start);
        let to = entry.byte_of(end - 1) + 1;
        self.file
            .seek(SeekFrom::Start(from))
            .map_err(|e| AlignError::io("<reference>", e))?;

        // Read up to the computed span rather than insisting on it. The `.fai` records one line
        // width per sequence, so a FASTA whose lines are not uniformly wrapped can put `to` past
        // the end of the file -- `read_exact` would abort the whole run over a short final line,
        // where a short read is exactly the "the reference does not have those bases" signal the
        // replay already knows how to charge.
        let mut buffer = vec![0u8; (to - from) as usize];
        let mut filled = 0;
        while filled < buffer.len() {
            match self.file.read(&mut buffer[filled..]) {
                Ok(0) => break,
                Ok(n) => filled += n,
                Err(e) => return Err(AlignError::io("<reference>", e)),
            }
        }
        buffer.truncate(filled);

        buffer.retain(|b| !b.is_ascii_whitespace());
        buffer.make_ascii_uppercase();
        Ok(buffer)
    }
}

/// Replays the retained alignments against the reference, filling in their `vnm`.
///
/// Records are fetched in reference order: the marker FASTA is 16.6 GB, and sorted seeks keep the
/// reads local. That makes the replay inherently sequential per file; parallelism belongs across
/// files, not inside one.
///
/// An alignment whose target is not in the reference at all is left unverified and counted
/// separately — that is a naming problem, not a bad alignment.
///
/// # Arguments
///
/// * `records` - Parsed alignments; only those retaining a CIGAR and SEQ are replayed
/// * `reference` - The opened reference
/// * `counters` - Counters to record unknown targets in
///
/// # Returns
///
/// How many alignments were replayed.
///
/// # Errors
///
/// Returns [`AlignError::Io`] when the reference cannot be read.
pub fn verify(
    records: &mut HashMap<ReadKey, AlnRecord>,
    reference: &mut Reference,
    counters: &mut ParseCounters,
) -> AlignResult<usize> {
    let mut order: Vec<ReadKey> = records
        .iter()
        .filter(|(_, record)| record.seq.is_some() && record.cigar.is_some())
        .map(|(key, _)| key.clone())
        .collect();
    order.sort_by(|a, b| {
        let (left, right) = (&records[a], &records[b]);
        (&left.target, left.pos0, a).cmp(&(&right.target, right.pos0, b))
    });

    let mut verified = 0;
    for key in order {
        let (target, pos0, span) = {
            let record = &records[&key];
            (
                record.target.clone(),
                record.pos0,
                record.counts.reference_span(),
            )
        };
        if !reference.contains(&target) {
            counters.unknown_ref += 1;
            continue;
        }
        let bases = reference.fetch(&target, pos0, pos0 + span)?;
        let record = records.get_mut(&key).expect("key came from this map");
        let (Some(cigar), Some(seq)) = (record.cigar.as_ref(), record.seq.as_ref()) else {
            continue;
        };
        record.vnm = Some(cigar::replay(cigar.as_bytes(), seq.as_bytes(), &bases));
        verified += 1;
    }

    Ok(verified)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::cigar::CigarCounts;
    use std::io::Write;

    fn temp_dir(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("benchpro_align_ref_{name}"));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn write_fasta(dir: &Path, content: &str) -> PathBuf {
        let path = dir.join("ref.fna");
        File::create(&path)
            .unwrap()
            .write_all(content.as_bytes())
            .unwrap();
        path
    }

    /// Two contigs, wrapped at 10 bases per line.
    const WRAPPED: &str =
        ">ctg1 some description\nACGTACGTAC\nGTGTGTGTGT\nAAAA\n>ctg2\nTTTTTTTTTT\n";

    #[test]
    fn builds_an_index_that_locates_every_base() {
        let dir = temp_dir("build");
        let fasta = write_fasta(&dir, WRAPPED);
        let mut reference = Reference::open(&fasta, None).unwrap();

        assert_eq!(reference.len(), 2);
        assert_eq!(reference.length_of("ctg1"), Some(24));
        assert_eq!(reference.length_of("ctg2"), Some(10));
        // The description after the first whitespace is not part of the name.
        assert!(reference.contains("ctg1"));

        assert_eq!(reference.fetch("ctg1", 0, 4).unwrap(), b"ACGT");
        // Spanning the line break must not pick up the newline.
        assert_eq!(reference.fetch("ctg1", 8, 12).unwrap(), b"ACGT");
        assert_eq!(reference.fetch("ctg1", 20, 24).unwrap(), b"AAAA");
        assert_eq!(reference.fetch("ctg2", 0, 10).unwrap(), b"TTTTTTTTTT");
    }

    #[test]
    fn fetching_past_the_end_returns_what_exists() {
        let dir = temp_dir("clip");
        let fasta = write_fasta(&dir, WRAPPED);
        let mut reference = Reference::open(&fasta, None).unwrap();

        // A short return is the signal the replay needs: the reference does not have those bases.
        assert_eq!(reference.fetch("ctg1", 20, 40).unwrap(), b"AAAA");
        assert!(reference.fetch("ctg1", 30, 40).unwrap().is_empty());
        assert!(reference.fetch("missing", 0, 10).unwrap().is_empty());
    }

    #[test]
    fn an_unwrapped_fasta_indexes_too() {
        let dir = temp_dir("unwrapped");
        let fasta = write_fasta(&dir, ">only\nACGTACGTACGTACGT\n");
        let mut reference = Reference::open(&fasta, None).unwrap();

        assert_eq!(reference.length_of("only"), Some(16));
        assert_eq!(reference.fetch("only", 4, 8).unwrap(), b"ACGT");
    }

    #[test]
    fn a_fasta_without_a_trailing_newline_still_indexes() {
        let dir = temp_dir("no_newline");
        let fasta = write_fasta(&dir, ">only\nACGTACGTAC\nGTGT");
        let mut reference = Reference::open(&fasta, None).unwrap();

        assert_eq!(reference.length_of("only"), Some(14));
        assert_eq!(reference.fetch("only", 10, 14).unwrap(), b"GTGT");
    }

    #[test]
    fn lowercase_reference_bases_are_normalised() {
        let dir = temp_dir("lowercase");
        let fasta = write_fasta(&dir, ">only\nacgtacgtac\n");
        let mut reference = Reference::open(&fasta, None).unwrap();
        assert_eq!(reference.fetch("only", 0, 4).unwrap(), b"ACGT");
    }

    #[test]
    fn the_index_is_reused_and_can_be_filtered() {
        let dir = temp_dir("reuse");
        let fasta = write_fasta(&dir, WRAPPED);
        Reference::open(&fasta, None).unwrap();
        assert!(
            fai_path(&fasta).exists(),
            "the index is written beside the FASTA"
        );

        let wanted: HashSet<Box<str>> = ["ctg2".into()].into_iter().collect();
        let filtered = Reference::open(&fasta, Some(&wanted)).unwrap();
        assert_eq!(
            filtered.len(),
            1,
            "only the wanted contig is held in memory"
        );
        assert!(!filtered.contains("ctg1"));
    }

    fn record(target: &str, pos0: u64, cigar: &str, seq: &str) -> AlnRecord {
        AlnRecord {
            target: target.into(),
            pos0,
            mapq: 60,
            nm: None,
            counts: cigar::count(cigar.as_bytes()),
            clip_ends: cigar::clip_ends(cigar.as_bytes()),
            malformed: false,
            proper_pair: false,
            cigar: Some(cigar.into()),
            seq: Some(seq.into()),
            vnm: None,
            offset: 0,
        }
    }

    fn key(id: &str) -> ReadKey {
        ReadKey {
            id: id.into(),
            mate: 1,
        }
    }

    #[test]
    fn verify_fills_in_the_replayed_edit_distance() {
        let dir = temp_dir("verify");
        let fasta = write_fasta(&dir, WRAPPED);
        let mut reference = Reference::open(&fasta, None).unwrap();
        let mut counters = ParseCounters::default();

        let mut records: HashMap<ReadKey, AlnRecord> = [
            // exact match against ACGTACGT
            (key("perfect"), record("ctg1", 0, "8M", "ACGTACGT")),
            // one mismatch in the middle
            (key("one_off"), record("ctg1", 0, "8M", "ACGTAAGT")),
            // a target the reference has never heard of
            (key("ghost"), record("nowhere", 0, "8M", "ACGTACGT")),
        ]
        .into_iter()
        .collect();

        let verified = verify(&mut records, &mut reference, &mut counters).unwrap();

        assert_eq!(verified, 2);
        assert_eq!(records[&key("perfect")].vnm, Some(0));
        assert_eq!(records[&key("one_off")].vnm, Some(1));
        assert_eq!(records[&key("ghost")].vnm, None);
        assert_eq!(
            counters.unknown_ref, 1,
            "an unknown target is a naming problem, not a bad alignment"
        );
    }

    #[test]
    fn records_without_a_retained_sequence_are_skipped() {
        let dir = temp_dir("no_seq");
        let fasta = write_fasta(&dir, WRAPPED);
        let mut reference = Reference::open(&fasta, None).unwrap();
        let mut counters = ParseCounters::default();

        let mut dropped = record("ctg1", 0, "8M", "ACGTACGT");
        dropped.cigar = None;
        dropped.seq = None;
        let mut records: HashMap<ReadKey, AlnRecord> =
            [(key("dropped"), dropped)].into_iter().collect();

        assert_eq!(
            verify(&mut records, &mut reference, &mut counters).unwrap(),
            0
        );
        assert_eq!(counters.unknown_ref, 0);
    }

    #[test]
    fn a_cigar_walking_past_the_contig_end_is_charged_for_it() {
        let dir = temp_dir("overrun");
        let fasta = write_fasta(&dir, ">only\nACGT\n");
        let mut reference = Reference::open(&fasta, None).unwrap();
        let mut counters = ParseCounters::default();

        // Claims 8 aligned bases on a 4-base contig: the 4 the reference cannot supply are
        // mismatches, not free.
        let mut records: HashMap<ReadKey, AlnRecord> =
            [(key("over"), record("only", 0, "8M", "ACGTACGT"))]
                .into_iter()
                .collect();
        verify(&mut records, &mut reference, &mut counters).unwrap();

        assert_eq!(records[&key("over")].vnm, Some(4));
    }

    #[test]
    fn a_ragged_final_line_does_not_abort_the_run() {
        // The .fai records the FIRST line's width, so a shorter final line makes the computed span
        // overshoot the file. That must degrade to a short read, not an error.
        let dir = temp_dir("ragged");
        let fasta = write_fasta(&dir, ">only\nACGTACGTAC\nACG\n");
        let mut reference = Reference::open(&fasta, None).unwrap();

        assert_eq!(reference.length_of("only"), Some(13));
        assert_eq!(reference.fetch("only", 0, 13).unwrap(), b"ACGTACGTACACG");
        assert_eq!(reference.fetch("only", 10, 13).unwrap(), b"ACG");
    }

    #[test]
    fn byte_offsets_follow_the_fai_formula() {
        let entry = FaiEntry {
            length: 100,
            offset: 10,
            linebases: 10,
            linewidth: 11,
        };
        assert_eq!(entry.byte_of(0), 10);
        assert_eq!(entry.byte_of(9), 19);
        assert_eq!(entry.byte_of(10), 21, "the newline is skipped");
        assert_eq!(entry.byte_of(20), 32);
    }

    #[test]
    fn counts_helper_is_consistent_with_the_record_span() {
        let counts: CigarCounts = cigar::count(b"10S80M5D");
        assert_eq!(counts.reference_span(), 85);
    }
}
