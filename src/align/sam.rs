//! Parsing an aligner's output into the records scoring needs.
//!
//! SAM goes through [`bioreader`]'s multithreaded reader; PAF has a small reader here, because
//! bioreader has none and PAF only ever feeds the mapping-level score.
//!
//! Two properties matter more than speed:
//!
//! * **Determinism.** The first record for a mate is its primary, and with several worker threads
//!   "first" is not a thread-local notion — so records are reconciled by their byte offset in the
//!   file. A benchmark whose numbers wobble between runs is a bug generator.
//! * **Bounded memory.** A SAM can hold tens of millions of records and their sequences dwarf
//!   everything else scoring needs, so CIGAR and SEQ are retained only for a reservoir sample of
//!   `verify_limit` records — reservoir, not the first N, because SAM order follows read order,
//!   which on simulated data is grouped by source genome.

use std::collections::{BinaryHeap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use bioreader::parallel::fastq::Merge;
use bioreader::parallel::sam::read_sam_state_par_ids;
use bioreader::sequence::sam_record::RefSamRecord;
use flate2::read::MultiGzDecoder;
use log::debug;

use super::cigar::{self, CigarCounts};
use super::error::{AlignError, AlignResult};
use super::meta::{is_gzipped, AlignmentFormat};
use super::truth::ReadKey;
use serde::{Deserialize, Serialize};

/// Batch size handed to each bioreader worker.
const BUFFER_SIZE: usize = 1 << 24;

/// One primary alignment, reduced to what scoring needs.
#[derive(Debug, Clone)]
pub struct AlnRecord {
    /// Reference sequence the read was placed on.
    pub target: Box<str>,
    /// 0-based leftmost position, converted from the file's 1-based POS.
    pub pos0: u64,
    /// Mapping quality.
    pub mapq: u8,
    /// The tool's own `NM` tag, or `None` when it emits none (protal).
    pub nm: Option<i64>,
    /// CIGAR base counts.
    pub counts: CigarCounts,
    /// Leading and trailing clip lengths, kept for every record because the clip geometry needs
    /// which end was clipped and the CIGAR itself is only retained for a sample.
    pub clip_ends: (u64, u64),
    /// Fingerprint of the CIGAR, kept for every record so an alignment can be compared against a
    /// gold standard without retaining the string. `None` when the format carries no CIGAR at all
    /// (PAF) — a sentinel value would make "no alignment to compare" indistinguishable from "an
    /// alignment that happens to hash to it".
    pub cigar_fingerprint: Option<u64>,
    /// SEQ length disagrees with what the CIGAR consumes: the record cannot be interpreted.
    pub malformed: bool,
    /// The record carries the proper-pair flag. Stored per record rather than counted as it is
    /// seen, because a record that loses the primary contest must not contribute to the totals.
    pub proper_pair: bool,
    /// CIGAR, retained only for records selected for the reference replay.
    pub cigar: Option<Box<str>>,
    /// SEQ, retained only for records selected for the reference replay.
    pub seq: Option<Box<str>>,
    /// Edit distance recomputed against the reference, filled in by the replay.
    pub vnm: Option<u64>,
    /// Byte offset of the record in its file, used to pick the primary deterministically.
    pub offset: u64,
}

impl AlnRecord {
    /// Query coverage: the share of the read that is actually aligned.
    ///
    /// # Returns
    ///
    /// `(matches + ins) / read_len`, or 0 for a zero-length read.
    pub fn coverage(&self) -> f64 {
        let read_len = self.counts.read_len();
        if read_len == 0 {
            0.0
        } else {
            (self.counts.matches + self.counts.ins) as f64 / read_len as f64
        }
    }

    /// Identity from an edit distance: `1 - nm/aln_len`, floored at 0.
    ///
    /// # Arguments
    ///
    /// * `nm` - An edit distance, reported or replayed
    ///
    /// # Returns
    ///
    /// The identity, or `None` when there is no edit distance or no alignment length.
    pub fn identity_of(&self, nm: Option<i64>) -> Option<f64> {
        let nm = nm?;
        let aln_len = self.counts.aln_len();
        if aln_len == 0 {
            return None;
        }
        Some((1.0 - nm as f64 / aln_len as f64).max(0.0))
    }

    /// Identity as the tool reports it.
    ///
    /// # Returns
    ///
    /// `1 - NM/aln_len`, or `None` when the record carries no `NM` tag.
    pub fn identity(&self) -> Option<f64> {
        self.identity_of(self.nm)
    }

    /// Identity from the edit distance recomputed against the reference.
    ///
    /// # Returns
    ///
    /// The replayed identity, or `None` when this record was not replayed.
    pub fn identity_verified(&self) -> Option<f64> {
        self.identity_of(self.vnm.map(|v| v as i64))
    }

    /// Whether the alignment contains an insertion or a deletion.
    ///
    /// # Returns
    ///
    /// `true` when the CIGAR has an `I`, `D` or `N` operation.
    pub fn has_indel(&self) -> bool {
        self.counts.ins > 0 || self.counts.del > 0
    }
}

/// Counts of everything seen while parsing, including what was skipped.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ParseCounters {
    /// Records read, header lines excluded.
    pub records: u64,
    /// Records flagged unmapped.
    pub unmapped: u64,
    /// Secondary or supplementary records.
    pub secondary: u64,
    /// Mapped records whose CIGAR is `*`.
    pub no_cigar: u64,
    /// Primary alignments without an `NM` tag.
    pub no_nm: u64,
    /// Primary alignments flagged as a proper pair.
    pub proper_pair: u64,
    /// Alignments whose target is absent from the reference; filled in by the replay.
    pub unknown_ref: u64,
}

impl ParseCounters {
    /// Folds another set of counters into this one, for pooled aggregation.
    ///
    /// # Arguments
    ///
    /// * `other` - The counters to add
    pub fn add(&mut self, other: &Self) {
        self.records += other.records;
        self.unmapped += other.unmapped;
        self.secondary += other.secondary;
        self.no_cigar += other.no_cigar;
        self.no_nm += other.no_nm;
        self.proper_pair += other.proper_pair;
        self.unknown_ref += other.unknown_ref;
    }
}

/// A parsed alignment file.
#[derive(Debug, Clone, Default)]
pub struct ParsedAlignment {
    /// Primary alignment per read mate.
    pub records: HashMap<ReadKey, AlnRecord>,
    /// What was seen and skipped on the way.
    pub counters: ParseCounters,
}

/// A record admitted to the replay sample, ordered by its hash so the worst is evicted first.
///
/// `Ord` is by hash alone descending is not needed: [`BinaryHeap`] is a max-heap, so the entry with
/// the *largest* hash sits on top and is the one a better candidate displaces.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Sampled {
    hash: u64,
    key: ReadKey,
}

impl Ord for Sampled {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // The key breaks ties so the ordering is total and does not depend on insertion order.
        self.hash
            .cmp(&other.hash)
            .then_with(|| self.key.cmp(&other.key))
    }
}

impl PartialOrd for Sampled {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// Worker-local parse state. `Merge` folds the workers together at the end.
#[derive(Debug, Clone, Default)]
struct ParseState {
    records: HashMap<ReadKey, AlnRecord>,
    counters: ParseCounters,
    /// How many records may retain a sequence; 0 means all of them.
    keep_limit: usize,
    /// Seed mixed into the sampling hash.
    seed: u64,
    /// The records currently holding a sequence, as a max-heap on the sampling hash.
    reservoir: BinaryHeap<Sampled>,
    /// The keys currently in `reservoir`. The heap cannot answer "is this key already in?" in less
    /// than linear time, and it has to: the same mate can reach `admit` twice (see it there).
    sampled: HashSet<ReadKey>,
}

impl Merge for ParseState {
    fn merge_from(&mut self, other: &mut Self) {
        self.counters.add(&other.counters);
        for (key, record) in other.records.drain() {
            match self.records.get(&key) {
                // Both workers saw this mate. The primary is the record that comes first in the
                // file, which is what makes the result independent of how batches were split.
                Some(existing) if existing.offset <= record.offset => {}
                _ => {
                    self.records.insert(key, record);
                }
            }
        }
        // Records merge first, so an eviction here always finds the record it must strip.
        for candidate in std::mem::take(&mut other.reservoir) {
            self.admit(candidate);
        }
    }
}

impl ParseState {
    /// Offers a record to the replay sample, evicting the current worst if it is beaten.
    ///
    /// Admission is decided by [`sample_hash`] alone, so the surviving set is the `keep_limit`
    /// records with the smallest hash — a pure function of the record set and the seed. That is
    /// what makes the sample independent of how batches happened to be split across workers.
    fn admit(&mut self, candidate: Sampled) {
        // A mate can be offered twice: two workers each saw a primary record for it (a malformed or
        // merged SAM), and both put it in their reservoir. Admitting it again would spend a second
        // slot on the same read, and worse, losing the second round would strip the bases off the
        // copy the heap is still holding -- so the record would sit in the sample unreplayable.
        if self.sampled.contains(&candidate.key) {
            return;
        }
        // The record this candidate speaks for may have lost the primary contest during the merge,
        // in which case the surviving record is a different one that may already have been
        // stripped. Sampling it would put an unreplayable entry in the heap and shrink the sample.
        if !self
            .records
            .get(&candidate.key)
            .is_some_and(|record| record.seq.is_some())
        {
            return;
        }
        if self.keep_limit == 0 || self.reservoir.len() < self.keep_limit {
            self.sampled.insert(candidate.key.clone());
            self.reservoir.push(candidate);
            return;
        }
        // The heap is full: keep the candidate only if it beats the current worst.
        match self.reservoir.peek() {
            Some(worst) if candidate < *worst => {
                let evicted = self.reservoir.pop().expect("heap is non-empty");
                self.sampled.remove(&evicted.key);
                self.strip(&evicted.key);
                self.sampled.insert(candidate.key.clone());
                self.reservoir.push(candidate);
            }
            _ => self.strip(&candidate.key),
        }
    }

    /// Drops a record's retained sequence, keeping its metrics.
    fn strip(&mut self, key: &ReadKey) {
        if let Some(record) = self.records.get_mut(key) {
            record.cigar = None;
            record.seq = None;
        }
    }
}

/// The sampling hash of a read mate.
///
/// A pure function of the key and the seed — never of arrival order, worker id or thread count.
/// Uses a fixed mixing function rather than [`std::collections::hash_map::DefaultHasher`], whose
/// `RandomState` is seeded per process and would vary run to run.
fn sample_hash(key: &ReadKey, seed: u64) -> u64 {
    let mut state = seed ^ 0x9E37_79B9_7F4A_7C15;
    for &byte in key.id.as_bytes() {
        state = mix(state ^ byte as u64);
    }
    mix(state ^ (key.mate as u64) << 56)
}

/// One round of the SplitMix64 finaliser.
fn mix(value: u64) -> u64 {
    let mut z = value.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Reads the `@SQ` reference-sequence lengths from a SAM header.
///
/// A fallback for contig lengths when no reference FASTA is given. The header is at the top of the
/// file, so this stops at the first alignment record rather than reading the whole SAM.
///
/// # Arguments
///
/// * `path` - SAM file, optionally gzipped
///
/// # Returns
///
/// Contig name to length, empty when the file carries no `@SQ` lines.
///
/// # Errors
///
/// Returns [`AlignError::Io`] when the file cannot be read.
pub fn header_lengths(path: &Path) -> AlignResult<HashMap<Box<str>, u64>> {
    let reader = BufReader::with_capacity(1 << 20, open(path)?);
    let mut lengths = HashMap::new();

    for line in reader.lines() {
        let line = line.map_err(|e| AlignError::io(path, e))?;
        if !line.starts_with('@') {
            break;
        }
        if !line.starts_with("@SQ") {
            continue;
        }
        let (mut name, mut length) = (None, None);
        for field in line.split('\t') {
            if let Some(value) = field.strip_prefix("SN:") {
                name = Some(value);
            } else if let Some(value) = field.strip_prefix("LN:") {
                length = value.parse::<u64>().ok();
            }
        }
        if let (Some(name), Some(length)) = (name, length) {
            lengths.insert(name.into(), length);
        }
    }

    Ok(lengths)
}

/// Parses an alignment file into primary records.
///
/// # Arguments
///
/// * `path` - Alignment file, optionally gzipped
/// * `format` - Which parser to use, taken from the samplesheet rather than probed
/// * `keep_seq` - Retain CIGAR and SEQ so records can be replayed against a reference
/// * `keep_limit` - How many records may retain a sequence (0 = all)
/// * `threads` - Worker threads for SAM parsing; 0 or 1 means single threaded
/// * `seed` - Seed for the reservoir, so a run is reproducible
///
/// # Returns
///
/// The primary alignment of each read mate, plus the counters describing what was skipped.
///
/// # Errors
///
/// Returns [`AlignError::Io`] when the file cannot be opened or read.
pub fn parse_alignment(
    path: &Path,
    format: AlignmentFormat,
    keep_seq: bool,
    keep_limit: usize,
    threads: usize,
    seed: u64,
) -> AlignResult<ParsedAlignment> {
    match format {
        AlignmentFormat::Sam => parse_sam(path, keep_seq, keep_limit, threads, seed),
        AlignmentFormat::Paf => parse_paf(path),
    }
}

/// Opens an alignment file, transparently decompressing a `.gz`.
///
/// `MultiGzDecoder` rather than `GzDecoder`: a bgzipped SAM is a multi-member gzip stream, and
/// `GzDecoder` stops after the first member — truncating the file without ever erroring.
fn open(path: &Path) -> AlignResult<Box<dyn Read + Send>> {
    let file = File::open(path).map_err(|e| AlignError::io(path, e))?;
    if is_gzipped(path) {
        Ok(Box::new(MultiGzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

/// Fills in the counters that can only be known once the primary of every mate is settled.
///
/// `no_nm` and `proper_pair` describe the records that *won* the primary contest, so counting them
/// as records stream past would count a duplicate that is later discarded — and which duplicate a
/// worker sees depends on how batches were split, which is not deterministic.
fn finish_counters(parsed: &mut ParsedAlignment) {
    parsed.counters.no_nm = parsed
        .records
        .values()
        .filter(|record| record.nm.is_none())
        .count() as u64;
    parsed.counters.proper_pair = parsed
        .records
        .values()
        .filter(|record| record.proper_pair)
        .count() as u64;
}

/// Resolves the requested thread count to the number of workers to spawn.
///
/// `0` means "all available", which is what the CLI documents — taking `max(1)` of it instead would
/// silently make the *default* single threaded, which is the one setting nobody passes explicitly
/// and therefore the one nobody would notice.
///
/// # Arguments
///
/// * `requested` - Threads asked for; 0 means all available
///
/// # Returns
///
/// At least 1, falling back to 1 when the available parallelism cannot be determined.
fn worker_threads(requested: usize) -> u32 {
    if requested > 0 {
        return requested as u32;
    }
    std::thread::available_parallelism()
        .map(|n| n.get() as u32)
        .unwrap_or(1)
}

/// Parses a SAM file through bioreader's multithreaded reader.
fn parse_sam(
    path: &Path,
    keep_seq: bool,
    keep_limit: usize,
    threads: usize,
    seed: u64,
) -> AlignResult<ParsedAlignment> {
    let reader = open(path)?;
    let num_threads = worker_threads(threads);

    let initial = ParseState {
        keep_limit: if keep_seq { keep_limit } else { 0 },
        seed,
        ..Default::default()
    };

    let state = read_sam_state_par_ids(
        reader,
        BUFFER_SIZE,
        num_threads,
        initial,
        move |record: &RefSamRecord, read_id, state: &mut ParseState| {
            consume_record(record, read_id.byte_offset, keep_seq, state);
        },
    )
    .map_err(|e| AlignError::io(path, e))?;

    debug!(
        "{}: {} primary alignments, {} unmapped, {} secondary",
        path.display(),
        state.records.len(),
        state.counters.unmapped,
        state.counters.secondary
    );

    let mut parsed = ParsedAlignment {
        records: state.records,
        counters: state.counters,
    };
    finish_counters(&mut parsed);
    Ok(parsed)
}

/// Folds one SAM record into the worker's state.
fn consume_record(record: &RefSamRecord, offset: u64, keep_seq: bool, state: &mut ParseState) {
    state.counters.records += 1;

    let Some(flag) = record.try_flag() else {
        return;
    };
    if flag & 0x4 != 0 {
        state.counters.unmapped += 1;
        return;
    }
    if flag & 0x100 != 0 || flag & 0x800 != 0 {
        state.counters.secondary += 1;
        return;
    }
    let cigar = record.cigar();
    if cigar == b"*" || cigar.is_empty() {
        state.counters.no_cigar += 1;
        return;
    }

    let key = ReadKey::from_sam(record.qname(), flag);
    // The first record for a mate is its primary; a later one at a higher offset never wins.
    if let Some(existing) = state.records.get(&key) {
        if existing.offset <= offset {
            return;
        }
    }

    let counts = cigar::count(cigar);
    let clip_ends = cigar::clip_ends(cigar);
    let cigar_fingerprint = Some(cigar::fingerprint(cigar));
    let seq = record.seq();
    // SAM requires SEQ length to equal the CIGAR's query consumption (soft clips count, hard do
    // not). A mismatch means the record cannot be interpreted: emitted, but broken.
    let malformed = seq != b"*" && seq.len() as u64 != counts.query_consumed();

    // `no_nm` and `proper_pair` are NOT counted here. This record may still lose the primary
    // contest to a lower-offset duplicate -- in another worker, whose batch this one cannot see --
    // and a counter incremented now could not be taken back. Both are derived from the surviving
    // record set once every worker has merged (see `finish_counters`).
    let nm = record.tag(b"NM").and_then(|tag| tag.as_int());

    let mut aln = AlnRecord {
        target: String::from_utf8_lossy(record.rname()).into_owned().into(),
        pos0: record.try_pos().unwrap_or(0).saturating_sub(1),
        mapq: record.try_mapq().unwrap_or(0),
        nm,
        counts,
        clip_ends,
        cigar_fingerprint,
        malformed,
        proper_pair: flag & 0x2 != 0,
        cigar: None,
        seq: None,
        vnm: None,
        offset,
    };

    // Retain the bases optimistically; `admit` strips them again if this record does not make the
    // sample. The record must already be in the map when that happens, so insert first.
    if keep_seq && seq != b"*" {
        aln.cigar = Some(String::from_utf8_lossy(cigar).into_owned().into());
        aln.seq = Some(String::from_utf8_lossy(seq).into_owned().into());
        let hash = sample_hash(&key, state.seed);
        state.records.insert(key.clone(), aln);
        state.admit(Sampled { hash, key });
        return;
    }

    state.records.insert(key, aln);
}

/// Parses a PAF file: mapping level only, since PAF carries no CIGAR, SEQ or NM.
///
/// Columns used are 1 (query name, mate from its `/N` suffix), 6 (target), 8 (target start, already
/// 0-based) and 12 (MAPQ). The first line for a read is its primary.
fn parse_paf(path: &Path) -> AlignResult<ParsedAlignment> {
    let reader = BufReader::new(open(path)?);
    let mut parsed = ParsedAlignment::default();

    for (i, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| AlignError::io(path, e))?;
        let line = line.trim_end_matches('\r');
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }
        parsed.counters.records += 1;

        let Some(key) = ReadKey::from_paf(fields[0]) else {
            continue;
        };
        if parsed.records.contains_key(&key) {
            continue;
        }
        let Ok(pos0) = fields[7].parse::<u64>() else {
            continue;
        };
        let mapq = fields
            .get(11)
            .and_then(|f| f.parse::<u8>().ok())
            .unwrap_or(0);

        parsed.records.insert(
            key,
            AlnRecord {
                target: fields[5].into(),
                pos0,
                mapq,
                nm: None,
                counts: CigarCounts::default(),
                clip_ends: (0, 0),
                cigar_fingerprint: None,
                malformed: false,
                proper_pair: false,
                cigar: None,
                seq: None,
                vnm: None,
                offset: i as u64,
            },
        );
    }

    finish_counters(&mut parsed);
    Ok(parsed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn temp_file(name: &str, content: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(format!("benchpro_align_sam_{name}"));
        let mut file = File::create(&path).unwrap();
        file.write_all(content.as_bytes()).unwrap();
        path
    }

    const SAM: &str = concat!(
        "@HD\tVN:1.6\n",
        "@SQ\tSN:ctg1\tLN:1000\n",
        // mate 1, clean 8M with NM:i:1
        "r1\t99\tctg1\t101\t60\t8M\t=\t201\t108\tACGTACGT\tIIIIIIII\tNM:i:1\n",
        // mate 2 of the same pair, soft clipped
        "r1\t147\tctg1\t201\t60\t2S6M\t=\t101\t-108\tTTACGTAC\tIIIIIIII\tNM:i:0\n",
        // a secondary record that must not be scored
        "r2\t256\tctg1\t301\t0\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tNM:i:0\n",
        // r2's primary, no NM tag at all (protal-style)
        "r2/1\t0\tctg1\t401\t30\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n",
        // unmapped
        "r3\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\n",
        // malformed: CIGAR consumes 8 query bases, SEQ has 4
        "r4\t0\tctg1\t501\t20\t8M\t*\t0\t0\tACGT\tIIII\tNM:i:0\n",
    );

    /// Each test writes its own file: they run concurrently, and a shared name means one test
    /// reads another's content.
    fn parse_named(name: &str, content: &str, threads: usize) -> ParsedAlignment {
        let path = temp_file(&format!("{name}.sam"), content);
        parse_alignment(&path, AlignmentFormat::Sam, true, 0, threads, 0).unwrap()
    }

    #[test]
    fn keeps_only_primary_mapped_records() {
        let parsed = parse_named("keeps_only", SAM, 1);

        assert_eq!(parsed.records.len(), 4); // r1/1, r1/2, r2/1, r4/1
        assert_eq!(parsed.counters.unmapped, 1);
        assert_eq!(parsed.counters.secondary, 1);
        assert_eq!(parsed.counters.no_nm, 1);
        assert_eq!(parsed.counters.proper_pair, 2);
    }

    #[test]
    fn positions_become_zero_based_and_mates_come_from_the_flag() {
        let parsed = parse_named("positions", SAM, 1);

        let mate1 = &parsed.records[&ReadKey {
            id: "r1".into(),
            mate: 1,
        }];
        assert_eq!(mate1.pos0, 100);
        assert_eq!(mate1.nm, Some(1));
        assert_eq!(mate1.mapq, 60);

        let mate2 = &parsed.records[&ReadKey {
            id: "r1".into(),
            mate: 2,
        }];
        assert_eq!(mate2.pos0, 200);
        assert_eq!(mate2.counts.soft, 2);
    }

    #[test]
    fn the_slash_suffix_convention_joins_with_the_stripped_one() {
        let parsed = parse_named("suffix", SAM, 1);
        // "r2/1" in the file must key as ("r2", 1), the same as a tool that strips the suffix.
        assert!(parsed.records.contains_key(&ReadKey {
            id: "r2".into(),
            mate: 1
        }));
    }

    #[test]
    fn missing_nm_is_recorded_rather_than_assumed_zero() {
        let parsed = parse_named("no_nm", SAM, 1);
        let r2 = &parsed.records[&ReadKey {
            id: "r2".into(),
            mate: 1,
        }];
        assert_eq!(r2.nm, None);
        assert_eq!(r2.identity(), None);
    }

    #[test]
    fn malformed_records_are_flagged_not_dropped() {
        let parsed = parse_named("malformed", SAM, 1);
        let r4 = &parsed.records[&ReadKey {
            id: "r4".into(),
            mate: 1,
        }];
        assert!(r4.malformed);
        assert!(
            !parsed.records[&ReadKey {
                id: "r1".into(),
                mate: 1
            }]
                .malformed
        );
    }

    #[test]
    fn identity_and_coverage_follow_the_python_definitions() {
        let parsed = parse_named("identity", SAM, 1);
        let mate1 = &parsed.records[&ReadKey {
            id: "r1".into(),
            mate: 1,
        }];
        assert_eq!(mate1.identity(), Some(1.0 - 1.0 / 8.0));
        assert_eq!(mate1.coverage(), 1.0);

        let mate2 = &parsed.records[&ReadKey {
            id: "r1".into(),
            mate: 2,
        }];
        assert_eq!(mate2.coverage(), 6.0 / 8.0);
    }

    #[test]
    fn thread_count_does_not_change_the_result() {
        // The whole point of reconciling by byte offset. Batches split differently per thread
        // count, and the record set must not notice.
        let single = parse_named("threads1", SAM, 1);
        let many = parse_named("threads4", SAM, 4);

        assert_eq!(single.records.len(), many.records.len());
        for (key, record) in &single.records {
            let other = &many.records[key];
            assert_eq!(record.pos0, other.pos0);
            assert_eq!(record.nm, other.nm);
            assert_eq!(record.target, other.target);
        }
    }

    #[test]
    fn later_duplicate_records_never_displace_the_primary() {
        let content = concat!(
            "r1\t0\tctgA\t101\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tNM:i:0\n",
            "r1\t0\tctgB\t201\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tNM:i:9\n",
        );
        let parsed = parse_named("duplicate", content, 1);
        assert_eq!(parsed.records.len(), 1);
        assert_eq!(
            &*parsed.records[&ReadKey {
                id: "r1".into(),
                mate: 1
            }]
                .target,
            "ctgA"
        );
    }

    #[test]
    fn reservoir_bounds_retained_sequences_without_losing_records() {
        let mut content = String::new();
        for i in 0..50 {
            content.push_str(&format!(
                "r{i}\t0\tctg1\t{}\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tNM:i:0\n",
                i + 1
            ));
        }
        let path = temp_file("reservoir.sam", &content);
        let parsed = parse_alignment(&path, AlignmentFormat::Sam, true, 10, 1, 0).unwrap();

        assert_eq!(parsed.records.len(), 50, "every record is still scored");
        let with_seq = parsed.records.values().filter(|r| r.seq.is_some()).count();
        assert_eq!(with_seq, 10, "only the budget retains a sequence");
    }

    #[test]
    fn reservoir_sampling_is_reproducible_for_a_seed() {
        let mut content = String::new();
        for i in 0..50 {
            content.push_str(&format!(
                "r{i}\t0\tctg1\t{}\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tNM:i:0\n",
                i + 1
            ));
        }
        let path = temp_file("seeded.sam", &content);
        let keys = |seed| {
            let parsed = parse_alignment(&path, AlignmentFormat::Sam, true, 10, 1, seed).unwrap();
            let mut keys: Vec<ReadKey> = parsed
                .records
                .iter()
                .filter(|(_, r)| r.seq.is_some())
                .map(|(k, _)| k.clone())
                .collect();
            keys.sort();
            keys
        };
        assert_eq!(keys(7), keys(7));
    }

    #[test]
    fn gzipped_sam_reads_the_same_as_plain() {
        use flate2::write::GzEncoder;
        use flate2::Compression;

        let path = std::env::temp_dir().join("benchpro_align_sam_gz.sam.gz");
        let mut encoder = GzEncoder::new(File::create(&path).unwrap(), Compression::fast());
        encoder.write_all(SAM.as_bytes()).unwrap();
        encoder.finish().unwrap();

        let parsed = parse_alignment(&path, AlignmentFormat::Sam, true, 0, 1, 0).unwrap();
        assert_eq!(parsed.records.len(), 4);
    }

    /// Inserts a record holding bases, then offers it to the sample -- the same order
    /// `consume_record` uses, since `admit` only samples records that are present and unstripped.
    fn insert_and_admit(state: &mut ParseState, id: &str, offset: u64) {
        let key = ReadKey {
            id: id.into(),
            mate: 1,
        };
        state.records.insert(
            key.clone(),
            AlnRecord {
                target: "ctg1".into(),
                pos0: 0,
                mapq: 60,
                nm: Some(0),
                counts: cigar::count(b"4M"),
                clip_ends: (0, 0),
                cigar_fingerprint: None,
                malformed: false,
                proper_pair: false,
                cigar: Some("4M".into()),
                seq: Some("ACGT".into()),
                vnm: None,
                offset,
            },
        );
        let hash = sample_hash(&key, state.seed);
        state.admit(Sampled { hash, key });
    }

    /// Runs `admit` over the given keys in the given order, returning the surviving sample.
    fn admit_all(order: &[u32], limit: usize) -> Vec<ReadKey> {
        let mut state = ParseState {
            keep_limit: limit,
            seed: 7,
            ..Default::default()
        };
        for i in order {
            insert_and_admit(&mut state, &format!("read{i}"), 0);
        }
        let mut kept: Vec<ReadKey> = state.reservoir.into_iter().map(|s| s.key).collect();
        kept.sort();
        kept
    }

    /// A state holding `keys` in its reservoir, each with a retained sequence.
    fn state_with(keys: &[&str], limit: usize, offset: u64) -> ParseState {
        let mut state = ParseState {
            keep_limit: limit,
            seed: 7,
            ..Default::default()
        };
        for id in keys {
            insert_and_admit(&mut state, id, offset);
        }
        state
    }

    #[test]
    fn a_mate_seen_by_two_workers_keeps_its_sequence() {
        // The same (id, mate) can appear as a primary twice in a malformed or merged SAM, and the
        // two copies can land in different workers. The merged sample must still hold it exactly
        // once, with its bases intact -- stripping it would silently shrink the replay set.
        let mut left = state_with(&["dup", "a"], 4, 0);
        let mut right = state_with(&["dup", "b"], 4, 100);

        left.merge_from(&mut right);

        let dup = ReadKey {
            id: "dup".into(),
            mate: 1,
        };
        assert!(
            left.records[&dup].seq.is_some(),
            "the duplicate lost its bases and can no longer be replayed"
        );
        let slots: Vec<&ReadKey> = left.reservoir.iter().map(|s| &s.key).collect();
        assert_eq!(
            slots.iter().filter(|k| ***k == dup).count(),
            1,
            "the duplicate occupies two of the {} sample slots",
            left.keep_limit
        );
    }

    #[test]
    fn a_duplicate_does_not_evict_a_legitimate_sample_member() {
        // keep_limit is 2 and there are only 2 distinct mates, so both must survive.
        let mut left = state_with(&["dup", "a"], 2, 0);
        let mut right = state_with(&["dup"], 2, 100);

        left.merge_from(&mut right);

        let retained = left.records.values().filter(|r| r.seq.is_some()).count();
        assert_eq!(retained, 2, "a duplicate cost a distinct mate its place");
    }

    #[test]
    fn the_replay_sample_does_not_depend_on_arrival_order() {
        // Which worker sees a record, and in what order, is decided by batch scheduling -- so the
        // sample has to be a function of the records alone, or the replayed metrics wobble between
        // runs at the same seed.
        let forward: Vec<u32> = (0..500).collect();
        let backward: Vec<u32> = (0..500).rev().collect();
        let interleaved: Vec<u32> = (0..500).map(|i| (i * 37) % 500).collect();

        let expected = admit_all(&forward, 50);
        assert_eq!(expected.len(), 50);
        assert_eq!(admit_all(&backward, 50), expected);
        assert_eq!(admit_all(&interleaved, 50), expected);
    }

    #[test]
    fn the_sample_is_the_lowest_hashing_records() {
        // Admission is bottom-k on the hash, which is what makes it order independent.
        let kept = admit_all(&(0..200).collect::<Vec<_>>(), 20);
        let threshold = kept
            .iter()
            .map(|k| sample_hash(k, 7))
            .max()
            .expect("non-empty");
        for i in 0..200u32 {
            let key = ReadKey {
                id: format!("read{i}").into(),
                mate: 1,
            };
            if sample_hash(&key, 7) < threshold {
                assert!(kept.contains(&key), "read{i} hashes low but was not kept");
            }
        }
    }

    #[test]
    fn a_different_seed_selects_a_different_sample() {
        let with_seed = |seed: u64| {
            let mut state = ParseState {
                keep_limit: 20,
                seed,
                ..Default::default()
            };
            for i in 0..200u32 {
                insert_and_admit(&mut state, &format!("read{i}"), 0);
            }
            let mut kept: Vec<ReadKey> = state.reservoir.into_iter().map(|s| s.key).collect();
            kept.sort();
            kept
        };
        assert_ne!(with_seed(1), with_seed(2));
        assert_eq!(with_seed(1), with_seed(1));
    }

    #[test]
    fn mates_of_one_read_hash_independently() {
        let mate1 = ReadKey {
            id: "r1".into(),
            mate: 1,
        };
        let mate2 = ReadKey {
            id: "r1".into(),
            mate: 2,
        };
        assert_ne!(sample_hash(&mate1, 0), sample_hash(&mate2, 0));
    }

    #[test]
    fn a_discarded_duplicate_does_not_reach_the_counters() {
        // The duplicate loses the primary contest, so its missing NM and its absent proper-pair
        // flag must not appear in the totals. Counting as records stream past would include them,
        // and which worker sees which copy is not deterministic.
        let content = concat!(
            "r1\t99\tctg1\t101\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tNM:i:0\n",
            "r1\t0\tctg1\t201\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n",
        );
        let parsed = parse_named("dup_counters", content, 1);

        assert_eq!(parsed.records.len(), 1);
        assert_eq!(parsed.counters.no_nm, 0, "the winning record carries an NM");
        assert_eq!(
            parsed.counters.proper_pair, 1,
            "the winning record is flagged 0x2; the duplicate is not"
        );
    }

    #[test]
    fn counters_describe_the_winning_records_only() {
        // Same mates, opposite order in the file: the record that wins changes, so the counters
        // must change with it rather than reflecting whatever streamed past.
        let content = concat!(
            "r1\t0\tctg1\t201\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n",
            "r1\t99\tctg1\t101\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tNM:i:0\n",
        );
        let parsed = parse_named("dup_counters_rev", content, 1);

        assert_eq!(parsed.counters.no_nm, 1, "now the NM-less record wins");
        assert_eq!(parsed.counters.proper_pair, 0);
    }

    #[test]
    fn thread_zero_means_all_available_not_one() {
        let available = std::thread::available_parallelism()
            .map(|n| n.get() as u32)
            .unwrap_or(1);
        assert_eq!(worker_threads(0), available);
        assert_eq!(worker_threads(1), 1);
        assert_eq!(worker_threads(8), 8);
    }

    #[test]
    fn header_lengths_come_from_the_sq_lines() {
        let path = temp_file("header.sam", SAM);
        let lengths = header_lengths(&path).unwrap();

        assert_eq!(lengths.len(), 1);
        assert_eq!(lengths[&Box::from("ctg1")], 1000);
    }

    #[test]
    fn a_sam_without_a_header_yields_no_lengths() {
        let path = temp_file(
            "headerless.sam",
            "r1\t0\tctg1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n",
        );
        assert!(header_lengths(&path).unwrap().is_empty());
    }

    #[test]
    fn paf_takes_the_mate_from_the_suffix_and_the_position_as_is() {
        let content = concat!(
            "r1/1\t150\t0\t150\t+\tctg1\t1000\t100\t250\t150\t150\t60\n",
            "r1/2\t150\t0\t150\t-\tctg1\t1000\t300\t450\t150\t150\t7\n",
            // no mate suffix: cannot be joined, so it is skipped rather than guessed
            "r2\t150\t0\t150\t+\tctg1\t1000\t500\t650\t150\t150\t60\n",
        );
        let path = temp_file("basic.paf", content);
        let parsed = parse_alignment(&path, AlignmentFormat::Paf, false, 0, 1, 0).unwrap();

        assert_eq!(parsed.records.len(), 2);
        let mate1 = &parsed.records[&ReadKey {
            id: "r1".into(),
            mate: 1,
        }];
        assert_eq!(mate1.pos0, 100, "PAF start is already 0-based");
        assert_eq!(mate1.mapq, 60);
        assert_eq!(
            parsed.records[&ReadKey {
                id: "r1".into(),
                mate: 2
            }]
                .mapq,
            7
        );
    }
}
