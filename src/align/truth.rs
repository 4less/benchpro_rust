//! Per-read ground truth, and the key that joins it to an alignment.
//!
//! The truth TSV is the one written by the flexalign benchmark's `build_truth.py` /
//! `build_truth_protal.py`: headerless, five tab separated columns
//!
//! ```text
//! read_id <TAB> mate(1|2) <TAB> true_contig <TAB> true_pos(1-based) <TAB> true_genome
//! ```
//!
//! On a marker database `true_genome` is protal's `TIID` and `true_contig`/`true_pos` refer to the
//! read's *source genome*, not to a marker — so only the genome field is meaningful there. That is
//! exactly what [`ScoringMode::Species`](crate::options::ScoringMode::Species) exists for.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::sync::Arc;

use crate::options::ScoringMode;
use log::warn;

use super::error::{AlignError, AlignResult};
use super::meta::AlignmentFormat;
use super::metrics::ScoringContext;
use super::sam::parse_alignment;

/// A read mate: the join key between an alignment record and its truth.
///
/// The FASTQ header convention is `<id>/<mate>`. Some aligners keep that suffix (minibwa, protal)
/// and some strip it (flexalign, bowtie2), so the id is always stored *without* it and the mate
/// number comes from the SAM FLAG. Get this wrong and every metric is silently zero.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ReadKey {
    /// Read id with any trailing `/1` or `/2` removed.
    pub id: Box<str>,
    /// Mate number, 1 or 2.
    pub mate: u8,
}

impl ReadKey {
    /// Builds a key from a raw QNAME and a SAM FLAG.
    ///
    /// # Arguments
    ///
    /// * `qname` - QNAME field, possibly with a `/1` or `/2` suffix
    /// * `flag` - SAM FLAG; `0x40` means mate 1, `0x80` mate 2, neither means unpaired (mate 1)
    ///
    /// # Returns
    ///
    /// The key an alignment record joins on.
    pub fn from_sam(qname: &[u8], flag: u16) -> Self {
        let mate = if flag & 0x80 != 0 { 2 } else { 1 };
        Self {
            id: String::from_utf8_lossy(strip_mate_suffix(qname))
                .into_owned()
                .into(),
            mate,
        }
    }

    /// Builds a key from a PAF query name, which carries the mate only in its `/N` suffix.
    ///
    /// # Arguments
    ///
    /// * `qname` - PAF column 1
    ///
    /// # Returns
    ///
    /// The key, or `None` when the name has no `/1` or `/2` suffix to take the mate from.
    pub fn from_paf(qname: &str) -> Option<Self> {
        let (id, mate) = qname.rsplit_once('/')?;
        let mate: u8 = mate.parse().ok()?;
        if mate != 1 && mate != 2 {
            return None;
        }
        Some(Self {
            id: id.into(),
            mate,
        })
    }
}

/// Removes a trailing `/1` or `/2` from a read name.
fn strip_mate_suffix(qname: &[u8]) -> &[u8] {
    match qname {
        [head @ .., b'/', b'1' | b'2'] => head,
        other => other,
    }
}

/// Where a read truly came from.
///
/// `contig` and `genome` are shared rather than owned per row: a truth file names a few hundred
/// distinct contigs and genomes across millions of reads, so interning them turns two allocations
/// per row into two per distinct value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TruthEntry {
    /// Source contig (of the source genome, on a marker DB).
    pub contig: Arc<str>,
    /// 0-based source position, converted from the file's 1-based column.
    pub pos0: u64,
    /// Source genome label — a genome accession, or a protal `TIID` on a marker DB.
    pub genome: Arc<str>,
    /// Fingerprint of the gold-standard CIGAR, when the truth came from a SAM. `None` for a truth
    /// TSV, which records where a read came from but not how it should align.
    pub cigar_fingerprint: Option<u64>,
    /// What the gold alignment looks like, when the truth came from a SAM. Lets recovery be
    /// reported for the hard reads separately from the easy ones.
    pub shape: Option<GoldShape>,
}

/// Properties of a gold-standard alignment that make a read harder to place.
///
/// A tool's overall accuracy is dominated by the easy majority — ungapped, unclipped reads that
/// almost anything aligns correctly. Reporting recovery over these subsets separately is what
/// shows whether a tool actually handles the hard cases or merely averages well.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct GoldShape {
    /// The true alignment contains an insertion or deletion.
    pub indel: bool,
    /// The true alignment is clipped at one or both ends.
    pub clipped: bool,
}

/// The per-read truth of one sample.
pub type Truth = HashMap<ReadKey, TruthEntry>;

/// Loads a truth TSV.
///
/// # Arguments
///
/// * `path` - Headerless five-column truth TSV
///
/// # Returns
///
/// A map from read mate to its true origin. Later rows win over earlier duplicates, matching the
/// Python loader.
///
/// # Errors
///
/// Returns [`AlignError::Io`] when the file cannot be read and [`AlignError::Parse`] naming the
/// line when a row has the wrong number of fields or an unparseable position.
pub fn load_truth(path: &Path) -> AlignResult<Truth> {
    // One read of the whole file, then byte slices into it. `BufReader::lines()` would allocate a
    // String per line, which on a ten-million-read truth is ten million allocations that are all
    // thrown away again.
    let mut file = File::open(path).map_err(|e| AlignError::io(path, e))?;
    let mut bytes = Vec::new();
    file.read_to_end(&mut bytes)
        .map_err(|e| AlignError::io(path, e))?;

    let mut truth = Truth::default();
    let mut pool: HashMap<&[u8], Arc<str>> = HashMap::new();

    for (i, raw) in bytes.split(|b| *b == b'\n').enumerate() {
        let line = match raw.split_last() {
            Some((b'\r', head)) => head,
            _ => raw,
        };
        if line.is_empty() {
            continue;
        }

        let mut fields = line.split(|b| *b == b'\t');
        let (Some(id), Some(mate), Some(contig), Some(pos), Some(genome)) = (
            fields.next(),
            fields.next(),
            fields.next(),
            fields.next(),
            fields.next(),
        ) else {
            // A `#` line that is not a record is a comment; a `#` line that IS a record is a read
            // whose id happens to start with `#`, and dropping it would lose data silently. The
            // format is headerless, so shape decides, not the first byte.
            if line[0] == b'#' {
                continue;
            }
            return Err(AlignError::parse(
                path,
                i + 1,
                format!(
                    "expected 5 tab separated fields (read_id, mate, contig, pos, genome), found {}",
                    line.split(|b| *b == b'\t').count()
                ),
            ));
        };

        let mate = match mate {
            b"1" => 1u8,
            b"2" => 2u8,
            // A commented-out header (`# read_id<TAB>mate<TAB>...`) has five fields but no valid
            // mate. It is still a comment, not a broken record.
            _ if line[0] == b'#' => continue,
            other => {
                return Err(AlignError::parse(
                    path,
                    i + 1,
                    format!("mate '{}' is not 1 or 2", String::from_utf8_lossy(other)),
                ))
            }
        };
        let Some(parsed_pos) = parse_u64(pos) else {
            if line[0] == b'#' {
                continue;
            }
            return Err(AlignError::parse(
                path,
                i + 1,
                format!(
                    "position '{}' is not a number",
                    String::from_utf8_lossy(pos)
                ),
            ));
        };
        let pos = parsed_pos;

        truth.insert(
            ReadKey {
                id: String::from_utf8_lossy(id).into_owned().into(),
                mate,
            },
            TruthEntry {
                cigar_fingerprint: None,
                shape: None,
                contig: intern(&mut pool, contig),
                // The file is 1-based; everything downstream is 0-based. A 0 means "no position",
                // so it stays 0 rather than wrapping to u64::MAX.
                pos0: pos.saturating_sub(1),
                genome: intern(&mut pool, genome),
            },
        );
    }

    Ok(truth)
}

/// Returns the shared copy of `value`, creating it on first sight.
///
/// The pool keys on the slice into the file buffer, so a repeat costs a lookup and an `Arc` clone
/// rather than an allocation.
fn intern<'a>(pool: &mut HashMap<&'a [u8], Arc<str>>, value: &'a [u8]) -> Arc<str> {
    if let Some(shared) = pool.get(value) {
        return shared.clone();
    }
    let shared: Arc<str> = String::from_utf8_lossy(value).into_owned().into();
    pool.insert(value, shared);
    pool[value].clone()
}

/// Parses an unsigned decimal, returning `None` on anything that is not all digits.
fn parse_u64(bytes: &[u8]) -> Option<u64> {
    if bytes.is_empty() {
        return None;
    }
    let mut value: u64 = 0;
    for &byte in bytes {
        if !byte.is_ascii_digit() {
            return None;
        }
        value = value.checked_mul(10)?.checked_add((byte - b'0') as u64)?;
    }
    Some(value)
}

/// Loads truth from a gold-standard SAM.
///
/// Read simulators emit one — `art_illumina -sam` is the usual source — recording where every read
/// really came from *and* how it really aligns. That last part is what a truth TSV cannot express,
/// and it is what makes the third benchmark possible: not just "did the read land in the right
/// place" but "is the reported alignment the one the simulator produced".
///
/// The genome label must be written in whatever vocabulary the scorer will compare against, or
/// every read is judged wrong:
///
/// * `full` scoring maps the *predicted* contig through `contig2genome`, falling back to `"NA"`, so
///   the truth label is that same lookup on the gold contig — which is exactly what
///   `build_truth.py` writes. A contig missing from the map yields `"NA"` on both sides, so a read
///   on an unmapped contig is scored as landing on the right genome rather than the wrong one.
/// * `species` scoring compares the predicted contig's prefix before `sep`, so the label is the
///   gold contig's prefix. A `contig2genome` is not consulted: the scorer never looks at one in
///   this mode, and labelling the truth from a map the scorer ignores is how the two vocabularies
///   drift apart.
///
/// # Arguments
///
/// * `path` - Gold-standard `.sam` or `.sam.gz`
/// * `format` - Which parser to use
/// * `contig2genome` - Optional contig-to-genome map, consulted only under `full` scoring
/// * `scoring` - Which vocabulary the genome label must be written in
/// * `sep` - Marker-contig separator, for `species` scoring
/// * `threads` - Worker threads for parsing
///
/// # Returns
///
/// A map from read mate to its true origin, carrying the gold CIGAR's fingerprint.
///
/// # Errors
///
/// Returns [`AlignError::Io`] when the file cannot be read.
pub fn load_truth_sam(
    path: &Path,
    format: AlignmentFormat,
    contig2genome: Option<&HashMap<Box<str>, Box<str>>>,
    scoring: ScoringMode,
    sep: &str,
    threads: usize,
) -> AlignResult<Truth> {
    let context = ScoringContext {
        scoring,
        sep,
        contig2genome,
        tolerance: 0,
    };
    let parsed = parse_alignment(path, format, false, 0, threads, 0)?;
    if parsed.counters.unmapped > 0 || parsed.counters.no_cigar > 0 {
        // build_truth.py keeps these; dropping them shrinks `total` and therefore every "of total"
        // percentage, so it must not happen quietly.
        warn!(
            "{}: {} unmapped and {} CIGAR-less gold record(s) carry no true locus and are not part              of the truth",
            path.display(),
            parsed.counters.unmapped,
            parsed.counters.no_cigar
        );
    }
    let mut truth = Truth::default();
    let mut pool: HashMap<Box<str>, Arc<str>> = HashMap::new();

    for (key, record) in parsed.records {
        // The scorer's own derivation, not a copy of it: keeping a second copy here is exactly
        // what let the truth and the scorer drift into different vocabularies before.
        let genome = context.label_of(&record.target);
        let genome = match pool.get(genome) {
            Some(shared) => shared.clone(),
            None => {
                let shared: Arc<str> = genome.into();
                pool.insert(genome.into(), shared.clone());
                shared
            }
        };
        let contig = match pool.get(&*record.target) {
            Some(shared) => shared.clone(),
            None => {
                let shared: Arc<str> = (*record.target).into();
                pool.insert(record.target.clone(), shared.clone());
                shared
            }
        };

        truth.insert(
            key,
            TruthEntry {
                contig,
                pos0: record.pos0,
                genome,
                cigar_fingerprint: record.cigar_fingerprint,
                shape: record.cigar_fingerprint.map(|_| GoldShape {
                    indel: record.has_indel(),
                    clipped: record.counts.clip() > 0,
                }),
            },
        );
    }

    Ok(truth)
}

/// Loads truth from whichever format the path names.
///
/// # Arguments
///
/// * `path` - Truth TSV, or a gold-standard `.sam`/`.paf`
/// * `contig2genome` - Optional contig-to-genome map, used only for a gold-standard SAM
/// * `scoring` - Scoring mode, used only for a gold-standard SAM
/// * `sep` - Marker-contig separator, used only for a gold-standard SAM
/// * `threads` - Worker threads, used only for a gold-standard SAM
///
/// # Returns
///
/// The sample's per-read truth.
///
/// # Errors
///
/// Returns [`AlignError::Io`] or [`AlignError::Parse`] as the underlying loader does.
pub fn load_truth_any(
    path: &Path,
    contig2genome: Option<&HashMap<Box<str>, Box<str>>>,
    scoring: ScoringMode,
    sep: &str,
    threads: usize,
) -> AlignResult<Truth> {
    match AlignmentFormat::from_path(path) {
        Some(format) => load_truth_sam(path, format, contig2genome, scoring, sep, threads),
        None => load_truth(path),
    }
}

/// Loads a `contig<TAB>genome` map.
///
/// # Arguments
///
/// * `path` - Headerless two-column TSV; extra columns are ignored
///
/// # Returns
///
/// A map from contig name to the genome it belongs to.
///
/// # Errors
///
/// Returns [`AlignError::Io`] when the file cannot be read and [`AlignError::Parse`] when a row
/// has fewer than two fields.
pub fn load_contig2genome(path: &Path) -> AlignResult<HashMap<Box<str>, Box<str>>> {
    let file = File::open(path).map_err(|e| AlignError::io(path, e))?;
    let mut map = HashMap::new();

    for (i, line) in BufReader::new(file).lines().enumerate() {
        let line = line.map_err(|e| AlignError::io(path, e))?;
        let line = line.trim_end_matches('\r');
        if line.is_empty() {
            continue;
        }
        // A contig name never begins with '#', so here the first byte is decisive -- unlike the
        // truth file, where a read id might.
        if line.starts_with('#') {
            continue;
        }
        let mut fields = line.split('\t');
        let (Some(contig), Some(genome)) = (fields.next(), fields.next()) else {
            return Err(AlignError::parse(
                path,
                i + 1,
                "expected 2 tab separated fields (contig, genome)",
            ));
        };
        map.insert(contig.into(), genome.into());
    }

    Ok(map)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn temp_file(name: &str, content: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(format!("benchpro_align_truth_{name}"));
        let mut file = File::create(&path).unwrap();
        file.write_all(content.as_bytes()).unwrap();
        path
    }

    #[test]
    fn sam_key_strips_the_mate_suffix_and_takes_the_mate_from_the_flag() {
        // The two conventions must land on the same key, or the join silently finds nothing.
        let stripped = ReadKey::from_sam(b"read1", 0x80);
        let suffixed = ReadKey::from_sam(b"read1/2", 0x80);
        assert_eq!(stripped, suffixed);
        assert_eq!(stripped.mate, 2);
        assert_eq!(&*stripped.id, "read1");
    }

    #[test]
    fn sam_key_treats_an_unpaired_read_as_mate_one() {
        assert_eq!(ReadKey::from_sam(b"read1", 0).mate, 1);
        assert_eq!(ReadKey::from_sam(b"read1", 0x40).mate, 1);
    }

    #[test]
    fn a_slash_that_is_not_a_mate_suffix_is_kept() {
        let key = ReadKey::from_sam(b"lane1/tile/read9", 0x40);
        assert_eq!(&*key.id, "lane1/tile/read9");
    }

    #[test]
    fn paf_key_comes_from_the_suffix_alone() {
        assert_eq!(
            ReadKey::from_paf("read1/2"),
            Some(ReadKey {
                id: "read1".into(),
                mate: 2
            })
        );
        assert_eq!(ReadKey::from_paf("read1"), None);
        assert_eq!(ReadKey::from_paf("read1/3"), None);
    }

    #[test]
    fn truth_positions_become_zero_based() {
        let path = temp_file(
            "basic.tsv",
            "r1\t1\tctg1\t100\tgenomeA\nr1\t2\tctg1\t250\tgenomeA\n",
        );
        let truth = load_truth(&path).unwrap();

        assert_eq!(truth.len(), 2);
        let entry = &truth[&ReadKey {
            id: "r1".into(),
            mate: 1,
        }];
        assert_eq!(entry.pos0, 99);
        assert_eq!(&*entry.contig, "ctg1");
        assert_eq!(&*entry.genome, "genomeA");
    }

    #[test]
    fn truth_position_zero_stays_zero() {
        let path = temp_file("zero.tsv", "r1\t1\t*\t0\tgenomeA\n");
        let truth = load_truth(&path).unwrap();
        assert_eq!(
            truth[&ReadKey {
                id: "r1".into(),
                mate: 1
            }]
                .pos0,
            0
        );
    }

    #[test]
    fn short_truth_row_names_its_line() {
        let path = temp_file("short.tsv", "r1\t1\tctg1\t100\tgenomeA\nr2\t1\tctg1\n");
        let err = load_truth(&path).unwrap_err();
        assert!(err.to_string().contains(":2:"), "{err}");
    }

    #[test]
    fn a_commented_out_header_is_not_a_broken_record() {
        // Five fields, so shape alone would call it a record -- but its mate column is the word
        // "mate", and erroring on an annotation someone added is not helpful.
        let path = temp_file(
            "commented_header.tsv",
            "# read_id\tmate\tcontig\tpos\tgenome\nr1\t1\tctg1\t100\tgenomeA\n",
        );
        let truth = load_truth(&path).unwrap();
        assert_eq!(truth.len(), 1);
    }

    #[test]
    fn a_read_id_starting_with_hash_is_a_read_not_a_comment() {
        // The truth format is headerless and has no comment convention, so a leading '#' cannot be
        // assumed to mean "ignore this line" -- dropping it would lose a read without a word.
        let path = temp_file(
            "hash.tsv",
            "# a real comment\n#oddread\t1\tctg1\t100\tgenomeA\nnormal\t1\tctg1\t5\tgenomeA\n",
        );
        let truth = load_truth(&path).unwrap();

        assert_eq!(truth.len(), 2, "the comment is skipped, the read is not");
        assert!(truth.contains_key(&ReadKey {
            id: "#oddread".into(),
            mate: 1
        }));
    }

    #[test]
    fn contig2genome_ignores_comments_and_extra_columns() {
        let path = temp_file(
            "c2g.tsv",
            "# header\nctg1\tgenomeA\tignored\nctg2\tgenomeB\n",
        );
        let map = load_contig2genome(&path).unwrap();
        assert_eq!(map.len(), 2);
        assert_eq!(&*map[&Box::from("ctg1")], "genomeA");
    }
}
