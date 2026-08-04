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

use super::error::{AlignError, AlignResult};

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
        if line.is_empty() || line[0] == b'#' {
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
            other => {
                return Err(AlignError::parse(
                    path,
                    i + 1,
                    format!("mate '{}' is not 1 or 2", String::from_utf8_lossy(other)),
                ))
            }
        };
        let Some(pos) = parse_u64(pos) else {
            return Err(AlignError::parse(
                path,
                i + 1,
                format!(
                    "position '{}' is not a number",
                    String::from_utf8_lossy(pos)
                ),
            ));
        };

        truth.insert(
            ReadKey {
                id: String::from_utf8_lossy(id).into_owned().into(),
                mate,
            },
            TruthEntry {
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
        if line.is_empty() || line.starts_with('#') {
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
