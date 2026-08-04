//! Reusing the previous run's results for contenders whose inputs have not changed.
//!
//! Scoring is dominated by parsing alignment files and replaying them against a reference, and a
//! benchmark is usually re-run because *one* tool changed. Re-reading the other twenty files to
//! recompute numbers that cannot have moved is the bulk of the wall time for none of the value.
//!
//! A row is reused when its [`Fingerprint`] matches the one stored beside the outputs. The
//! fingerprint covers every input the result depends on — the files, and the options that change
//! what the numbers mean — so a row is only skipped when its answer is already known.
//!
//! `--force` ignores the cache entirely.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};

use log::{debug, warn};
use serde::{Deserialize, Serialize};

use crate::options::AlignArgs;

use super::error::{AlignError, AlignResult};
use super::meta::AlignRow;
use super::report::SampleResult;

/// Exact `f64` serialization, as bit patterns.
///
/// `serde_json` does not round-trip every `f64`: it writes `98.33333333333333` and reads back
/// `98.33333333333331`, one ULP low, because its parser is not correctly rounded. Storing the bits
/// makes a reused result bit-identical to a recomputed one, which the determinism the outputs
/// promise depends on.
pub mod f64_bits {
    use serde::{Deserialize, Deserializer, Serializer};

    /// Writes a float as its bit pattern.
    ///
    /// # Errors
    ///
    /// Propagates the serializer's own error.
    pub fn serialize<S: Serializer>(value: &f64, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_u64(value.to_bits())
    }

    /// Reads a float back from its bit pattern.
    ///
    /// # Errors
    ///
    /// Propagates the deserializer's own error.
    pub fn deserialize<'de, D: Deserializer<'de>>(deserializer: D) -> Result<f64, D::Error> {
        u64::deserialize(deserializer).map(f64::from_bits)
    }
}

/// The same, for an optional float.
pub mod opt_f64_bits {
    use serde::{Deserialize, Deserializer, Serializer};

    /// Writes an optional float as an optional bit pattern.
    ///
    /// # Errors
    ///
    /// Propagates the serializer's own error.
    pub fn serialize<S: Serializer>(value: &Option<f64>, serializer: S) -> Result<S::Ok, S::Error> {
        match value {
            Some(value) => serializer.serialize_some(&value.to_bits()),
            None => serializer.serialize_none(),
        }
    }

    /// Reads an optional float back.
    ///
    /// # Errors
    ///
    /// Propagates the deserializer's own error.
    pub fn deserialize<'de, D: Deserializer<'de>>(
        deserializer: D,
    ) -> Result<Option<f64>, D::Error> {
        Ok(Option::<u64>::deserialize(deserializer)?.map(f64::from_bits))
    }
}

/// The cache format's identity. Entries written under a different one are ignored.
///
/// It carries the crate version, because the risk this guards is not a changed struct — `serde`
/// catches those, since every cached type denies unknown fields — but a changed *computation*.
/// Three audits running have found denominator defects: fixes that alter what a number means
/// without touching a single field. A warm run after such a fix would serve the old wrong number
/// for every unchanged contender, and nothing would prompt anyone to bump a hand-written integer.
/// Tying it to the version means a release cannot serve results computed by a different one.
///
/// **Bump `FORMAT_REVISION` when a metric's computation changes within a version.**
const FORMAT_REVISION: u32 = 2;

/// The full format identity stored in the file.
fn format_identity() -> String {
    format!("{FORMAT_REVISION}+{}", env!("CARGO_PKG_VERSION"))
}

/// What a file looked like when a result was computed from it.
///
/// Size and modification time, not a content hash: the alignment files run to gigabytes, and
/// hashing them every run would cost exactly what the cache exists to save. The consequence is
/// that a file edited in place, to the same length, within the filesystem's timestamp granularity
/// looks unchanged — `--force` is the answer when that matters.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FileStamp {
    /// Path as the samplesheet gave it.
    pub path: String,
    /// Size in bytes.
    pub len: u64,
    /// Modification time in nanoseconds since the epoch, when the platform reports one. Seconds
    /// would leave a one-second window in which a same-size rewrite looks unchanged, and the finer
    /// resolution costs nothing.
    pub modified: Option<u128>,
}

impl FileStamp {
    /// Stamps a file.
    ///
    /// # Arguments
    ///
    /// * `path` - File to stamp
    ///
    /// # Returns
    ///
    /// The stamp, or `None` when the file cannot be read — an unreadable input is not a cache hit.
    pub fn of(path: &Path) -> Option<Self> {
        let meta = std::fs::metadata(path).ok()?;
        Some(Self {
            // Canonical, so two runs from different working directories do not compare a relative
            // path against an equal-looking one that means a different file.
            path: std::fs::canonicalize(path)
                .unwrap_or_else(|_| path.to_path_buf())
                .to_string_lossy()
                .into_owned(),
            len: meta.len(),
            modified: meta
                .modified()
                .ok()
                .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
                .map(|d| d.as_nanos()),
        })
    }
}

/// Everything a contender's result depends on.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Fingerprint {
    /// The alignment file.
    pub alignment: FileStamp,
    /// The truth file.
    pub truth: FileStamp,
    /// The contig map, when one is named.
    pub contig2genome: Option<FileStamp>,
    /// The reference, when one is named — it changes the replayed metrics.
    pub reference: Option<FileStamp>,
    /// The peer, and every input that determines the peer's *records* — its alignment and its
    /// reference. The head-to-head compares against those records, so this row's numbers move when
    /// either changes, though none of its own inputs did. The reference matters as much as the
    /// alignment: whether the peer has one decides whether the comparison uses replayed or
    /// self-reported edit distances.
    pub peer: Option<PeerStamp>,
    /// The options that change what the numbers mean, rendered as one string.
    pub options: String,
}

impl Fingerprint {
    /// Builds the fingerprint of one samplesheet row under the current options.
    ///
    /// # Arguments
    ///
    /// * `row` - The samplesheet row
    /// * `args` - The run's options
    /// * `peer` - The peer's row, when this row names a peer that is present on the sample
    ///
    /// # Returns
    ///
    /// The fingerprint, or `None` when an input file cannot be stamped.
    pub fn of(row: &AlignRow, args: &AlignArgs, peer: Option<&AlignRow>) -> Option<Self> {
        // Only the options that change a *result*. `--per-read`, `--outprefix` and the thread count
        // do not: the first changes what is written, the others nothing at all.
        let options = format!(
            "scoring={:?};sep={};tolerance={};verify_sample={};no_replay={};seed={};clip={}",
            row.scoring,
            row.sep,
            args.tolerance,
            args.verify_sample,
            args.no_replay,
            args.seed,
            args.clip_geometry,
        );
        Some(Self {
            alignment: FileStamp::of(&row.alignment)?,
            truth: FileStamp::of(&row.truth)?,
            // A named-but-unstampable input is a miss, not an absent input: `None` here would
            // make a deleted reference look like "no reference was ever configured".
            contig2genome: stamp_optional(row.contig2genome.as_deref())?,
            reference: stamp_optional(row.reference.as_deref())?,
            peer: match (&row.peer, peer) {
                (Some(name), Some(peer)) => Some(PeerStamp {
                    tool: name.clone(),
                    alignment: Some(FileStamp::of(&peer.alignment)?),
                    reference: stamp_optional(peer.reference.as_deref())?,
                }),
                // A named peer that is not on this sample produces no head-to-head at all, which
                // is a different result from one that does -- so the name still counts.
                (Some(name), None) => Some(PeerStamp {
                    tool: name.clone(),
                    alignment: None,
                    reference: None,
                }),
                (None, _) => None,
            },
            options,
        })
    }
}

/// What a peer contributed to a head-to-head.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PeerStamp {
    /// The peer's tool name.
    pub tool: String,
    /// Its alignment file, or `None` when the named peer is not on this sample at all — which
    /// produces no head-to-head, a different result from one that does.
    pub alignment: Option<FileStamp>,
    /// Its reference, which decides whether the comparison is replayed or self-reported.
    pub reference: Option<FileStamp>,
}

/// Stamps an optional input: `Some(None)` when it is not configured, `None` when it is configured
/// but cannot be read.
fn stamp_optional(path: Option<&Path>) -> Option<Option<FileStamp>> {
    match path {
        None => Some(None),
        Some(path) => FileStamp::of(path).map(Some),
    }
}

/// One contender's stored result.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CacheEntry {
    /// What the inputs looked like.
    pub fingerprint: Fingerprint,
    /// The scored result. `mappable` is not trusted on load — see [`Cache::reuse`].
    pub result: SampleResult,
    /// Per-genome counts, keyed by genome name.
    pub per_genome: HashMap<String, super::metrics::SubsetScore>,
}

/// Results from previous runs, keyed by `(dataset, sample, tool)`.
#[derive(Debug, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Cache {
    /// Format identity of the entries below.
    format: String,
    /// One entry per scored contender.
    entries: HashMap<String, CacheEntry>,
}

/// The cache key of a samplesheet row.
fn key(row: &AlignRow) -> String {
    format!("{}\u{1}{}\u{1}{}", row.id, row.sample, row.tool)
}

impl Cache {
    /// The path the cache lives at, beside the outputs.
    ///
    /// # Arguments
    ///
    /// * `outprefix` - The run's output prefix
    ///
    /// # Returns
    ///
    /// `<outprefix>.align_cache.json`.
    pub fn path_for(outprefix: &str) -> PathBuf {
        PathBuf::from(format!("{outprefix}.align_cache.json"))
    }

    /// Loads the cache beside `outprefix`, or an empty one.
    ///
    /// A cache that cannot be read or was written by another format version is discarded with a
    /// warning rather than failing the run: a stale cache must never be able to break a benchmark.
    ///
    /// # Arguments
    ///
    /// * `outprefix` - The run's output prefix
    ///
    /// # Returns
    ///
    /// The stored results, or an empty cache.
    pub fn load(outprefix: &str) -> Self {
        let path = Self::path_for(outprefix);
        let Ok(file) = File::open(&path) else {
            return Self::empty();
        };
        match serde_json::from_reader::<_, Self>(BufReader::new(file)) {
            Ok(cache) if cache.format == format_identity() => {
                debug!(
                    "{}: {} cached result(s)",
                    path.display(),
                    cache.entries.len()
                );
                cache
            }
            Ok(cache) => {
                warn!(
                    "{} was written by cache format {} (this is {}); recomputing everything",
                    path.display(),
                    cache.format,
                    format_identity()
                );
                Self::empty()
            }
            Err(err) => {
                warn!(
                    "{} could not be read ({err}); recomputing everything",
                    path.display()
                );
                Self::empty()
            }
        }
    }

    /// An empty cache of the current format.
    ///
    /// # Returns
    ///
    /// A cache with no entries.
    pub fn empty() -> Self {
        Self {
            format: format_identity(),
            entries: HashMap::new(),
        }
    }

    /// The stored result for a row whose inputs are unchanged.
    ///
    /// `mappable` is deliberately zeroed: it is the *field's* best result on that sample, so it
    /// moves when a contender joins or leaves the samplesheet even though nothing about this row
    /// changed. The caller recomputes it from the whole group.
    ///
    /// # Arguments
    ///
    /// * `row` - The samplesheet row
    /// * `fingerprint` - Its fingerprint under the current options
    ///
    /// # Returns
    ///
    /// The stored result and per-genome counts, or `None` on a miss.
    pub fn reuse(
        &self,
        row: &AlignRow,
        fingerprint: &Fingerprint,
    ) -> Option<(SampleResult, super::metrics::PerGenome)> {
        let entry = self.entries.get(&key(row))?;
        if &entry.fingerprint != fingerprint {
            return None;
        }
        let mut result = entry.result.clone();
        result.mappable = 0;
        let per_genome = entry
            .per_genome
            .iter()
            .map(|(genome, counts)| (std::sync::Arc::from(genome.as_str()), *counts))
            .collect();
        Some((result, per_genome))
    }

    /// Stores a freshly computed result.
    ///
    /// # Arguments
    ///
    /// * `row` - The samplesheet row it came from
    /// * `fingerprint` - Its fingerprint
    /// * `result` - The scored result
    /// * `per_genome` - Its per-genome counts
    pub fn store(
        &mut self,
        row: &AlignRow,
        fingerprint: Fingerprint,
        result: &SampleResult,
        per_genome: &super::metrics::PerGenome,
    ) {
        self.entries.insert(
            key(row),
            CacheEntry {
                fingerprint,
                result: result.clone(),
                per_genome: per_genome
                    .iter()
                    .map(|(genome, counts)| (genome.to_string(), *counts))
                    .collect(),
            },
        );
    }

    /// Writes the cache beside the outputs.
    ///
    /// # Arguments
    ///
    /// * `outprefix` - The run's output prefix
    ///
    /// # Returns
    ///
    /// `Ok(())` once written.
    ///
    /// # Errors
    ///
    /// Returns [`AlignError::Output`] when the file cannot be written.
    pub fn save(&self, outprefix: &str) -> AlignResult<()> {
        let path = Self::path_for(outprefix);
        // Written under a scratch name and renamed, as the `.fai` builder and the per-read writer
        // do: a run killed mid-write must not destroy the previous good cache.
        let mut scratch = path.as_os_str().to_os_string();
        scratch.push(format!(".tmp.{}", std::process::id()));
        let scratch = PathBuf::from(scratch);

        let file = File::create(&scratch).map_err(|e| AlignError::Output {
            path: scratch.clone(),
            message: e.to_string(),
        })?;
        serde_json::to_writer(BufWriter::new(file), self).map_err(|e| AlignError::Output {
            path: scratch.clone(),
            message: e.to_string(),
        })?;
        std::fs::rename(&scratch, &path).map_err(|e| AlignError::Output {
            path,
            message: e.to_string(),
        })
    }

    /// Drops entries for rows the samplesheet no longer names.
    ///
    /// Without this a long-lived `--outprefix` accumulates every contender ever scored under it.
    ///
    /// # Arguments
    ///
    /// * `rows` - The samplesheet's current rows
    pub fn retain(&mut self, rows: &[AlignRow]) {
        let live: std::collections::HashSet<String> = rows.iter().map(key).collect();
        self.entries.retain(|stored, _| live.contains(stored));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::meta::AlignmentFormat;
    use crate::options::ScoringMode;
    use std::io::Write;

    fn temp_dir(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("benchpro_align_cache_{name}"));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn write(dir: &Path, name: &str, content: &str) -> PathBuf {
        let path = dir.join(name);
        File::create(&path)
            .unwrap()
            .write_all(content.as_bytes())
            .unwrap();
        path
    }

    fn row(alignment: PathBuf, truth: PathBuf) -> AlignRow {
        AlignRow {
            id: "ds".into(),
            sample: "s1".into(),
            tool: "t".into(),
            alignment,
            format: AlignmentFormat::Sam,
            truth,
            reference: None,
            contig2genome: None,
            peer: None,
            scoring: ScoringMode::Species,
            sep: "_".into(),
        }
    }

    fn args() -> AlignArgs {
        AlignArgs {
            meta: PathBuf::from("m.tsv"),
            outprefix: "out".into(),
            scoring: ScoringMode::Species,
            sep: "_".into(),
            tolerance: 100,
            verify_sample: 100_000,
            no_replay: false,
            threads: 0,
            per_read: false,
            clip_geometry: false,
            seed: 0,
            validate_meta: false,
            force: false,
        }
    }

    #[test]
    fn an_unchanged_row_fingerprints_the_same() {
        let dir = temp_dir("stable");
        let r = row(write(&dir, "a.sam", "x"), write(&dir, "t.tsv", "y"));
        assert_eq!(
            Fingerprint::of(&r, &args(), None).unwrap(),
            Fingerprint::of(&r, &args(), None).unwrap()
        );
    }

    #[test]
    fn a_changed_peer_changes_the_fingerprint() {
        // The head-to-head is computed against the peer's records, so this row's numbers move when
        // the peer's file does -- reusing it would report a stale comparison as current.
        let dir = temp_dir("peer");
        let mut r = row(write(&dir, "a.sam", "x"), write(&dir, "t.tsv", "y"));
        r.peer = Some("other".into());
        let peer_row = row(write(&dir, "peer.sam", "short"), write(&dir, "t.tsv", "y"));

        let before = Fingerprint::of(&r, &args(), Some(&peer_row)).unwrap();
        write(&dir, "peer.sam", "a considerably longer file");
        let after = Fingerprint::of(&r, &args(), Some(&peer_row)).unwrap();

        assert_ne!(before, after);
    }

    #[test]
    fn the_peers_reference_is_part_of_the_fingerprint() {
        // Whether the PEER has a reference decides whether the head-to-head compares replayed or
        // self-reported edit distances, so it changes this row's numbers without touching any of
        // this row's own inputs.
        let dir = temp_dir("peer_reference");
        let mut r = row(write(&dir, "a.sam", "x"), write(&dir, "t.tsv", "y"));
        r.peer = Some("other".into());

        let mut peer_row = row(write(&dir, "peer.sam", "p"), write(&dir, "t.tsv", "y"));
        let without = Fingerprint::of(&r, &args(), Some(&peer_row)).unwrap();

        peer_row.reference = Some(write(&dir, "ref.fna", ">c\nACGT\n"));
        let with = Fingerprint::of(&r, &args(), Some(&peer_row)).unwrap();
        assert_ne!(
            without, with,
            "adding a reference to the peer must invalidate"
        );

        write(&dir, "ref.fna", ">c\nACGTACGTACGT\n");
        let changed = Fingerprint::of(&r, &args(), Some(&peer_row)).unwrap();
        assert_ne!(
            with, changed,
            "changing the peer's reference must invalidate"
        );
    }

    #[test]
    fn naming_a_peer_that_is_not_present_differs_from_naming_none() {
        let dir = temp_dir("absent_peer");
        let r = row(write(&dir, "a.sam", "x"), write(&dir, "t.tsv", "y"));
        let mut with_peer = r.clone();
        with_peer.peer = Some("ghost".into());

        assert_ne!(
            Fingerprint::of(&r, &args(), None).unwrap(),
            Fingerprint::of(&with_peer, &args(), None).unwrap()
        );
    }

    #[test]
    fn a_changed_alignment_changes_the_fingerprint() {
        let dir = temp_dir("changed");
        let alignment = write(&dir, "a.sam", "x");
        let truth = write(&dir, "t.tsv", "y");
        let before =
            Fingerprint::of(&row(alignment.clone(), truth.clone()), &args(), None).unwrap();

        write(&dir, "a.sam", "much longer content");
        let after = Fingerprint::of(&row(alignment, truth), &args(), None).unwrap();
        assert_ne!(before, after);
    }

    #[test]
    fn options_that_change_the_numbers_change_the_fingerprint() {
        let dir = temp_dir("options");
        let r = row(write(&dir, "a.sam", "x"), write(&dir, "t.tsv", "y"));
        let base = Fingerprint::of(&r, &args(), None).unwrap();

        let mut tolerance = args();
        tolerance.tolerance = 50;
        assert_ne!(base, Fingerprint::of(&r, &tolerance, None).unwrap());

        let mut seed = args();
        seed.seed = 7;
        assert_ne!(base, Fingerprint::of(&r, &seed, None).unwrap());

        let mut replay = args();
        replay.no_replay = true;
        assert_ne!(base, Fingerprint::of(&r, &replay, None).unwrap());
    }

    #[test]
    fn options_that_only_change_what_is_written_do_not() {
        // --per-read and --threads cannot move a number, so they must not invalidate a result.
        let dir = temp_dir("harmless");
        let r = row(write(&dir, "a.sam", "x"), write(&dir, "t.tsv", "y"));
        let base = Fingerprint::of(&r, &args(), None).unwrap();

        let mut per_read = args();
        per_read.per_read = true;
        assert_eq!(base, Fingerprint::of(&r, &per_read, None).unwrap());

        let mut threads = args();
        threads.threads = 4;
        assert_eq!(base, Fingerprint::of(&r, &threads, None).unwrap());
    }

    #[test]
    fn a_missing_input_has_no_fingerprint_so_cannot_hit() {
        let dir = temp_dir("missing");
        let r = row(dir.join("gone.sam"), write(&dir, "t.tsv", "y"));
        assert!(Fingerprint::of(&r, &args(), None).is_none());
    }

    #[test]
    fn a_cache_from_another_format_version_is_discarded() {
        let dir = temp_dir("version");
        let prefix = dir.join("out").to_string_lossy().into_owned();
        std::fs::write(Cache::path_for(&prefix), r#"{"format":999,"entries":{}}"#).unwrap();
        // Loads as empty rather than failing the run.
        assert!(Cache::load(&prefix).entries.is_empty());
    }

    #[test]
    fn an_unreadable_cache_is_discarded_rather_than_fatal() {
        let dir = temp_dir("corrupt");
        let prefix = dir.join("out").to_string_lossy().into_owned();
        std::fs::write(Cache::path_for(&prefix), "not json at all").unwrap();
        assert!(Cache::load(&prefix).entries.is_empty());
    }
}
