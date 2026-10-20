//! Taxonomy inference for meta rows that omit the `Taxonomy` column, and consistency checks
//! for rows that provide it.
//!
//! Two layers of evidence are combined:
//!
//! * **Text evidence** (always available): features that only occur in GTDB names (`_A`
//!   suffixes, placeholder names such as `UBA1234`, GTDB accessions) or only in NCBI names
//!   (numeric TAXPATHs, `[Clostridium]` brackets, `sp.`, `Candidatus`).
//! * **Index evidence** (optional): the `taxdet` presence indexes from
//!   `4less/ncbi-taxdump-analyzer`, which list every name in every GTDB release and NCBI
//!   taxdump. They confirm the taxonomy and identify the release.
//!
//! Every decision keeps its evidence as text so it can be reported to the user.

use std::{
    collections::{HashMap, HashSet},
    fmt,
    fs::File,
    io::{self, BufRead, BufReader},
    path::{Path, PathBuf},
};

use lazy_static::lazy_static;
use log::{info, warn};
use polars::prelude::{NamedFrom, Series};
use regex::Regex;

use crate::meta::{Meta, MetaColumn, MetaColumnStrings};

/// Maximum number of data lines scanned per profile.
const MAX_SCAN_LINES: usize = 200_000;
/// Minimum number of distinct names needed before index evidence counts.
const MIN_INDEX_NAMES: usize = 5;

lazy_static! {
    static ref GTDB_SUFFIX: Regex =
        Regex::new(r"\b[pcofgs]__[A-Z][A-Za-z0-9-]*_[A-Z]{1,2}\b").unwrap();
    static ref GTDB_PLACEHOLDER: Regex =
        Regex::new(r"\b[gs]__[A-Z]{2,}-?\d{2,}\b| sp\d{9}\b").unwrap();
    static ref GTDB_ACCESSION: Regex = Regex::new(r"\b(RS|GB)_GC[AF]_\d{9}").unwrap();
    static ref GTDB_DOMAIN: Regex = Regex::new(r"(^|[\t|;])d__(Bacteria|Archaea)\b").unwrap();
    static ref NCBI_TAXPATH: Regex = Regex::new(r"(^|\t)\d+(\|\d*)+(\t|$)").unwrap();
    static ref NCBI_BRACKETS: Regex = Regex::new(r"\[[A-Z][a-z]+\]").unwrap();
    static ref NCBI_NAMING: Regex = Regex::new(r" sp\. |\buncultured |\bCandidatus ").unwrap();
    static ref NCBI_KINGDOM: Regex = Regex::new(r"(^|[\t|;])(k|sk)__(Bacteria|Archaea)\b").unwrap();
    static ref RANK_PREFIXED: Regex = Regex::new(r"\b[a-z]{1,2}__([^|;\t]+)").unwrap();
}

/// Taxonomy family that determines how names are parsed and compared.
///
/// # Examples
///
/// ```ignore
/// assert_eq!(TaxonomyFamily::from_label("GTDB r214"), Some(TaxonomyFamily::GTDB));
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TaxonomyFamily {
    GTDB,
    NCBI,
}

impl TaxonomyFamily {
    /// Maps a meta `Taxonomy` value to its family; `None` for custom taxonomies.
    ///
    /// # Arguments
    ///
    /// * `label` - Value from the meta `Taxonomy` column
    ///
    /// # Returns
    ///
    /// Family, or `None` when the label is not GTDB/NCBI-based.
    pub fn from_label(label: &str) -> Option<Self> {
        if label.starts_with("GTDB") {
            Some(Self::GTDB)
        } else if label.starts_with("NCBI") || label.starts_with("ChocoPhlAn") {
            Some(Self::NCBI)
        } else {
            None
        }
    }
}

impl fmt::Display for TaxonomyFamily {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::GTDB => write!(f, "GTDB"),
            Self::NCBI => write!(f, "NCBI"),
        }
    }
}

/// How strongly the evidence supports an inferred taxonomy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Confidence {
    None,
    Low,
    Medium,
    High,
}

impl fmt::Display for Confidence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            Self::None => "none",
            Self::Low => "low",
            Self::Medium => "medium",
            Self::High => "high",
        };
        write!(f, "{}", s)
    }
}

/// Feature counts and candidate names collected from one profile file.
#[derive(Debug, Clone, Default)]
pub struct TextEvidence {
    /// `(feature, number of lines)` pairs that indicate GTDB.
    pub gtdb: Vec<(&'static str, usize)>,
    /// `(feature, number of lines)` pairs that indicate NCBI.
    pub ncbi: Vec<(&'static str, usize)>,
    /// Distinct taxon names without rank prefixes, used for index lookups.
    pub names: HashSet<String>,
}

impl TextEvidence {
    fn total(features: &[(&'static str, usize)]) -> usize {
        features.iter().map(|(_, n)| n).sum()
    }

    fn describe(features: &[(&'static str, usize)]) -> String {
        let parts = features
            .iter()
            .filter(|(_, n)| *n > 0)
            .map(|(name, n)| format!("{} {}", n, name))
            .collect::<Vec<_>>();
        if parts.is_empty() {
            "none".to_string()
        } else {
            parts.join(", ")
        }
    }
}

/// Result of looking up one profile's names in one taxonomy index.
#[derive(Debug, Clone)]
pub struct IndexMatch {
    pub family: TaxonomyFamily,
    pub found: usize,
    pub total: usize,
    /// Versions containing the most names (ties included), in index order.
    pub best_versions: Vec<String>,
    pub best_count: usize,
    /// All versions of this index and how many of the names each contains.
    pub versions: Vec<String>,
    pub per_version: Vec<usize>,
}

impl IndexMatch {
    fn fraction(&self) -> f64 {
        if self.total == 0 {
            0.0
        } else {
            self.found as f64 / self.total as f64
        }
    }

    fn describe(&self) -> String {
        let versions = describe_versions(&self.best_versions);
        format!(
            "{} index: {}/{} names ({:.0}%), best {} ({} names)",
            self.family,
            self.found,
            self.total,
            100.0 * self.fraction(),
            versions,
            self.best_count
        )
    }
}

/// Inferred taxonomy for one profile file, with the evidence behind it.
#[derive(Debug, Clone)]
pub struct FileInference {
    pub family: Option<TaxonomyFamily>,
    pub confidence: Confidence,
    /// Human-readable summary of all evidence.
    pub evidence: String,
    /// Release label such as `GTDB r214`, when the index identifies a single best release.
    pub release: Option<String>,
    /// All best-matching versions from the index of the inferred family.
    pub release_candidates: Vec<String>,
    /// Index match for the inferred family, if an index was used.
    pub index_match: Option<IndexMatch>,
}

/// Name-to-bitset presence data for one taxonomy, restricted to the queried names.
struct FamilyIndex {
    family: TaxonomyFamily,
    versions: Vec<String>,
    names: HashMap<String, Vec<u64>>,
}

/// Presence indexes produced by `taxdet` (`gtdb_index/` and `ncbi_index/` subfolders).
///
/// # Examples
///
/// ```ignore
/// let names: HashSet<String> = ["Akkermansia muciniphila".to_string()].into();
/// let index = TaxonomyIndex::load(Path::new("~/.taxdet/index"), &names)?;
/// ```
pub struct TaxonomyIndex {
    pub root: PathBuf,
    families: Vec<FamilyIndex>,
}

impl TaxonomyIndex {
    /// Loads the indexes below `root`, keeping only rows for `names`.
    ///
    /// # Arguments
    ///
    /// * `root` - Folder containing `gtdb_index/` and/or `ncbi_index/`
    /// * `names` - Names to keep (without rank prefixes)
    ///
    /// # Returns
    ///
    /// Loaded index; families whose folder is missing are skipped.
    ///
    /// # Errors
    ///
    /// Returns an I/O error when an index folder exists but cannot be read, or when no
    /// family index is found below `root`.
    pub fn load(root: &Path, names: &HashSet<String>) -> io::Result<Self> {
        let mut families = Vec::new();
        for (family, folder) in [
            (TaxonomyFamily::GTDB, "gtdb_index"),
            (TaxonomyFamily::NCBI, "ncbi_index"),
        ] {
            let dir = root.join(folder);
            if dir.join("version_columns.tsv").is_file() {
                families.push(Self::load_family(family, &dir, names)?);
            }
        }
        if families.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!(
                    "no gtdb_index/ or ncbi_index/ folder in '{}'",
                    root.display()
                ),
            ));
        }
        Ok(Self {
            root: root.to_path_buf(),
            families,
        })
    }

    fn load_family(
        family: TaxonomyFamily,
        dir: &Path,
        names: &HashSet<String>,
    ) -> io::Result<FamilyIndex> {
        let mut versions = Vec::new();
        for line in BufReader::new(File::open(dir.join("version_columns.tsv"))?)
            .lines()
            .skip(1)
        {
            let line = line?;
            let mut fields = line.split('\t');
            if let (Some(_), Some(version)) = (fields.next(), fields.next()) {
                versions.push(version.trim().to_string());
            }
        }
        let chunks = versions.len().div_ceil(64);

        let mut found = HashMap::new();
        for line in BufReader::new(File::open(dir.join("scientific_name_matrix.tsv"))?)
            .lines()
            .skip(1)
        {
            let line = line?;
            let Some((name, bits)) = line.rsplit_once('\t') else {
                continue;
            };
            if names.contains(name) {
                found.insert(name.to_string(), parse_bitset_hex(bits, chunks));
            }
        }

        Ok(FamilyIndex {
            family,
            versions,
            names: found,
        })
    }

    /// Families available in this index.
    ///
    /// # Returns
    ///
    /// Loaded taxonomy families.
    pub fn families(&self) -> Vec<TaxonomyFamily> {
        self.families.iter().map(|f| f.family).collect()
    }

    fn lookup(&self, names: &HashSet<String>) -> Vec<IndexMatch> {
        self.families
            .iter()
            .map(|index| {
                let mut per_version = vec![0usize; index.versions.len()];
                let mut found = 0;
                for name in names {
                    if let Some(bits) = index.names.get(name) {
                        found += 1;
                        for (v, count) in per_version.iter_mut().enumerate() {
                            if (bits[v / 64] >> (v % 64)) & 1 == 1 {
                                *count += 1;
                            }
                        }
                    }
                }
                let best_count = per_version.iter().copied().max().unwrap_or(0);
                let best_versions = if best_count == 0 {
                    Vec::new()
                } else {
                    index
                        .versions
                        .iter()
                        .zip(per_version.iter())
                        .filter(|(_, &c)| c == best_count)
                        .map(|(v, _)| v.clone())
                        .collect()
                };
                IndexMatch {
                    family: index.family,
                    found,
                    total: names.len(),
                    best_versions,
                    best_count,
                    versions: index.versions.clone(),
                    per_version,
                }
            })
            .collect()
    }
}

fn parse_bitset_hex(s: &str, chunks: usize) -> Vec<u64> {
    (0..chunks)
        .map(|i| {
            s.get(i * 16..i * 16 + 16)
                .and_then(|chunk| u64::from_str_radix(chunk, 16).ok())
                .unwrap_or(0)
        })
        .collect()
}

/// Converts an index version id to a meta `Taxonomy` label.
///
/// # Arguments
///
/// * `version` - Version id such as `g-r214` or `new_taxdump_2024-03-01`
///
/// # Returns
///
/// Label such as `GTDB r214` or `NCBI 2024-03-01`.
pub fn release_label(version: &str) -> String {
    if let Some(release) = version.strip_prefix("g-").or(version.strip_prefix("gtdb-")) {
        return format!("GTDB {}", release);
    }
    let date = version.rsplit('_').next().unwrap_or(version);
    format!("NCBI {}", date)
}

/// Returns the file to scan for a meta profile reference (`matrix:<path>::<column>` aware).
fn scan_path(reference: &Path) -> PathBuf {
    let s = reference.to_string_lossy();
    match s.strip_prefix("matrix:") {
        Some(rest) => PathBuf::from(rest.split("::").next().unwrap_or(rest)),
        None => reference.to_path_buf(),
    }
}

/// Collects text evidence from a profile file.
///
/// # Arguments
///
/// * `path` - Profile path (or `matrix:<path>::<column>` reference)
///
/// # Returns
///
/// Feature counts and candidate names.
///
/// # Errors
///
/// Returns an I/O error when the file cannot be read.
///
/// # Examples
///
/// ```ignore
/// let evidence = scan_file(Path::new("profile.tsv"))?;
/// ```
pub fn scan_file(path: &Path) -> io::Result<TextEvidence> {
    let reader = BufReader::new(File::open(scan_path(path))?);
    let mut counts = [0usize; 8];
    let mut names = HashSet::new();

    for line in reader
        .lines()
        .map_while(Result::ok)
        .filter(|l| !l.starts_with('#') && !l.starts_with('@') && !l.trim().is_empty())
        .take(MAX_SCAN_LINES)
    {
        let checks = [
            &*GTDB_SUFFIX,
            &*GTDB_PLACEHOLDER,
            &*GTDB_ACCESSION,
            &*GTDB_DOMAIN,
            &*NCBI_TAXPATH,
            &*NCBI_BRACKETS,
            &*NCBI_NAMING,
            &*NCBI_KINGDOM,
        ];
        for (count, regex) in counts.iter_mut().zip(checks) {
            if regex.is_match(&line) {
                *count += 1;
            }
        }
        collect_names(&line, &mut names);
    }

    Ok(TextEvidence {
        gtdb: vec![
            ("lines with _A-style suffixes", counts[0]),
            ("lines with placeholder names", counts[1]),
            ("lines with GTDB accessions", counts[2]),
            ("lines with d__ domain", counts[3]),
        ],
        ncbi: vec![
            ("lines with numeric TAXPATH", counts[4]),
            ("lines with [Bracketed] names", counts[5]),
            ("lines with sp./uncultured/Candidatus", counts[6]),
            ("lines with k__ kingdom", counts[7]),
        ],
        names,
    })
}

/// Adds the taxon names found on one line (rank-prefixed tokens, or `|`-separated name paths).
fn collect_names(line: &str, names: &mut HashSet<String>) {
    let mut push = |name: &str| {
        let name = name.trim();
        if !name.is_empty()
            && !name.contains('(')
            && name.chars().next().is_some_and(|c| c.is_alphabetic())
        {
            names.insert(name.to_string());
        }
    };
    if line.contains("__") {
        for capture in RANK_PREFIXED.captures_iter(line) {
            push(&capture[1]);
        }
        return;
    }
    for field in line.split('\t') {
        if field.contains('|') && field.chars().any(|c| c.is_alphabetic()) {
            field.split('|').for_each(&mut push);
        }
    }
}

/// Decides a file's taxonomy from text and (optional) index evidence.
///
/// # Arguments
///
/// * `text` - Text evidence from [`scan_file`]
/// * `index` - Optional presence index
///
/// # Returns
///
/// Inferred family, confidence and an evidence summary.
pub fn infer(text: &TextEvidence, index: Option<&TaxonomyIndex>) -> FileInference {
    let (text_family, text_confidence) = decide_from_text(text);
    let mut evidence = format!(
        "text: GTDB features [{}]; NCBI features [{}]",
        TextEvidence::describe(&text.gtdb),
        TextEvidence::describe(&text.ncbi)
    );

    let matches = index.map(|i| i.lookup(&text.names)).unwrap_or_default();
    for m in &matches {
        evidence.push_str("; ");
        evidence.push_str(&m.describe());
    }
    let index_family = decide_from_index(&matches);

    let (family, confidence) = match (text_family, index_family) {
        (Some(t), Some(i)) if t == i => (Some(t), Confidence::High),
        (Some(_), Some(i)) => (Some(i), Confidence::Medium),
        (None, Some(i)) => (Some(i), Confidence::High),
        (t, None) => (t, text_confidence),
    };

    let index_match = family.and_then(|f| matches.iter().find(|m| m.family == f).cloned());
    let candidates = index_match
        .as_ref()
        .map(|m| m.best_versions.clone())
        .unwrap_or_default();
    let release = match candidates.as_slice() {
        [one] => Some(release_label(one)),
        _ => None,
    };

    FileInference {
        family,
        confidence,
        evidence,
        release,
        release_candidates: candidates,
        index_match,
    }
}

/// Counts prediction names that the gold standard's best release does not contain.
///
/// # Arguments
///
/// * `prediction` - Index match of the prediction
/// * `gold_std` - Index match of the gold standard (same family)
///
/// # Returns
///
/// `(missing, gold_release)`: names present in the prediction's best release but absent
/// from the gold standard's best release, and that release's version id.
///
/// # Examples
///
/// ```ignore
/// let (missing, release) = release_gap(&pred_match, &gold_match)?;
/// ```
pub fn release_gap(prediction: &IndexMatch, gold_std: &IndexMatch) -> Option<(usize, String)> {
    if prediction.family != gold_std.family {
        return None;
    }
    let (count, version) = gold_std
        .best_versions
        .iter()
        .filter_map(|v| {
            let i = prediction.versions.iter().position(|x| x == v)?;
            Some((prediction.per_version[i], v.clone()))
        })
        .max_by_key(|(count, _)| *count)?;
    Some((prediction.best_count.saturating_sub(count), version))
}

fn decide_from_text(text: &TextEvidence) -> (Option<TaxonomyFamily>, Confidence) {
    let gtdb = TextEvidence::total(&text.gtdb);
    let ncbi = TextEvidence::total(&text.ncbi);
    let pick = |major: usize, minor: usize, family| {
        if minor == 0 {
            let confidence = if major >= 3 {
                Confidence::High
            } else {
                Confidence::Low
            };
            (Some(family), confidence)
        } else if major >= 10 * minor && major >= 5 {
            (Some(family), Confidence::Medium)
        } else {
            (Some(family), Confidence::Low)
        }
    };
    match gtdb.cmp(&ncbi) {
        _ if gtdb == 0 && ncbi == 0 => (None, Confidence::None),
        std::cmp::Ordering::Greater => pick(gtdb, ncbi, TaxonomyFamily::GTDB),
        std::cmp::Ordering::Less => pick(ncbi, gtdb, TaxonomyFamily::NCBI),
        std::cmp::Ordering::Equal => (None, Confidence::Low),
    }
}

fn decide_from_index(matches: &[IndexMatch]) -> Option<TaxonomyFamily> {
    let best = matches
        .iter()
        .filter(|m| m.total >= MIN_INDEX_NAMES)
        .max_by(|a, b| a.fraction().total_cmp(&b.fraction()))?;
    if best.fraction() < 0.8 {
        return None;
    }
    // With both indexes loaded, the winner must clearly beat the other family.
    let runner_up = matches
        .iter()
        .filter(|m| m.family != best.family)
        .map(IndexMatch::fraction)
        .fold(0.0, f64::max);
    let clear = matches.len() < 2 || best.fraction() - runner_up >= 0.1;
    clear.then_some(best.family)
}

/// Options controlling taxonomy inference.
#[derive(Debug, Clone, Default)]
pub struct InferenceOptions {
    /// Folder with `taxdet` indexes, if any.
    pub index_dir: Option<PathBuf>,
    /// Continue when prediction and gold standard appear to use different taxonomies.
    pub allow_taxonomy_mismatch: bool,
}

/// Taxonomy decision for one meta row, written to `<outprefix>_inputs.tsv`.
#[derive(Debug, Clone, Default)]
pub struct RowTaxonomy {
    pub id: String,
    /// Final `Taxonomy` value used for parsing.
    pub taxonomy: String,
    /// `given` (from the meta file) or `inferred`.
    pub source: String,
    pub profile: Option<FileInference>,
    pub goldstd: Option<FileInference>,
    /// Warnings raised for this row.
    pub warnings: Vec<RowWarning>,
}

/// One warning about a meta row, with a kind for grouping in logs.
#[derive(Debug, Clone)]
pub struct RowWarning {
    pub kind: &'static str,
    pub message: String,
}

/// Prediction and gold standard look like different taxonomies.
pub const WARN_MISMATCH: &str = "prediction and gold standard use different taxonomies";
/// The meta `Taxonomy` value contradicts the evidence in a profile.
pub const WARN_GIVEN_CONFLICT: &str = "Taxonomy column contradicts the profile contents";
/// Taxonomy was inferred with medium confidence only.
pub const WARN_MEDIUM: &str = "taxonomy inferred with medium confidence";
/// Prediction names are missing from the gold standard's release.
pub const WARN_RELEASE_GAP: &str = "prediction names missing from the gold standard's release";

/// Resolves taxonomy for every meta row: infers missing values and checks given ones.
///
/// Writes the final values into the `Taxonomy` and `TaxonomySource` columns of
/// `meta.raw` and into each entry.
///
/// # Arguments
///
/// * `meta` - Parsed meta table; updated in place
/// * `options` - Index location and mismatch handling
///
/// # Returns
///
/// One [`RowTaxonomy`] per meta entry.
///
/// # Errors
///
/// Returns a message listing every row whose taxonomy cannot be determined, or whose
/// prediction and gold standard use different taxonomies.
///
/// # Examples
///
/// ```ignore
/// let rows = resolve_meta(&mut meta, &InferenceOptions::default())?;
/// ```
pub fn resolve_meta(
    meta: &mut Meta,
    options: &InferenceOptions,
) -> Result<Vec<RowTaxonomy>, String> {
    let mut evidence: HashMap<PathBuf, Option<TextEvidence>> = HashMap::new();
    for entry in &meta.entries {
        for path in [&entry.profile, &entry.goldstd] {
            evidence
                .entry(path.clone())
                .or_insert_with(|| scan_file(path).ok());
        }
    }

    let index = options.index_dir.as_ref().and_then(|dir| {
        let names = evidence
            .values()
            .flatten()
            .flat_map(|e| e.names.iter().cloned())
            .collect::<HashSet<_>>();
        match TaxonomyIndex::load(dir, &names) {
            Ok(index) => {
                let families = index
                    .families()
                    .iter()
                    .map(|f| f.to_string())
                    .collect::<Vec<_>>();
                info!(
                    "Taxonomy index: {} ({})",
                    index.root.display(),
                    families.join(", ")
                );
                Some(index)
            }
            Err(e) => {
                warn!("Taxonomy index not used: {}", e);
                None
            }
        }
    });

    let infer_path = |path: &PathBuf| {
        evidence
            .get(path)
            .and_then(|e| e.as_ref())
            .map(|e| infer(e, index.as_ref()))
    };

    let mut rows = Vec::with_capacity(meta.entries.len());
    let mut errors = Vec::new();
    for entry in &meta.entries {
        let profile = infer_path(&entry.profile);
        let goldstd = infer_path(&entry.goldstd);
        let mut row = RowTaxonomy {
            id: entry.id.clone(),
            ..Default::default()
        };

        let confident = |f: &Option<FileInference>| {
            f.as_ref()
                .filter(|f| f.confidence >= Confidence::Medium)
                .and_then(|f| f.family)
        };

        if let (Some(p), Some(g)) = (confident(&profile), confident(&goldstd)) {
            if p != g {
                let message = format!(
                    "row '{}': prediction looks like {} but gold standard looks like {}; \
                     taxon names cannot match. Prediction: {}. Gold standard: {}.",
                    entry.id,
                    p,
                    g,
                    profile.as_ref().unwrap().evidence,
                    goldstd.as_ref().unwrap().evidence
                );
                if options.allow_taxonomy_mismatch {
                    row.warnings.push(RowWarning {
                        kind: WARN_MISMATCH,
                        message,
                    });
                } else {
                    errors.push(message);
                }
            }
        }

        if !entry.taxonomy.is_empty() {
            row.taxonomy = entry.taxonomy.clone();
            row.source = "given".to_string();
            let given = TaxonomyFamily::from_label(&entry.taxonomy);
            for (what, inference) in [("prediction", &profile), ("gold standard", &goldstd)] {
                if let (Some(given), Some(inferred)) = (given, confident(inference)) {
                    if given != inferred {
                        row.warnings.push(RowWarning {
                            kind: WARN_GIVEN_CONFLICT,
                            message: format!(
                                "row '{}': Taxonomy is '{}' but the {} looks like {} ({}). \
                             Using '{}'.",
                                entry.id,
                                entry.taxonomy,
                                what,
                                inferred,
                                inference.as_ref().unwrap().evidence,
                                entry.taxonomy
                            ),
                        });
                    }
                }
            }
        } else {
            row.source = "inferred".to_string();
            match (&profile, &goldstd) {
                (Some(p), Some(g))
                    if p.family.is_some()
                        && p.confidence >= Confidence::Medium
                        && g.family.is_some()
                        && g.confidence >= Confidence::Medium =>
                {
                    let family = p.family.unwrap();
                    row.taxonomy = match (&p.release, &g.release) {
                        (Some(a), Some(b)) if a == b => a.clone(),
                        _ => family.to_string(),
                    };
                    for (what, f) in [("prediction", p), ("gold standard", g)] {
                        if f.confidence == Confidence::Medium {
                            row.warnings.push(RowWarning {
                                kind: WARN_MEDIUM,
                                message: format!(
                                    "row '{}': {} taxonomy inferred as {} with medium confidence \
                                 ({}). Add a Taxonomy column to confirm.",
                                    entry.id, what, family, f.evidence
                                ),
                            });
                        }
                    }
                }
                _ => {
                    let describe = |what: &str, path: &Path, f: &Option<FileInference>| match f {
                        None => format!("{} '{}' could not be read", what, path.display()),
                        Some(f) => format!(
                            "{} '{}': {} ({} confidence; {})",
                            what,
                            path.display(),
                            f.family.map_or("unknown".to_string(), |x| x.to_string()),
                            f.confidence,
                            f.evidence
                        ),
                    };
                    errors.push(format!(
                        "row '{}': cannot tell the taxonomy. {}. {}. Add a Taxonomy column \
                         (GTDB or NCBI) for this row.",
                        entry.id,
                        describe("prediction", &entry.profile, &profile),
                        describe("gold standard", &entry.goldstd, &goldstd)
                    ));
                }
            }
        }

        let gap = match (&profile, &goldstd) {
            (Some(p), Some(g)) => match (&p.index_match, &g.index_match) {
                (Some(pm), Some(gm)) => release_gap(pm, gm).map(|gap| (gap, pm.best_count)),
                _ => None,
            },
            _ => None,
        };
        if let Some(((missing, gold_release), found)) = gap {
            // Small gaps are normal between neighbouring NCBI taxdumps.
            if missing >= 2 && missing * 100 >= found {
                row.warnings.push(RowWarning {
                    kind: WARN_RELEASE_GAP,
                    message: format!(
                    "row '{}': prediction best matches {}, but {} of its names ({:.1}%) do not \
                     exist in the gold standard's release {}; those taxa can only count as \
                     FP/FN.",
                    entry.id,
                    describe_versions(&profile.as_ref().unwrap().release_candidates),
                    missing,
                    100.0 * missing as f64 / found.max(1) as f64,
                    release_label(&gold_release)
                ),
                });
            }
        }

        row.profile = profile;
        row.goldstd = goldstd;
        rows.push(row);
    }

    if !errors.is_empty() {
        return Err(errors.join("\n"));
    }

    apply_to_meta(meta, &rows);
    Ok(rows)
}

/// Describes a set of tied versions, e.g. `GTDB r214` or `NCBI 2018-09-01 .. 2021-03-01 (57)`.
///
/// # Arguments
///
/// * `versions` - Version ids in index order
///
/// # Returns
///
/// Short human-readable description.
pub fn describe_versions(versions: &[String]) -> String {
    match versions {
        [] => "no release".to_string(),
        [one] => release_label(one),
        [first, .., last] => format!(
            "{} .. {} ({} releases)",
            release_label(first),
            release_label(last),
            versions.len()
        ),
    }
}

fn apply_to_meta(meta: &mut Meta, rows: &[RowTaxonomy]) {
    let by_id = rows
        .iter()
        .map(|r| (r.id.as_str(), r))
        .collect::<HashMap<_, _>>();
    for entry in meta.entries.iter_mut() {
        if let Some(row) = by_id.get(entry.id.as_str()) {
            entry.taxonomy = row.taxonomy.clone();
        }
    }

    let Ok(ids) = meta
        .raw
        .column(MetaColumn::ID.to_str())
        .and_then(|c| c.str().cloned())
    else {
        return;
    };
    let lookup = |id: Option<&str>, f: fn(&RowTaxonomy) -> String| {
        id.and_then(|id| by_id.get(id))
            .map(|r| f(r))
            .unwrap_or_default()
    };
    let taxonomy = ids
        .iter()
        .map(|id| lookup(id, |r| r.taxonomy.clone()))
        .collect::<Vec<_>>();
    let source = ids
        .iter()
        .map(|id| lookup(id, |r| r.source.clone()))
        .collect::<Vec<_>>();
    let _ = meta
        .raw
        .with_column(Series::new(MetaColumnStrings::TAXONOMY.into(), taxonomy));
    let _ = meta
        .raw
        .with_column(Series::new("TaxonomySource".into(), source));
}

#[cfg(test)]
mod tests {
    use super::*;

    fn evidence_from(lines: &[&str]) -> TextEvidence {
        static COUNTER: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        let n = COUNTER.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
        let dir = std::env::temp_dir().join(format!("benchpro_taxinf_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(format!("{}.tsv", n));
        std::fs::write(&path, lines.join("\n")).unwrap();
        let evidence = scan_file(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        evidence
    }

    #[test]
    fn test_gtdb_profile_is_detected_with_high_confidence() {
        let text = evidence_from(&[
            "d__Bacteria;p__Firmicutes_A;c__Clostridia;g__Clostridium_K;s__Clostridium_K a\t0.1",
            "d__Bacteria;p__Bacteroidota;g__CAG-110;s__CAG-110 sp000432435\t0.2",
            "d__Bacteria;p__Bacillota_A;g__Ruminococcus_D;s__Ruminococcus_D albus\t0.3",
        ]);
        let result = infer(&text, None);
        assert_eq!(result.family, Some(TaxonomyFamily::GTDB));
        assert_eq!(result.confidence, Confidence::High);
        assert!(text.names.contains("Clostridium_K a"));
    }

    #[test]
    fn test_ncbi_cami_profile_is_detected_with_high_confidence() {
        let text = evidence_from(&[
            "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE",
            "2\tsuperkingdom\t2\tBacteria\t100.0",
            "1239\tphylum\t2|1239\tBacteria|Firmicutes\t49.7",
            "31979\tfamily\t2|1239|186801|186802|31979\tBacteria|Firmicutes|Clostridia|\
             Clostridiales|Clostridiaceae\t21.2",
        ]);
        let result = infer(&text, None);
        assert_eq!(result.family, Some(TaxonomyFamily::NCBI));
        assert_eq!(result.confidence, Confidence::Low);
        assert!(text.names.contains("Clostridiaceae"));
    }

    #[test]
    fn test_prefixed_ncbi_profile_is_not_mistaken_for_gtdb() {
        let text = evidence_from(&[
            "k__Bacteria|p__Firmicutes|c__Clostridia|g__[Clostridium]|s__[Clostridium] scindens",
            "k__Bacteria|p__Bacteroidetes|g__Prevotella|s__Prevotella sp. CAG:279",
            "k__Bacteria|p__Firmicutes|g__Faecalibacterium|s__Faecalibacterium prausnitzii",
        ]);
        let result = infer(&text, None);
        assert_eq!(result.family, Some(TaxonomyFamily::NCBI));
        assert_eq!(result.confidence, Confidence::High);
    }

    #[test]
    fn test_no_evidence_gives_no_family() {
        let text = evidence_from(&["Bacteria\t1.0", "Firmicutes\t0.5"]);
        let result = infer(&text, None);
        assert_eq!(result.family, None);
        assert_eq!(result.confidence, Confidence::None);
    }

    #[test]
    fn test_release_label() {
        assert_eq!(release_label("g-r214"), "GTDB r214");
        assert_eq!(release_label("new_taxdump_2024-03-01"), "NCBI 2024-03-01");
        assert_eq!(
            TaxonomyFamily::from_label("GTDB r207"),
            Some(TaxonomyFamily::GTDB)
        );
        assert_eq!(
            TaxonomyFamily::from_label("ChocoPhlAn"),
            Some(TaxonomyFamily::NCBI)
        );
        assert_eq!(TaxonomyFamily::from_label("Custom"), None);
    }

    #[test]
    fn test_parse_bitset_hex_bit_order() {
        // Version 0 is the lowest bit of the first 16-hex-digit chunk.
        let bits = parse_bitset_hex("0000000000000005", 1);
        assert_eq!(bits, vec![0b101]);
    }

    const GTDB_LINES: &[&str] = &[
        "d__Bacteria;p__Firmicutes_A;c__Clostridia;g__Clostridium_K;s__Clostridium_K a\t0.4",
        "d__Bacteria;p__Bacteroidota;g__CAG-110;s__CAG-110 sp000432435\t0.3",
        "d__Bacteria;p__Bacillota_A;g__Ruminococcus_D;s__Ruminococcus_D albus\t0.3",
    ];
    const NCBI_LINES: &[&str] = &[
        "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE",
        "1239\tphylum\t2|1239\tBacteria|Firmicutes\t50",
        "186801\tclass\t2|1239|186801\tBacteria|Firmicutes|Clostridia\t50",
        "39491\tspecies\t2|1239|186801|39491\tBacteria|Firmicutes|Clostridia|[Eubacterium] rectale\t50",
        "853\tspecies\t2|1239|186801|853\tBacteria|Firmicutes|Clostridia|Faecalibacterium sp. A\t50",
    ];

    fn write_profile(name: &str, lines: &[&str]) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("benchpro_resolve_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(name);
        std::fs::write(&path, lines.join("\n")).unwrap();
        path
    }

    fn meta_for(rows: &[(&str, Option<&str>, &Path, &Path)]) -> Meta {
        let column = |f: &dyn Fn(&(&str, Option<&str>, &Path, &Path)) -> String| {
            rows.iter().map(f).collect::<Vec<_>>()
        };
        let mut columns = vec![
            Series::new("ID".into(), column(&|r| r.0.to_string())),
            Series::new("Sample".into(), column(&|r| r.0.to_string())),
            Series::new("Tool".into(), column(&|_| "tool".to_string())),
            Series::new("Profile".into(), column(&|r| r.2.display().to_string())),
            Series::new("GoldStd".into(), column(&|r| r.3.display().to_string())),
        ];
        if rows.iter().any(|r| r.1.is_some()) {
            columns.push(Series::new(
                "Taxonomy".into(),
                column(&|r| r.1.unwrap_or("").to_string()),
            ));
        }
        Meta::from_polars_df(polars::frame::DataFrame::new(columns).unwrap()).unwrap()
    }

    #[test]
    fn test_resolve_meta_infers_missing_taxonomy_and_records_source() {
        let gtdb_a = write_profile("infer_a.tsv", GTDB_LINES);
        let gtdb_b = write_profile("infer_b.tsv", GTDB_LINES);
        let mut meta = meta_for(&[("r1", None, &gtdb_a, &gtdb_b)]);

        let rows = resolve_meta(&mut meta, &InferenceOptions::default()).unwrap();

        assert_eq!(rows[0].taxonomy, "GTDB");
        assert_eq!(rows[0].source, "inferred");
        assert_eq!(meta.entries[0].taxonomy, "GTDB");
        let source = meta
            .raw
            .column("TaxonomySource")
            .unwrap()
            .str()
            .unwrap()
            .get(0);
        assert_eq!(source, Some("inferred"));
    }

    #[test]
    fn test_resolve_meta_rejects_prediction_gold_mismatch_unless_allowed() {
        let gtdb = write_profile("mismatch_gtdb.tsv", GTDB_LINES);
        let ncbi = write_profile("mismatch_ncbi.tsv", NCBI_LINES);

        let mut meta = meta_for(&[("r1", None, &gtdb, &ncbi)]);
        let err = resolve_meta(&mut meta, &InferenceOptions::default()).unwrap_err();
        assert!(err.contains("prediction looks like GTDB but gold standard looks like NCBI"));

        let mut meta = meta_for(&[("r1", Some("GTDB"), &gtdb, &ncbi)]);
        let options = InferenceOptions {
            allow_taxonomy_mismatch: true,
            ..Default::default()
        };
        let rows = resolve_meta(&mut meta, &options).unwrap();
        assert!(rows[0].warnings.iter().any(|w| w.kind == WARN_MISMATCH));
    }

    #[test]
    fn test_resolve_meta_keeps_given_taxonomy_but_warns_on_conflict() {
        let gtdb_a = write_profile("given_a.tsv", GTDB_LINES);
        let gtdb_b = write_profile("given_b.tsv", GTDB_LINES);
        let mut meta = meta_for(&[("r1", Some("NCBI"), &gtdb_a, &gtdb_b)]);

        let rows = resolve_meta(&mut meta, &InferenceOptions::default()).unwrap();

        assert_eq!(rows[0].taxonomy, "NCBI");
        assert_eq!(rows[0].source, "given");
        assert_eq!(rows[0].warnings.len(), 2);
        assert!(rows[0]
            .warnings
            .iter()
            .all(|w| w.kind == WARN_GIVEN_CONFLICT));
    }

    #[test]
    fn test_resolve_meta_stops_without_evidence() {
        let plain = write_profile("plain.tsv", &["Bacteria\t1.0", "Firmicutes\t0.5"]);
        let gtdb = write_profile("plain_gold.tsv", GTDB_LINES);
        let mut meta = meta_for(&[("r1", None, &plain, &gtdb)]);

        let err = resolve_meta(&mut meta, &InferenceOptions::default()).unwrap_err();
        assert!(
            err.contains("row 'r1': cannot tell the taxonomy"),
            "{}",
            err
        );
        assert!(err.contains("Add a Taxonomy column"), "{}", err);
    }

    #[test]
    fn test_index_identifies_release_and_release_gap() {
        let root = std::env::temp_dir().join(format!("benchpro_index_{}", std::process::id()));
        let dir = root.join("gtdb_index");
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(
            dir.join("version_columns.tsv"),
            "col_idx\tversion_id\n0\tg-r207\n1\tg-r214\n",
        )
        .unwrap();
        // Bit 0 = r207, bit 1 = r214. "Firmicutes_A" was renamed in r214.
        std::fs::write(
            dir.join("scientific_name_matrix.tsv"),
            "scientific_name\tbitset_hex\n\
             Firmicutes_A\t0000000000000001\n\
             Bacillota_A\t0000000000000002\n\
             Clostridia\t0000000000000003\n\
             Bacteria\t0000000000000003\n\
             Clostridium_K\t0000000000000003\n\
             Clostridium_K a\t0000000000000003\n",
        )
        .unwrap();

        let old = [
            "Firmicutes_A",
            "Clostridia",
            "Bacteria",
            "Clostridium_K",
            "Clostridium_K a",
        ];
        let new = [
            "Bacillota_A",
            "Clostridia",
            "Bacteria",
            "Clostridium_K",
            "Clostridium_K a",
        ];
        let to_set = |names: &[&str]| names.iter().map(|n| n.to_string()).collect::<HashSet<_>>();
        let all = to_set(&old).union(&to_set(&new)).cloned().collect();
        let index = TaxonomyIndex::load(&root, &all).unwrap();

        let text = |names: &[&str]| TextEvidence {
            names: to_set(names),
            ..Default::default()
        };
        let pred = infer(&text(&old), Some(&index));
        let gold = infer(&text(&new), Some(&index));
        std::fs::remove_dir_all(&root).unwrap();

        assert_eq!(pred.family, Some(TaxonomyFamily::GTDB));
        assert_eq!(pred.release.as_deref(), Some("GTDB r207"));
        assert_eq!(gold.release.as_deref(), Some("GTDB r214"));
        let (missing, release) = release_gap(
            pred.index_match.as_ref().unwrap(),
            gold.index_match.as_ref().unwrap(),
        )
        .unwrap();
        assert_eq!((missing, release.as_str()), (1, "g-r214"));
    }
}
