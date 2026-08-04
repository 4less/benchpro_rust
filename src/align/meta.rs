//! Samplesheet parsing for the `align` subcommand.
//!
//! A deliberately separate parser from [`crate::meta`]: that one requires `Profile`/`GoldStd`
//! columns and expands taxonomies, neither of which applies to an alignment run. What is shared is
//! the loading path (`.xlsx` via calamine, `.csv`/`.tsv` via Polars) and the convention that the
//! whole file is validated up front — a benchmark that dies on row 40 after an hour of work has
//! wasted the hour.

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use calamine::Xlsx;
use log::debug;
use polars::prelude::{CsvParseOptions, CsvReadOptions, DataFrame, SerReader};

use crate::options::{AlignArgs, ScoringMode};
use crate::utils::workbook_to_dataframe;

use super::error::{AlignError, AlignResult};

/// Column names recognised in an alignment samplesheet.
pub mod columns {
    /// Dataset label; results are grouped by it.
    pub const ID: &str = "ID";
    /// Sample id within the dataset; aggregation is across samples.
    pub const SAMPLE: &str = "Sample";
    /// Contender name.
    pub const TOOL: &str = "Tool";
    /// Path to the tool's `.sam`/`.sam.gz`/`.paf`/`.paf.gz`.
    pub const ALIGNMENT: &str = "Alignment";
    /// Path to the per-read truth TSV.
    pub const TRUTH: &str = "Truth";
    /// Reference FASTA; present means alignments are replayed against it.
    pub const REFERENCE: &str = "Reference";
    /// `contig<TAB>genome` TSV, required for `full` scoring.
    pub const CONTIG2GENOME: &str = "Contig2Genome";
    /// Tool to run the base-level head-to-head against.
    pub const PEER: &str = "Peer";
    /// Per-row scoring mode override.
    pub const SCORING: &str = "Scoring";
    /// Per-row marker-contig separator override.
    pub const SEP: &str = "Sep";

    /// Columns without which a row cannot be scored.
    pub const REQUIRED: &[&str] = &[ID, SAMPLE, TOOL, ALIGNMENT, TRUTH];
}

/// Which parser an alignment file needs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentFormat {
    /// SAM: full records, so base-level metrics are available.
    Sam,
    /// PAF: mapping-level only (no CIGAR, no SEQ, no NM).
    Paf,
}

impl AlignmentFormat {
    /// Detects the format from a path's extension, ignoring a trailing `.gz`.
    ///
    /// # Arguments
    ///
    /// * `path` - Alignment file path
    ///
    /// # Returns
    ///
    /// The detected format, or `None` when the extension is neither `.sam` nor `.paf`.
    pub fn from_path(path: &Path) -> Option<Self> {
        let name = path.file_name()?.to_str()?.to_ascii_lowercase();
        let stem = name.strip_suffix(".gz").unwrap_or(&name);
        if stem.ends_with(".sam") {
            Some(Self::Sam)
        } else if stem.ends_with(".paf") {
            Some(Self::Paf)
        } else {
            None
        }
    }

    /// Whether this format carries the CIGAR and SEQ that base-level metrics need.
    ///
    /// # Returns
    ///
    /// `true` for SAM, `false` for PAF.
    pub fn has_alignments(self) -> bool {
        matches!(self, Self::Sam)
    }
}

/// Is this alignment file gzip compressed?
///
/// # Arguments
///
/// * `path` - Alignment file path
///
/// # Returns
///
/// `true` when the path ends in `.gz`.
pub fn is_gzipped(path: &Path) -> bool {
    path.file_name()
        .and_then(|n| n.to_str())
        .is_some_and(|n| n.to_ascii_lowercase().ends_with(".gz"))
}

/// One scored contender on one sample.
#[derive(Debug, Clone)]
pub struct AlignRow {
    /// Dataset label.
    pub id: String,
    /// Sample id.
    pub sample: String,
    /// Contender name.
    pub tool: String,
    /// Alignment file.
    pub alignment: PathBuf,
    /// Format implied by the alignment file's extension.
    pub format: AlignmentFormat,
    /// Per-read truth file.
    pub truth: PathBuf,
    /// Reference FASTA, when the alignments should be replayed.
    pub reference: Option<PathBuf>,
    /// `contig -> genome` map, for `full` scoring.
    pub contig2genome: Option<PathBuf>,
    /// Tool this one is compared against at base level.
    pub peer: Option<String>,
    /// Effective scoring mode (row override, else the CLI default).
    pub scoring: ScoringMode,
    /// Effective marker-contig separator (row override, else the CLI default).
    pub sep: String,
}

impl AlignRow {
    /// The key that groups rows into one comparison: dataset and sample.
    ///
    /// # Returns
    ///
    /// A `(dataset, sample)` pair borrowed from the row.
    pub fn group_key(&self) -> (&str, &str) {
        (&self.id, &self.sample)
    }
}

/// A validated alignment samplesheet.
#[derive(Debug, Clone)]
pub struct AlignMeta {
    /// Rows in file order.
    pub rows: Vec<AlignRow>,
}

impl AlignMeta {
    /// Loads and fully validates a samplesheet.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to an `.xlsx`, `.csv` or `.tsv` samplesheet
    /// * `args` - CLI arguments supplying the defaults a row may override
    ///
    /// # Returns
    ///
    /// The parsed rows, in file order.
    ///
    /// # Errors
    ///
    /// Returns [`AlignError::Meta`] when the file cannot be read or a required column is missing,
    /// and [`AlignError::MetaValidation`] listing *every* problem when the rows themselves do not
    /// describe a runnable benchmark.
    pub fn from_path(path: &Path, args: &AlignArgs) -> AlignResult<Self> {
        let df = load_dataframe(path)?;
        let meta = Self::from_dataframe(path, df, args)?;
        meta.validate(path)?;
        Ok(meta)
    }

    /// Builds the row list from an already loaded DataFrame, without touching the filesystem.
    ///
    /// # Arguments
    ///
    /// * `path` - Samplesheet path, used only in error messages
    /// * `df` - The samplesheet as a DataFrame
    /// * `args` - CLI arguments supplying per-row defaults
    ///
    /// # Returns
    ///
    /// The parsed rows, in file order.
    ///
    /// # Errors
    ///
    /// Returns [`AlignError::Meta`] when a required column is missing or a cell cannot be read.
    pub fn from_dataframe(path: &Path, df: DataFrame, args: &AlignArgs) -> AlignResult<Self> {
        let present: HashMap<String, usize> = df
            .get_column_names()
            .iter()
            .enumerate()
            .map(|(i, name)| (name.to_lowercase(), i))
            .collect();

        let missing: Vec<&str> = columns::REQUIRED
            .iter()
            .copied()
            .filter(|c| !present.contains_key(&c.to_lowercase()))
            .collect();
        if !missing.is_empty() {
            return Err(AlignError::Meta {
                path: path.to_path_buf(),
                message: format!(
                    "missing required column(s): {}. Found: {}",
                    missing.join(", "),
                    df.get_column_names()
                        .iter()
                        .map(|n| n.to_string())
                        .collect::<Vec<_>>()
                        .join(", ")
                ),
            });
        }

        let column = |name: &str| -> Option<usize> { present.get(&name.to_lowercase()).copied() };
        let cell = |row: usize, name: &str| -> Option<String> {
            let idx = column(name)?;
            let value = df
                .get_columns()
                .get(idx)?
                .get(row)
                .ok()
                .map(|v| v.to_string())?;
            // Polars renders string cells quoted and nulls as "null"; both need undoing before
            // the value is used as a path.
            let trimmed = value.trim().trim_matches('"').trim();
            if trimmed.is_empty() || trimmed.eq_ignore_ascii_case("null") {
                None
            } else {
                Some(trimmed.to_string())
            }
        };

        let mut rows = Vec::with_capacity(df.height());
        let mut problems = Vec::new();

        for i in 0..df.height() {
            let row_no = i + 2; // 1-based, and the header occupies line 1.
            let mut required = |name: &str| -> Option<String> {
                match cell(i, name) {
                    Some(v) => Some(v),
                    None => {
                        problems.push(format!("row {row_no}: column '{name}' is empty"));
                        None
                    }
                }
            };
            let (id, sample, tool, alignment, truth) = (
                required(columns::ID),
                required(columns::SAMPLE),
                required(columns::TOOL),
                required(columns::ALIGNMENT),
                required(columns::TRUTH),
            );
            let (Some(id), Some(sample), Some(tool), Some(alignment), Some(truth)) =
                (id, sample, tool, alignment, truth)
            else {
                continue;
            };

            let alignment = PathBuf::from(alignment);
            let Some(format) = AlignmentFormat::from_path(&alignment) else {
                problems.push(format!(
                    "row {row_no}: cannot tell the format of '{}' — expected .sam, .paf, \
                     .sam.gz or .paf.gz",
                    alignment.display()
                ));
                continue;
            };

            let scoring = match cell(i, columns::SCORING).as_deref() {
                None => args.scoring,
                Some(s) if s.eq_ignore_ascii_case("full") => ScoringMode::Full,
                Some(s) if s.eq_ignore_ascii_case("species") => ScoringMode::Species,
                Some(other) => {
                    problems.push(format!(
                        "row {row_no}: Scoring '{other}' is not 'full' or 'species'"
                    ));
                    continue;
                }
            };

            rows.push(AlignRow {
                id,
                sample,
                tool,
                alignment,
                format,
                truth: PathBuf::from(truth),
                reference: cell(i, columns::REFERENCE).map(PathBuf::from),
                contig2genome: cell(i, columns::CONTIG2GENOME).map(PathBuf::from),
                peer: cell(i, columns::PEER),
                scoring,
                sep: cell(i, columns::SEP).unwrap_or_else(|| args.sep.clone()),
            });
        }

        if !problems.is_empty() {
            return Err(AlignError::MetaValidation {
                path: path.to_path_buf(),
                problems,
            });
        }

        debug!("Parsed {} alignment meta row(s)", rows.len());
        Ok(Self { rows })
    }

    /// Checks everything that can be checked before any parsing starts, reporting all failures.
    ///
    /// Covers: input files exist; `(ID, Sample, Tool)` is unique; a `Peer` names a tool present on
    /// the same `(ID, Sample)` and is not the row itself; `full` scoring has a `Contig2Genome`.
    ///
    /// # Arguments
    ///
    /// * `path` - Samplesheet path, for the error message
    ///
    /// # Returns
    ///
    /// `Ok(())` when every row is runnable.
    ///
    /// # Errors
    ///
    /// Returns [`AlignError::MetaValidation`] listing every problem found.
    pub fn validate(&self, path: &Path) -> AlignResult<()> {
        let mut problems = Vec::new();
        let mut seen: HashSet<(&str, &str, &str)> = HashSet::new();

        let mut tools_per_group: HashMap<(&str, &str), HashSet<&str>> = HashMap::new();
        for row in &self.rows {
            tools_per_group
                .entry(row.group_key())
                .or_default()
                .insert(row.tool.as_str());
        }

        for (i, row) in self.rows.iter().enumerate() {
            let row_no = i + 2;
            let mut check = |label: &str, p: &Path| {
                if !p.exists() {
                    problems.push(format!(
                        "row {row_no}: {label} '{}' does not exist",
                        p.display()
                    ));
                }
            };
            check(columns::ALIGNMENT, &row.alignment);
            check(columns::TRUTH, &row.truth);
            if let Some(reference) = &row.reference {
                check(columns::REFERENCE, reference);
            }
            if let Some(c2g) = &row.contig2genome {
                check(columns::CONTIG2GENOME, c2g);
            }

            if !seen.insert((&row.id, &row.sample, &row.tool)) {
                problems.push(format!(
                    "row {row_no}: duplicate (ID, Sample, Tool) = ({}, {}, {})",
                    row.id, row.sample, row.tool
                ));
            }

            if row.scoring == ScoringMode::Full && row.contig2genome.is_none() {
                problems.push(format!(
                    "row {row_no}: scoring 'full' needs a {} column (the genome stratum has no \
                     other source); use scoring 'species' for a marker reference",
                    columns::CONTIG2GENOME
                ));
            }

            if let Some(peer) = &row.peer {
                if peer == &row.tool {
                    problems.push(format!(
                        "row {row_no}: Peer '{peer}' is the row's own tool; a head-to-head needs \
                         two different tools"
                    ));
                } else if !tools_per_group
                    .get(&row.group_key())
                    .is_some_and(|tools| tools.contains(peer.as_str()))
                {
                    problems.push(format!(
                        "row {row_no}: Peer '{peer}' is not a Tool on (ID, Sample) = ({}, {})",
                        row.id, row.sample
                    ));
                }
            }
        }

        if self.rows.is_empty() {
            problems.push("the samplesheet has no rows".to_string());
        }

        if problems.is_empty() {
            Ok(())
        } else {
            Err(AlignError::MetaValidation {
                path: path.to_path_buf(),
                problems,
            })
        }
    }

    /// Distinct dataset labels, in first-seen order.
    ///
    /// # Returns
    ///
    /// Every distinct `ID`.
    pub fn datasets(&self) -> Vec<&str> {
        unique_in_order(self.rows.iter().map(|r| r.id.as_str()))
    }

    /// Distinct tool names, in first-seen order.
    ///
    /// # Returns
    ///
    /// Every distinct `Tool`.
    pub fn tools(&self) -> Vec<&str> {
        unique_in_order(self.rows.iter().map(|r| r.tool.as_str()))
    }

    /// Rows grouped into the units that are scored together, in first-seen order.
    ///
    /// A head-to-head needs every contender on a `(dataset, sample)` in memory at once, which is
    /// what makes the group — not the individual row — the unit of work.
    ///
    /// # Returns
    ///
    /// One entry per `(dataset, sample)`, each holding that group's rows in file order.
    pub fn groups(&self) -> Vec<((&str, &str), Vec<&AlignRow>)> {
        let mut order: Vec<(&str, &str)> = Vec::new();
        let mut grouped: HashMap<(&str, &str), Vec<&AlignRow>> = HashMap::new();
        for row in &self.rows {
            let key = row.group_key();
            if !grouped.contains_key(&key) {
                order.push(key);
            }
            grouped.entry(key).or_default().push(row);
        }
        order
            .into_iter()
            .map(|key| {
                let rows = grouped.remove(&key).unwrap_or_default();
                (key, rows)
            })
            .collect()
    }
}

/// Distinct items of an iterator, keeping first-seen order.
fn unique_in_order<'a>(items: impl Iterator<Item = &'a str>) -> Vec<&'a str> {
    let mut seen = HashSet::new();
    items.filter(|item| seen.insert(*item)).collect()
}

/// Reads a samplesheet into a DataFrame, dispatching on the file extension.
fn load_dataframe(path: &Path) -> AlignResult<DataFrame> {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_ascii_lowercase());

    match ext.as_deref() {
        Some("xlsx") => {
            let mut workbook: Xlsx<_> =
                calamine::open_workbook(path).map_err(|e| AlignError::Meta {
                    path: path.to_path_buf(),
                    message: format!("cannot open workbook: {e}"),
                })?;
            workbook_to_dataframe(&mut workbook).map_err(|e| AlignError::Meta {
                path: path.to_path_buf(),
                message: format!("cannot read workbook: {e}"),
            })
        }
        Some("tsv") | Some("csv") => {
            let separator = if ext.as_deref() == Some("tsv") {
                b'\t'
            } else {
                b','
            };
            CsvReadOptions::default()
                .with_parse_options(
                    CsvParseOptions::default()
                        .with_separator(separator)
                        .with_comment_prefix(Some("#")),
                )
                .try_into_reader_with_file_path(Some(path.to_path_buf()))
                .and_then(|reader| reader.finish())
                .map_err(|e| AlignError::Meta {
                    path: path.to_path_buf(),
                    message: format!("cannot read table: {e}"),
                })
        }
        other => Err(AlignError::Meta {
            path: path.to_path_buf(),
            message: format!(
                "unsupported extension {:?}; expected .xlsx, .csv or .tsv",
                other.unwrap_or("<none>")
            ),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use polars::prelude::{NamedFrom, Series};

    fn args() -> AlignArgs {
        AlignArgs {
            meta: PathBuf::from("meta.tsv"),
            outprefix: "out".to_string(),
            scoring: ScoringMode::Species,
            sep: "_".to_string(),
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

    fn df(columns: &[(&str, Vec<&str>)]) -> DataFrame {
        DataFrame::new(
            columns
                .iter()
                .map(|(name, values)| Series::new((*name).into(), values.clone()))
                .collect(),
        )
        .expect("test frame")
    }

    fn minimal() -> DataFrame {
        df(&[
            ("ID", vec!["ds", "ds"]),
            ("Sample", vec!["s1", "s1"]),
            ("Tool", vec!["flexalign", "minibwa"]),
            ("Alignment", vec!["a.sam", "b.sam.gz"]),
            ("Truth", vec!["t.tsv", "t.tsv"]),
        ])
    }

    #[test]
    fn parses_required_columns_and_applies_cli_defaults() {
        let meta = AlignMeta::from_dataframe(Path::new("meta.tsv"), minimal(), &args()).unwrap();

        assert_eq!(meta.rows.len(), 2);
        assert_eq!(meta.rows[0].tool, "flexalign");
        assert_eq!(meta.rows[0].format, AlignmentFormat::Sam);
        assert_eq!(meta.rows[0].scoring, ScoringMode::Species);
        assert_eq!(meta.rows[0].sep, "_");
        assert!(meta.rows[0].reference.is_none());
        assert_eq!(meta.tools(), vec!["flexalign", "minibwa"]);
        assert_eq!(meta.datasets(), vec!["ds"]);
    }

    #[test]
    fn column_names_are_case_insensitive() {
        let frame = df(&[
            ("id", vec!["ds"]),
            ("sample", vec!["s1"]),
            ("TOOL", vec!["flexalign"]),
            ("alignment", vec!["a.sam"]),
            ("truth", vec!["t.tsv"]),
        ]);
        let meta = AlignMeta::from_dataframe(Path::new("meta.tsv"), frame, &args()).unwrap();
        assert_eq!(meta.rows[0].id, "ds");
    }

    #[test]
    fn missing_required_column_is_reported_with_what_was_found() {
        let frame = df(&[("ID", vec!["ds"]), ("Sample", vec!["s1"])]);
        let err = AlignMeta::from_dataframe(Path::new("meta.tsv"), frame, &args()).unwrap_err();
        let text = err.to_string();
        assert!(text.contains("Tool"), "{text}");
        assert!(text.contains("Alignment"), "{text}");
        assert!(text.contains("Found: ID, Sample"), "{text}");
    }

    #[test]
    fn row_overrides_beat_cli_defaults() {
        let frame = df(&[
            ("ID", vec!["ds"]),
            ("Sample", vec!["s1"]),
            ("Tool", vec!["flexalign"]),
            ("Alignment", vec!["a.sam"]),
            ("Truth", vec!["t.tsv"]),
            ("Scoring", vec!["full"]),
            ("Sep", vec!["|"]),
        ]);
        let meta = AlignMeta::from_dataframe(Path::new("meta.tsv"), frame, &args()).unwrap();
        assert_eq!(meta.rows[0].scoring, ScoringMode::Full);
        assert_eq!(meta.rows[0].sep, "|");
    }

    #[test]
    fn unknown_alignment_extension_is_rejected() {
        let frame = df(&[
            ("ID", vec!["ds"]),
            ("Sample", vec!["s1"]),
            ("Tool", vec!["flexalign"]),
            ("Alignment", vec!["a.bam"]),
            ("Truth", vec!["t.tsv"]),
        ]);
        let err = AlignMeta::from_dataframe(Path::new("meta.tsv"), frame, &args()).unwrap_err();
        assert!(err.to_string().contains("cannot tell the format"));
    }

    #[test]
    fn validation_reports_missing_files_duplicates_and_bad_peers_together() {
        let frame = df(&[
            ("ID", vec!["ds", "ds"]),
            ("Sample", vec!["s1", "s1"]),
            ("Tool", vec!["flexalign", "flexalign"]),
            ("Alignment", vec!["nope.sam", "nope2.sam"]),
            ("Truth", vec!["nope.tsv", "nope.tsv"]),
            ("Peer", vec!["ghost", "flexalign"]),
        ]);
        let meta = AlignMeta::from_dataframe(Path::new("meta.tsv"), frame, &args()).unwrap();
        let err = meta.validate(Path::new("meta.tsv")).unwrap_err();
        let text = err.to_string();

        assert!(text.contains("nope.sam' does not exist"), "{text}");
        assert!(text.contains("duplicate (ID, Sample, Tool)"), "{text}");
        assert!(text.contains("Peer 'ghost' is not a Tool"), "{text}");
        assert!(text.contains("is the row's own tool"), "{text}");
    }

    #[test]
    fn full_scoring_without_contig2genome_is_rejected() {
        let frame = df(&[
            ("ID", vec!["ds"]),
            ("Sample", vec!["s1"]),
            ("Tool", vec!["flexalign"]),
            ("Alignment", vec!["a.sam"]),
            ("Truth", vec!["t.tsv"]),
            ("Scoring", vec!["full"]),
        ]);
        let meta = AlignMeta::from_dataframe(Path::new("meta.tsv"), frame, &args()).unwrap();
        let err = meta.validate(Path::new("meta.tsv")).unwrap_err();
        assert!(err.to_string().contains("needs a Contig2Genome column"));
    }

    #[test]
    fn groups_keep_file_order_and_hold_every_contender() {
        let frame = df(&[
            ("ID", vec!["ds", "ds", "ds"]),
            ("Sample", vec!["s1", "s2", "s1"]),
            ("Tool", vec!["a", "a", "b"]),
            ("Alignment", vec!["a1.sam", "a2.sam", "b1.sam"]),
            ("Truth", vec!["t1.tsv", "t2.tsv", "t1.tsv"]),
        ]);
        let meta = AlignMeta::from_dataframe(Path::new("meta.tsv"), frame, &args()).unwrap();
        let groups = meta.groups();

        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].0, ("ds", "s1"));
        assert_eq!(groups[0].1.len(), 2);
        assert_eq!(groups[1].0, ("ds", "s2"));
        assert_eq!(groups[1].1.len(), 1);
    }

    #[test]
    fn format_detection_ignores_gz_and_case() {
        assert_eq!(
            AlignmentFormat::from_path(Path::new("x/y.SAM.gz")),
            Some(AlignmentFormat::Sam)
        );
        assert_eq!(
            AlignmentFormat::from_path(Path::new("x/y.paf")),
            Some(AlignmentFormat::Paf)
        );
        assert_eq!(AlignmentFormat::from_path(Path::new("x/y.bam")), None);
        assert!(AlignmentFormat::Sam.has_alignments());
        assert!(!AlignmentFormat::Paf.has_alignments());
        assert!(is_gzipped(Path::new("a.sam.gz")));
        assert!(!is_gzipped(Path::new("a.sam")));
    }
}
