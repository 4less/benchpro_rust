use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::path::{Path, PathBuf};

use log::{info, warn};
use phylotree::tree::{NodeId, Tree};
use polars::prelude::{CsvWriter, DataFrame, NamedFrom, SerWriter, Series};
use rayon::prelude::*;
use thiserror::Error;

use crate::meta::Meta;
use crate::msa::{read_fasta_alignment, MsaError};
use crate::options::StrainArgs;

/// Errors produced while processing strain inputs.
#[derive(Debug, Error)]
pub enum StrainError {
    /// IO errors when reading input or writing output.
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    /// Errors from Polars when reading or writing tables.
    #[error("Polars error: {0}")]
    Polars(#[from] polars::error::PolarsError),
    /// Errors from parsing Newick trees.
    #[error("Tree parse error: {0}")]
    Tree(#[from] phylotree::tree::NewickParseError),
    /// Errors from input meta files or missing data.
    #[error("{0}")]
    Meta(String),
    /// Errors from MSA file processing.
    #[error("MSA error: {0}")]
    Msa(#[from] MsaError),
}

struct StrainMetaColumns;
impl StrainMetaColumns {
    const ID: &'static str = "ID";
    const SPECIES: &'static str = "Species";
    const TREE: &'static str = "Tree";
    const META: &'static str = "Meta";
    const MSA: &'static str = "MSA";
    const GOLD_MSA: &'static str = "GoldMSA";
}

struct SampleMetaColumns;
impl SampleMetaColumns {
    const ID: &'static str = "ID";
    const GENOME: &'static str = "genome";
    const COVERAGE: &'static str = "coverage";
}

struct SampleInfo {
    id: String,
    coverage: f64,
}

struct StrainJob {
    id: String,
    species: String,
    tree_path: PathBuf,
    meta_path: PathBuf,
    msa_path: Option<PathBuf>,
    gold_msa_path: Option<PathBuf>,
}

struct MonophylyRow {
    id: String,
    species: String,
    genome: String,
    tip_count: usize,
    lca_tip_count: usize,
    monophyly_score: f64,
    pw_dist_min: f64,
    pw_dist_max: f64,
    pw_dist_mean: f64,
    pw_dist_median: f64,
}

struct TipRow {
    id: String,
    species: String,
    genome: String,
    sample: String,
    coverage: f64,
    tip_count: usize,
    lca_tip_count: usize,
    monophyly_score: f64,
    pw_dist_min: f64,
    pw_dist_max: f64,
    pw_dist_mean: f64,
    pw_dist_median: f64,
}

/// Minimum valid pairs needed to compute error scores for a sample.
const MIN_PAIRS: usize = 2;

/// Small constant added to distances before log/weight operations.
const MSA_EPSILON: f64 = 1e-6;

struct SampleErrorRow {
    id: String,
    species: String,
    sample: String,
    genome: String,
    coverage: f64,
    pair_count: usize,
    mean_gold_dist: f64,
    mean_sample_dist: f64,
    /// Weighted Pearson r of log(d_gold) vs log(d_sample), weight = 1/(d_gold+ε).
    weighted_log_pearson_r: f64,
    /// 1 − weighted_log_pearson_r (0 = perfect, larger = worse).
    weighted_log_pearson_error: f64,
    /// Spearman rank correlation of d_gold vs d_sample.
    spearman_r: f64,
    /// 1 − spearman_r.
    spearman_error: f64,
    /// Weighted Kendall τ: concordance of pairwise orderings, weight = 1/(min_gold_dist+ε).
    weighted_kendall_tau: f64,
    /// 1 − weighted_kendall_tau.
    weighted_kendall_error: f64,
}

/// Run the `strain` subcommand.
///
/// Loads each row from the samplesheet, reads the associated phylogenetic tree and
/// sample-to-genome mapping, then computes per-genome monophyly scores and pairwise
/// branch-length distances.  Results are written to `{outprefix}.monophyly.tsv`.
///
/// # Arguments
///
/// * `args` - Parsed strain CLI arguments
///
/// # Errors
///
/// Returns [`StrainError`] if any input file cannot be read or parsed.
pub fn run_strain(args: &StrainArgs) -> Result<(), StrainError> {
    let jobs = load_strain_jobs(&args.meta)?;
    let mut rows = Vec::new();
    let mut tip_rows = Vec::new();
    let mut sample_error_rows: Vec<SampleErrorRow> = Vec::new();

    for job in &jobs {
        info!(
            "Processing '{}' (tree: {})",
            job.id,
            job.tree_path.display()
        );
        let tree = load_tree(&job.tree_path)?;
        let genome_groups = load_genome_groups(&job.meta_path)?;

        let name2id: HashMap<String, NodeId> = tree
            .get_leaves()
            .iter()
            .filter_map(|id| {
                tree.get(id)
                    .ok()
                    .and_then(|n| n.name.as_ref().map(|name| (name.clone(), *id)))
            })
            .collect();

        for (genome, samples) in &genome_groups {
            // Resolve samples to tree node IDs, keeping coverage alongside.
            let resolved: Vec<(&SampleInfo, NodeId)> = samples
                .iter()
                .filter_map(|s| {
                    let node_id = name2id.get(&s.id);
                    if node_id.is_none() {
                        warn!(
                            "Sample '{}' not found as tip in tree '{}'",
                            s.id,
                            job.tree_path.display()
                        );
                    }
                    node_id.map(|nid| (s, *nid))
                })
                .collect();

            let tip_count = resolved.len();
            if tip_count == 0 {
                warn!(
                    "No tree tips found for genome '{}' in '{}', skipping.",
                    genome,
                    job.tree_path.display()
                );
                continue;
            }

            let tip_ids: Vec<NodeId> = resolved.iter().map(|(_, nid)| *nid).collect();

            let (monophyly_score, lca_tip_count) = if tip_count == 1 {
                (1.0_f64, 1_usize)
            } else {
                let lca = find_lca(&tree, &tip_ids)?;
                let lca_leaves = tree
                    .get_subtree_leaves(&lca)
                    .map_err(|e| StrainError::Meta(e.to_string()))?;
                let lca_tip_count = lca_leaves.len();
                (tip_count as f64 / lca_tip_count as f64, lca_tip_count)
            };

            let (pw_min, pw_max, pw_mean, pw_median) =
                compute_pairwise_distances(&tree, &tip_ids)?;

            rows.push(MonophylyRow {
                id: job.id.clone(),
                species: job.species.clone(),
                genome: genome.clone(),
                tip_count,
                lca_tip_count,
                monophyly_score,
                pw_dist_min: pw_min,
                pw_dist_max: pw_max,
                pw_dist_mean: pw_mean,
                pw_dist_median: pw_median,
            });

            // Per-tip rows: distances from this tip to all others in the group.
            for (i, (sample, nid)) in resolved.iter().enumerate() {
                let others: Vec<NodeId> = tip_ids
                    .iter()
                    .enumerate()
                    .filter(|(j, _)| *j != i)
                    .map(|(_, id)| *id)
                    .collect();

                let (tip_min, tip_max, tip_mean, tip_median) = if others.is_empty() {
                    (0.0, 0.0, 0.0, 0.0)
                } else {
                    compute_distances_from_one(&tree, nid, &others)?
                };

                tip_rows.push(TipRow {
                    id: job.id.clone(),
                    species: job.species.clone(),
                    genome: genome.clone(),
                    sample: sample.id.clone(),
                    coverage: sample.coverage,
                    tip_count,
                    lca_tip_count,
                    monophyly_score,
                    pw_dist_min: tip_min,
                    pw_dist_max: tip_max,
                    pw_dist_mean: tip_mean,
                    pw_dist_median: tip_median,
                });
            }
        }

        match (&job.msa_path, &job.gold_msa_path) {
            (Some(msa_path), Some(gold_msa_path)) => {
                info!("Computing MSA distance errors for '{}'...", job.id);
                match compute_msa_error_rows(job, msa_path, gold_msa_path, &genome_groups) {
                    Ok(error_rows) => sample_error_rows.extend(error_rows),
                    Err(e) => warn!(
                        "MSA error computation failed for '{}': {}",
                        job.id, e
                    ),
                }
            }
            (Some(_), None) | (None, Some(_)) => {
                warn!(
                    "'{}' has only one of MSA/GoldMSA defined; skipping distance error analysis.",
                    job.id
                );
            }
            (None, None) => {}
        }
    }

    let output = PathBuf::from(format!("{}.monophyly.tsv", args.outprefix));
    write_monophyly_output(&rows, &output)?;
    info!("Wrote monophyly stats to '{}'.", output.display());

    let tip_output = PathBuf::from(format!("{}.monophyly_tips.tsv", args.outprefix));
    write_tip_output(&tip_rows, &tip_output)?;
    info!("Wrote per-tip monophyly stats to '{}'.", tip_output.display());

    if !sample_error_rows.is_empty() {
        let error_output = PathBuf::from(format!("{}.sample_error.tsv", args.outprefix));
        write_sample_error_output(&sample_error_rows, &error_output)?;
        info!("Wrote sample error stats to '{}'.", error_output.display());
    }

    Ok(())
}

fn find_lca(tree: &Tree, tip_ids: &[NodeId]) -> Result<NodeId, StrainError> {
    let mut lca = tip_ids[0];
    for &tid in &tip_ids[1..] {
        lca = tree
            .get_common_ancestor(&lca, &tid)
            .map_err(|e| StrainError::Meta(e.to_string()))?;
    }
    Ok(lca)
}

fn compute_pairwise_distances(
    tree: &Tree,
    tip_ids: &[NodeId],
) -> Result<(f64, f64, f64, f64), StrainError> {
    if tip_ids.len() == 1 {
        return Ok((0.0, 0.0, 0.0, 0.0));
    }

    let mut distances = Vec::new();
    for i in 0..tip_ids.len() {
        for j in (i + 1)..tip_ids.len() {
            let (dist, _) = tree
                .get_distance(&tip_ids[i], &tip_ids[j])
                .map_err(|e| StrainError::Meta(e.to_string()))?;
            if let Some(d) = dist {
                distances.push(d);
            }
        }
    }

    if distances.is_empty() {
        return Ok((0.0, 0.0, 0.0, 0.0));
    }

    let min = distances.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = distances.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mean = distances.iter().sum::<f64>() / distances.len() as f64;

    let mut sorted = distances.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = if sorted.len() % 2 == 0 {
        let mid = sorted.len() / 2;
        (sorted[mid - 1] + sorted[mid]) / 2.0
    } else {
        sorted[sorted.len() / 2]
    };

    Ok((min, max, mean, median))
}

fn compute_distances_from_one(
    tree: &Tree,
    source: &NodeId,
    targets: &[NodeId],
) -> Result<(f64, f64, f64, f64), StrainError> {
    let mut distances = Vec::with_capacity(targets.len());
    for target in targets {
        let (dist, _) = tree
            .get_distance(source, target)
            .map_err(|e| StrainError::Meta(e.to_string()))?;
        if let Some(d) = dist {
            distances.push(d);
        }
    }

    if distances.is_empty() {
        return Ok((0.0, 0.0, 0.0, 0.0));
    }

    let min = distances.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = distances.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mean = distances.iter().sum::<f64>() / distances.len() as f64;

    let mut sorted = distances.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = if sorted.len() % 2 == 0 {
        let mid = sorted.len() / 2;
        (sorted[mid - 1] + sorted[mid]) / 2.0
    } else {
        sorted[sorted.len() / 2]
    };

    Ok((min, max, mean, median))
}

fn load_tree(path: &Path) -> Result<Tree, StrainError> {
    let raw = fs::read_to_string(path)?;
    Ok(Tree::from_newick(raw.trim())?)
}

fn normalize_header(name: &str) -> String {
    name.trim()
        .trim_start_matches('\u{feff}')
        .to_ascii_lowercase()
}

fn find_column_name(df: &DataFrame, expected: &str) -> Option<String> {
    let expected_norm = normalize_header(expected);
    df.get_column_names()
        .iter()
        .find(|name| normalize_header(name) == expected_norm)
        .map(|name| name.to_string())
}

fn load_strain_jobs(path: &Path) -> Result<Vec<StrainJob>, StrainError> {
    let df = Meta::polars_from_path(path).ok_or_else(|| {
        StrainError::Meta(format!("Failed to read meta file '{}'", path.display()))
    })?;

    let id_col = find_column_name(&df, StrainMetaColumns::ID).ok_or_else(|| {
        StrainError::Meta(format!(
            "Missing column '{}' in '{}'",
            StrainMetaColumns::ID,
            path.display()
        ))
    })?;
    let species_col = find_column_name(&df, StrainMetaColumns::SPECIES).ok_or_else(|| {
        StrainError::Meta(format!(
            "Missing column '{}' in '{}'",
            StrainMetaColumns::SPECIES,
            path.display()
        ))
    })?;
    let tree_col = find_column_name(&df, StrainMetaColumns::TREE).ok_or_else(|| {
        StrainError::Meta(format!(
            "Missing column '{}' in '{}'",
            StrainMetaColumns::TREE,
            path.display()
        ))
    })?;
    let meta_col = find_column_name(&df, StrainMetaColumns::META).ok_or_else(|| {
        StrainError::Meta(format!(
            "Missing column '{}' in '{}'",
            StrainMetaColumns::META,
            path.display()
        ))
    })?;

    let msa_col = find_column_name(&df, StrainMetaColumns::MSA);
    let gold_msa_col = find_column_name(&df, StrainMetaColumns::GOLD_MSA);

    let ids = df.column(&id_col)?.str().map_err(StrainError::Polars)?;
    let species_vals = df
        .column(&species_col)?
        .str()
        .map_err(StrainError::Polars)?;
    let trees = df.column(&tree_col)?.str().map_err(StrainError::Polars)?;
    let metas = df.column(&meta_col)?.str().map_err(StrainError::Polars)?;

    let mut jobs = Vec::with_capacity(df.height());
    for row in 0..df.height() {
        let id = ids
            .get(row)
            .ok_or_else(|| StrainError::Meta("Missing ID value".to_string()))?
            .to_string();
        let species = species_vals
            .get(row)
            .ok_or_else(|| StrainError::Meta("Missing Species value".to_string()))?
            .to_string();
        let tree_path = trees
            .get(row)
            .ok_or_else(|| StrainError::Meta("Missing Tree value".to_string()))?
            .to_string();
        let meta_path = metas
            .get(row)
            .ok_or_else(|| StrainError::Meta("Missing Meta value".to_string()))?
            .to_string();

        let msa_path = msa_col
            .as_ref()
            .and_then(|col| df.column(col).ok())
            .and_then(|s| s.str().ok().map(|ca| ca.get(row).map(PathBuf::from)))
            .flatten();
        let gold_msa_path = gold_msa_col
            .as_ref()
            .and_then(|col| df.column(col).ok())
            .and_then(|s| s.str().ok().map(|ca| ca.get(row).map(PathBuf::from)))
            .flatten();

        jobs.push(StrainJob {
            id,
            species,
            tree_path: PathBuf::from(tree_path),
            meta_path: PathBuf::from(meta_path),
            msa_path,
            gold_msa_path,
        });
    }

    Ok(jobs)
}

fn load_genome_groups(
    path: &Path,
) -> Result<BTreeMap<String, Vec<SampleInfo>>, StrainError> {
    let df = Meta::polars_from_path(path).ok_or_else(|| {
        StrainError::Meta(format!(
            "Failed to read sample meta file '{}'",
            path.display()
        ))
    })?;

    let id_col = find_column_name(&df, SampleMetaColumns::ID).ok_or_else(|| {
        StrainError::Meta(format!(
            "Missing column '{}' in '{}'",
            SampleMetaColumns::ID,
            path.display()
        ))
    })?;
    let genome_col = find_column_name(&df, SampleMetaColumns::GENOME).ok_or_else(|| {
        StrainError::Meta(format!(
            "Missing column '{}' in '{}'",
            SampleMetaColumns::GENOME,
            path.display()
        ))
    })?;
    let coverage_col = find_column_name(&df, SampleMetaColumns::COVERAGE).ok_or_else(|| {
        StrainError::Meta(format!(
            "Missing column '{}' in '{}'",
            SampleMetaColumns::COVERAGE,
            path.display()
        ))
    })?;

    let ids = df.column(&id_col)?.str().map_err(StrainError::Polars)?;
    let genomes = df.column(&genome_col)?.str().map_err(StrainError::Polars)?;
    let coverages = df
        .column(&coverage_col)?
        .cast(&polars::prelude::DataType::Float64)
        .map_err(StrainError::Polars)?;
    let coverages = coverages.f64().map_err(StrainError::Polars)?;

    let mut groups: BTreeMap<String, Vec<SampleInfo>> = BTreeMap::new();
    for row in 0..df.height() {
        let id = ids
            .get(row)
            .ok_or_else(|| StrainError::Meta("Missing sample ID value".to_string()))?
            .to_string();
        let genome = genomes
            .get(row)
            .ok_or_else(|| StrainError::Meta("Missing genome value".to_string()))?
            .to_string();
        let coverage = coverages.get(row).unwrap_or(f64::NAN);
        groups.entry(genome).or_default().push(SampleInfo { id, coverage });
    }

    Ok(groups)
}

fn write_monophyly_output(rows: &[MonophylyRow], output: &Path) -> Result<(), StrainError> {
    let mut ids = Vec::with_capacity(rows.len());
    let mut species = Vec::with_capacity(rows.len());
    let mut genomes = Vec::with_capacity(rows.len());
    let mut tip_counts = Vec::with_capacity(rows.len());
    let mut lca_tip_counts = Vec::with_capacity(rows.len());
    let mut monophyly_scores = Vec::with_capacity(rows.len());
    let mut pw_mins = Vec::with_capacity(rows.len());
    let mut pw_maxs = Vec::with_capacity(rows.len());
    let mut pw_means = Vec::with_capacity(rows.len());
    let mut pw_medians = Vec::with_capacity(rows.len());

    for row in rows {
        ids.push(row.id.clone());
        species.push(row.species.clone());
        genomes.push(row.genome.clone());
        tip_counts.push(row.tip_count as u64);
        lca_tip_counts.push(row.lca_tip_count as u64);
        monophyly_scores.push(row.monophyly_score);
        pw_mins.push(row.pw_dist_min);
        pw_maxs.push(row.pw_dist_max);
        pw_means.push(row.pw_dist_mean);
        pw_medians.push(row.pw_dist_median);
    }

    let mut df = DataFrame::new(vec![
        Series::new("ID".into(), ids),
        Series::new("Species".into(), species),
        Series::new("Genome".into(), genomes),
        Series::new("TipCount".into(), tip_counts),
        Series::new("LcaTipCount".into(), lca_tip_counts),
        Series::new("MonophylyScore".into(), monophyly_scores),
        Series::new("PwDistMin".into(), pw_mins),
        Series::new("PwDistMax".into(), pw_maxs),
        Series::new("PwDistMean".into(), pw_means),
        Series::new("PwDistMedian".into(), pw_medians),
    ])
    .map_err(StrainError::Polars)?;

    let mut file = std::fs::File::create(output)?;
    CsvWriter::new(&mut file)
        .include_header(true)
        .with_separator(b'\t')
        .finish(&mut df)
        .map_err(StrainError::Polars)?;

    Ok(())
}

fn write_tip_output(rows: &[TipRow], output: &Path) -> Result<(), StrainError> {
    let mut ids = Vec::with_capacity(rows.len());
    let mut species = Vec::with_capacity(rows.len());
    let mut genomes = Vec::with_capacity(rows.len());
    let mut samples = Vec::with_capacity(rows.len());
    let mut coverages = Vec::with_capacity(rows.len());
    let mut tip_counts = Vec::with_capacity(rows.len());
    let mut lca_tip_counts = Vec::with_capacity(rows.len());
    let mut monophyly_scores = Vec::with_capacity(rows.len());
    let mut pw_mins = Vec::with_capacity(rows.len());
    let mut pw_maxs = Vec::with_capacity(rows.len());
    let mut pw_means = Vec::with_capacity(rows.len());
    let mut pw_medians = Vec::with_capacity(rows.len());

    for row in rows {
        ids.push(row.id.clone());
        species.push(row.species.clone());
        genomes.push(row.genome.clone());
        samples.push(row.sample.clone());
        coverages.push(row.coverage);
        tip_counts.push(row.tip_count as u64);
        lca_tip_counts.push(row.lca_tip_count as u64);
        monophyly_scores.push(row.monophyly_score);
        pw_mins.push(row.pw_dist_min);
        pw_maxs.push(row.pw_dist_max);
        pw_means.push(row.pw_dist_mean);
        pw_medians.push(row.pw_dist_median);
    }

    let mut df = DataFrame::new(vec![
        Series::new("ID".into(), ids),
        Series::new("Species".into(), species),
        Series::new("Genome".into(), genomes),
        Series::new("Sample".into(), samples),
        Series::new("Coverage".into(), coverages),
        Series::new("TipCount".into(), tip_counts),
        Series::new("LcaTipCount".into(), lca_tip_counts),
        Series::new("MonophylyScore".into(), monophyly_scores),
        Series::new("PwDistMin".into(), pw_mins),
        Series::new("PwDistMax".into(), pw_maxs),
        Series::new("PwDistMean".into(), pw_means),
        Series::new("PwDistMedian".into(), pw_medians),
    ])
    .map_err(StrainError::Polars)?;

    let mut file = std::fs::File::create(output)?;
    CsvWriter::new(&mut file)
        .include_header(true)
        .with_separator(b'\t')
        .finish(&mut df)
        .map_err(StrainError::Polars)?;

    Ok(())
}

fn p_distance(seq_a: &[u8], seq_b: &[u8]) -> Option<f64> {
    let mut informative = 0u64;
    let mut differ = 0u64;
    for (&a, &b) in seq_a.iter().zip(seq_b.iter()) {
        let a_info = a != b'-' && a != b'N';
        let b_info = b != b'-' && b != b'N';
        if a_info && b_info {
            informative += 1;
            if a != b {
                differ += 1;
            }
        }
    }
    if informative == 0 {
        None
    } else {
        Some(differ as f64 / informative as f64)
    }
}

fn compute_pairwise_msa_distances(
    sequences: &HashMap<String, Vec<u8>>,
) -> HashMap<(String, String), f64> {
    let mut keys: Vec<&String> = sequences.keys().collect();
    keys.sort();
    let n = keys.len();

    let pairs: Vec<(usize, usize)> = (0..n)
        .flat_map(|i| (i + 1..n).map(move |j| (i, j)))
        .collect();

    pairs
        .par_iter()
        .filter_map(|&(i, j)| {
            let a = keys[i];
            let b = keys[j];
            let seq_a = &sequences[a];
            let seq_b = &sequences[b];
            p_distance(seq_a, seq_b).map(|dist| ((a.clone(), b.clone()), dist))
        })
        .collect()
}

fn lookup_msa_dist(map: &HashMap<(String, String), f64>, a: &str, b: &str) -> Option<f64> {
    let key = if a <= b {
        (a.to_string(), b.to_string())
    } else {
        (b.to_string(), a.to_string())
    };
    map.get(&key).copied()
}

fn rank_vector(values: &[f64]) -> Vec<f64> {
    let n = values.len();
    let mut indexed: Vec<(usize, f64)> = values.iter().copied().enumerate().collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    let mut ranks = vec![0.0f64; n];
    let mut i = 0;
    while i < n {
        let mut j = i;
        while j < n && (indexed[j].1 - indexed[i].1).abs() < 1e-15 {
            j += 1;
        }
        let avg_rank = (i as f64 + (j - 1) as f64) / 2.0 + 1.0;
        for k in i..j {
            ranks[indexed[k].0] = avg_rank;
        }
        i = j;
    }
    ranks
}

fn pearson_correlation(x: &[f64], y: &[f64]) -> Option<f64> {
    let n = x.len();
    if n < 2 {
        return None;
    }
    let mx = x.iter().sum::<f64>() / n as f64;
    let my = y.iter().sum::<f64>() / n as f64;
    let cov: f64 = x.iter().zip(y.iter()).map(|(xi, yi)| (xi - mx) * (yi - my)).sum();
    let var_x: f64 = x.iter().map(|xi| (xi - mx).powi(2)).sum();
    let var_y: f64 = y.iter().map(|yi| (yi - my).powi(2)).sum();
    let denom = (var_x * var_y).sqrt();
    if denom < 1e-15 {
        return None;
    }
    Some(cov / denom)
}

fn weighted_pearson_correlation(x: &[f64], y: &[f64], w: &[f64]) -> Option<f64> {
    let n = x.len();
    if n < 2 {
        return None;
    }
    let w_sum: f64 = w.iter().sum();
    if w_sum < 1e-15 {
        return None;
    }
    let mx: f64 = w.iter().zip(x.iter()).map(|(wi, xi)| wi * xi).sum::<f64>() / w_sum;
    let my: f64 = w.iter().zip(y.iter()).map(|(wi, yi)| wi * yi).sum::<f64>() / w_sum;
    let cov: f64 = w
        .iter()
        .zip(x.iter())
        .zip(y.iter())
        .map(|((wi, xi), yi)| wi * (xi - mx) * (yi - my))
        .sum::<f64>()
        / w_sum;
    let var_x: f64 = w.iter().zip(x.iter()).map(|(wi, xi)| wi * (xi - mx).powi(2)).sum::<f64>()
        / w_sum;
    let var_y: f64 = w.iter().zip(y.iter()).map(|(wi, yi)| wi * (yi - my).powi(2)).sum::<f64>()
        / w_sum;
    let denom = (var_x * var_y).sqrt();
    if denom < 1e-15 {
        return None;
    }
    Some(cov / denom)
}

struct SampleScores {
    pair_count: usize,
    mean_gold_dist: f64,
    mean_sample_dist: f64,
    weighted_log_pearson_r: f64,
    weighted_log_pearson_error: f64,
    spearman_r: f64,
    spearman_error: f64,
    weighted_kendall_tau: f64,
    weighted_kendall_error: f64,
}

fn compute_sample_scores(
    sample_id: &str,
    sample_genome: &str,
    all_samples: &[(String, String, f64)],
    sample_dists: &HashMap<(String, String), f64>,
    genome_dists: &HashMap<(String, String), f64>,
) -> SampleScores {
    let mut gold_vec: Vec<f64> = Vec::new();
    let mut pred_vec: Vec<f64> = Vec::new();

    for (other_id, other_genome, _) in all_samples {
        if other_id == sample_id {
            continue;
        }
        let gd = if sample_genome == other_genome.as_str() {
            Some(0.0_f64)
        } else {
            lookup_msa_dist(genome_dists, sample_genome, other_genome)
        };
        let sd = lookup_msa_dist(sample_dists, sample_id, other_id);
        if let (Some(gd), Some(sd)) = (gd, sd) {
            gold_vec.push(gd);
            pred_vec.push(sd);
        }
    }

    let pair_count = gold_vec.len();
    let nan_row = SampleScores {
        pair_count,
        mean_gold_dist: f64::NAN,
        mean_sample_dist: f64::NAN,
        weighted_log_pearson_r: f64::NAN,
        weighted_log_pearson_error: f64::NAN,
        spearman_r: f64::NAN,
        spearman_error: f64::NAN,
        weighted_kendall_tau: f64::NAN,
        weighted_kendall_error: f64::NAN,
    };

    if pair_count < MIN_PAIRS {
        return nan_row;
    }

    let mean_gold = gold_vec.iter().sum::<f64>() / pair_count as f64;
    let mean_sample = pred_vec.iter().sum::<f64>() / pair_count as f64;

    // Score 1: Weighted Log Pearson — scale-invariant, emphasises close pairs.
    let log_gold: Vec<f64> = gold_vec.iter().map(|&d| (d + MSA_EPSILON).ln()).collect();
    let log_pred: Vec<f64> = pred_vec.iter().map(|&d| (d + MSA_EPSILON).ln()).collect();
    let weights: Vec<f64> = gold_vec.iter().map(|&d| 1.0 / (d + MSA_EPSILON)).collect();
    let wlp_r = weighted_pearson_correlation(&log_gold, &log_pred, &weights).unwrap_or(f64::NAN);
    let wlp_error = 1.0 - wlp_r;

    // Score 2: Spearman — pure rank correlation, no absolute values.
    let rank_gold = rank_vector(&gold_vec);
    let rank_pred = rank_vector(&pred_vec);
    let sp_r = pearson_correlation(&rank_gold, &rank_pred).unwrap_or(f64::NAN);
    let sp_error = 1.0 - sp_r;

    // Score 3: Weighted Kendall τ — concordance of pairwise orderings,
    // weight = 1/min(d_gold_j, d_gold_k)+ε so near-neighbour pairs matter most.
    let mut concordant_w = 0.0_f64;
    let mut discordant_w = 0.0_f64;
    for i in 0..pair_count {
        for j in (i + 1)..pair_count {
            let sign_gold = (gold_vec[i] - gold_vec[j]).signum();
            if sign_gold == 0.0 {
                continue; // tied in gold: no ordering signal
            }
            let sign_pred = (pred_vec[i] - pred_vec[j]).signum();
            let w = 1.0 / (gold_vec[i].min(gold_vec[j]) + MSA_EPSILON);
            if sign_pred == sign_gold {
                concordant_w += w;
            } else {
                // discordant or tied in pred: prediction fails to order correctly
                discordant_w += w;
            }
        }
    }
    let total_w = concordant_w + discordant_w;
    let kendall_tau = if total_w < 1e-15 {
        f64::NAN
    } else {
        (concordant_w - discordant_w) / total_w
    };
    let kendall_error = 1.0 - kendall_tau;

    SampleScores {
        pair_count,
        mean_gold_dist: mean_gold,
        mean_sample_dist: mean_sample,
        weighted_log_pearson_r: wlp_r,
        weighted_log_pearson_error: wlp_error,
        spearman_r: sp_r,
        spearman_error: sp_error,
        weighted_kendall_tau: kendall_tau,
        weighted_kendall_error: kendall_error,
    }
}

fn compute_msa_error_rows(
    job: &StrainJob,
    msa_path: &Path,
    gold_msa_path: &Path,
    genome_groups: &BTreeMap<String, Vec<SampleInfo>>,
) -> Result<Vec<SampleErrorRow>, StrainError> {
    let msa_alignment = read_fasta_alignment(msa_path)?;
    let gold_alignment = read_fasta_alignment(gold_msa_path)?;

    info!(
        "  Sample MSA: {} sequences, {} positions",
        msa_alignment.sequences.len(),
        msa_alignment.length
    );
    info!(
        "  Gold MSA:   {} sequences, {} positions",
        gold_alignment.sequences.len(),
        gold_alignment.length
    );

    let sample_dists = compute_pairwise_msa_distances(&msa_alignment.sequences);
    let genome_dists = compute_pairwise_msa_distances(&gold_alignment.sequences);

    let all_samples: Vec<(String, String, f64)> = genome_groups
        .iter()
        .flat_map(|(genome, samples)| {
            samples
                .iter()
                .map(|s| (s.id.clone(), genome.clone(), s.coverage))
        })
        .collect();

    let mut error_rows = Vec::new();
    for (genome, samples) in genome_groups {
        for sample in samples {
            if !msa_alignment.sequences.contains_key(&sample.id) {
                warn!(
                    "Sample '{}' not found in MSA '{}', skipping.",
                    sample.id,
                    msa_path.display()
                );
                continue;
            }
            let sc = compute_sample_scores(
                &sample.id,
                genome,
                &all_samples,
                &sample_dists,
                &genome_dists,
            );
            error_rows.push(SampleErrorRow {
                id: job.id.clone(),
                species: job.species.clone(),
                sample: sample.id.clone(),
                genome: genome.clone(),
                coverage: sample.coverage,
                pair_count: sc.pair_count,
                mean_gold_dist: sc.mean_gold_dist,
                mean_sample_dist: sc.mean_sample_dist,
                weighted_log_pearson_r: sc.weighted_log_pearson_r,
                weighted_log_pearson_error: sc.weighted_log_pearson_error,
                spearman_r: sc.spearman_r,
                spearman_error: sc.spearman_error,
                weighted_kendall_tau: sc.weighted_kendall_tau,
                weighted_kendall_error: sc.weighted_kendall_error,
            });
        }
    }

    Ok(error_rows)
}

fn write_sample_error_output(
    rows: &[SampleErrorRow],
    output: &Path,
) -> Result<(), StrainError> {
    let mut ids = Vec::with_capacity(rows.len());
    let mut species = Vec::with_capacity(rows.len());
    let mut samples = Vec::with_capacity(rows.len());
    let mut genomes = Vec::with_capacity(rows.len());
    let mut coverages = Vec::with_capacity(rows.len());
    let mut pair_counts = Vec::with_capacity(rows.len());
    let mut mean_gold_dists = Vec::with_capacity(rows.len());
    let mut mean_sample_dists = Vec::with_capacity(rows.len());
    let mut wlp_rs = Vec::with_capacity(rows.len());
    let mut wlp_errors = Vec::with_capacity(rows.len());
    let mut sp_rs = Vec::with_capacity(rows.len());
    let mut sp_errors = Vec::with_capacity(rows.len());
    let mut kt_taus = Vec::with_capacity(rows.len());
    let mut kt_errors = Vec::with_capacity(rows.len());

    for row in rows {
        ids.push(row.id.clone());
        species.push(row.species.clone());
        samples.push(row.sample.clone());
        genomes.push(row.genome.clone());
        coverages.push(row.coverage);
        pair_counts.push(row.pair_count as u64);
        mean_gold_dists.push(row.mean_gold_dist);
        mean_sample_dists.push(row.mean_sample_dist);
        wlp_rs.push(row.weighted_log_pearson_r);
        wlp_errors.push(row.weighted_log_pearson_error);
        sp_rs.push(row.spearman_r);
        sp_errors.push(row.spearman_error);
        kt_taus.push(row.weighted_kendall_tau);
        kt_errors.push(row.weighted_kendall_error);
    }

    let mut df = DataFrame::new(vec![
        Series::new("ID".into(), ids),
        Series::new("Species".into(), species),
        Series::new("Sample".into(), samples),
        Series::new("Genome".into(), genomes),
        Series::new("Coverage".into(), coverages),
        Series::new("PairCount".into(), pair_counts),
        Series::new("MeanGoldDist".into(), mean_gold_dists),
        Series::new("MeanSampleDist".into(), mean_sample_dists),
        Series::new("WeightedLogPearsonR".into(), wlp_rs),
        Series::new("WeightedLogPearsonError".into(), wlp_errors),
        Series::new("SpearmanR".into(), sp_rs),
        Series::new("SpearmanError".into(), sp_errors),
        Series::new("WeightedKendallTau".into(), kt_taus),
        Series::new("WeightedKendallError".into(), kt_errors),
    ])
    .map_err(StrainError::Polars)?;

    let mut file = std::fs::File::create(output)?;
    CsvWriter::new(&mut file)
        .include_header(true)
        .with_separator(b'\t')
        .finish(&mut df)
        .map_err(StrainError::Polars)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{compute_pairwise_distances, find_lca};
    use phylotree::tree::Tree;

    fn sample_tree() -> Tree {
        Tree::from_newick("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);").unwrap()
    }

    #[test]
    fn find_lca_two_leaves() {
        let tree = sample_tree();
        let a = tree.get_by_name("A").unwrap().id;
        let b = tree.get_by_name("B").unwrap().id;
        let lca = find_lca(&tree, &[a, b]).unwrap();
        let lca_leaves = tree.get_subtree_leaves(&lca).unwrap();
        let lca_names: Vec<_> = lca_leaves
            .iter()
            .filter_map(|id| tree.get(id).unwrap().name.clone())
            .collect();
        assert!(lca_names.contains(&"A".to_string()));
        assert!(lca_names.contains(&"B".to_string()));
        assert!(!lca_names.contains(&"C".to_string()));
    }

    #[test]
    fn find_lca_all_leaves_is_root() {
        let tree = sample_tree();
        let a = tree.get_by_name("A").unwrap().id;
        let c = tree.get_by_name("C").unwrap().id;
        let lca = find_lca(&tree, &[a, c]).unwrap();
        let lca_leaves = tree.get_subtree_leaves(&lca).unwrap();
        assert_eq!(lca_leaves.len(), 4);
    }

    #[test]
    fn pairwise_distances_single_tip() {
        let tree = sample_tree();
        let a = tree.get_by_name("A").unwrap().id;
        let (min, max, mean, median) = compute_pairwise_distances(&tree, &[a]).unwrap();
        assert_eq!((min, max, mean, median), (0.0, 0.0, 0.0, 0.0));
    }

    #[test]
    fn pairwise_distances_two_tips() {
        let tree = sample_tree();
        let a = tree.get_by_name("A").unwrap().id;
        let b = tree.get_by_name("B").unwrap().id;
        let (min, max, mean, median) = compute_pairwise_distances(&tree, &[a, b]).unwrap();
        // A:0.1 + B:0.2 = 0.3 (one pair)
        assert!((min - 0.3).abs() < 1e-9);
        assert!((max - 0.3).abs() < 1e-9);
        assert!((mean - 0.3).abs() < 1e-9);
        assert!((median - 0.3).abs() < 1e-9);
    }

    #[test]
    fn monophyly_score_perfect() {
        let tree = sample_tree();
        let a = tree.get_by_name("A").unwrap().id;
        let b = tree.get_by_name("B").unwrap().id;
        let lca = find_lca(&tree, &[a, b]).unwrap();
        let lca_leaves = tree.get_subtree_leaves(&lca).unwrap();
        let score = 2.0 / lca_leaves.len() as f64;
        assert!((score - 1.0).abs() < 1e-9);
    }

    #[test]
    fn monophyly_score_imperfect() {
        let tree = sample_tree();
        let a = tree.get_by_name("A").unwrap().id;
        let c = tree.get_by_name("C").unwrap().id;
        let lca = find_lca(&tree, &[a, c]).unwrap();
        let lca_leaves = tree.get_subtree_leaves(&lca).unwrap();
        let score = 2.0 / lca_leaves.len() as f64;
        assert!((score - 0.5).abs() < 1e-9);
    }
}
