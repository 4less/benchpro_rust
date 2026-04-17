use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::path::{Path, PathBuf};

use log::{info, warn};
use phylotree::tree::{NodeId, Tree};
use polars::prelude::{CsvWriter, DataFrame, NamedFrom, SerWriter, Series};
use thiserror::Error;

use crate::meta::Meta;
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
}

struct StrainMetaColumns;
impl StrainMetaColumns {
    const ID: &'static str = "ID";
    const SPECIES: &'static str = "Species";
    const TREE: &'static str = "Tree";
    const META: &'static str = "Meta";
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
    }

    let output = PathBuf::from(format!("{}.monophyly.tsv", args.outprefix));
    write_monophyly_output(&rows, &output)?;
    info!("Wrote monophyly stats to '{}'.", output.display());

    let tip_output = PathBuf::from(format!("{}.monophyly_tips.tsv", args.outprefix));
    write_tip_output(&tip_rows, &tip_output)?;
    info!("Wrote per-tip monophyly stats to '{}'.", tip_output.display());

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

        jobs.push(StrainJob {
            id,
            species,
            tree_path: PathBuf::from(tree_path),
            meta_path: PathBuf::from(meta_path),
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
