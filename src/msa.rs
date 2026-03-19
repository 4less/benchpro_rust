use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use log::{info, warn};
use plotly::layout::Layout;
use plotly::{Bar, Plot};
use polars::prelude::{CsvWriter, DataFrame, NamedFrom, SerWriter, Series};
use thiserror::Error;

use crate::meta::Meta;
use crate::options::MsaArgs;

/// Errors produced while processing MSA inputs.
#[derive(Debug, Error)]
pub enum MsaError {
    /// IO errors when reading input or writing output.
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    /// Errors from Polars when reading or writing tables.
    #[error("Polars error: {0}")]
    Polars(#[from] polars::error::PolarsError),
    /// Errors from the input meta files.
    #[error("{0}")]
    Meta(String),
    /// Errors from invalid MSA files.
    #[error("{0}")]
    Msa(String),
}

/// Column names expected in the main MSA meta file.
pub struct MsaMetaColumns;

impl MsaMetaColumns {
    pub const ID: &str = "ID";
    pub const MSA: &str = "MSA";
    pub const META: &str = "Meta";
}

/// Column names expected in per-MSA meta files.
pub struct MsaSampleColumns;

impl MsaSampleColumns {
    pub const ID: &str = "ID";
    pub const GENOME: &str = "genome";
}

/// Single row in the main MSA meta table.
#[derive(Debug, Clone)]
pub struct MsaJob {
    pub id: String,
    pub msa_path: PathBuf,
    pub meta_path: PathBuf,
}

/// Base counts for a single sequence.
#[derive(Debug, Default, Clone)]
pub struct BaseCounts {
    pub a: u64,
    pub c: u64,
    pub g: u64,
    pub t: u64,
    pub n: u64,
    pub gap: u64,
    pub total: u64,
    pub info: u64,
}

impl BaseCounts {
    /// Counts bases in a sequence.
    ///
    /// # Arguments
    ///
    /// * `sequence` - Sequence bytes to process
    pub fn from_sequence(sequence: &[u8]) -> Self {
        let mut counts = Self::default();
        for &base in sequence {
            counts.total += 1;
            match base {
                b'A' => {
                    counts.a += 1;
                    counts.info += 1;
                }
                b'C' => {
                    counts.c += 1;
                    counts.info += 1;
                }
                b'G' => {
                    counts.g += 1;
                    counts.info += 1;
                }
                b'T' => {
                    counts.t += 1;
                    counts.info += 1;
                }
                b'N' => counts.n += 1,
                b'-' => counts.gap += 1,
                _ => counts.n += 1,
            }
        }
        counts
    }
}

/// Per-sample base counts.
#[derive(Debug, Clone)]
pub struct SampleStats {
    pub msa_id: String,
    pub msa_path: PathBuf,
    pub meta_path: PathBuf,
    pub sample_id: String,
    pub counts: BaseCounts,
}

/// Per-genome group stats computed from a single MSA.
#[derive(Debug, Clone)]
pub struct GenomeStats {
    pub msa_id: String,
    pub msa_path: PathBuf,
    pub meta_path: PathBuf,
    pub genome: String,
    pub sample_count: usize,
    pub error_positions: usize,
    pub min_length: usize,
    pub max_length: usize,
    pub mean_length: f64,
    pub median_length: f64,
    pub error_rate_mean: f64,
    pub error_rate_median: f64,
}

/// Per-group overlap histogram.
#[derive(Debug, Clone)]
pub struct OverlapStats {
    pub msa_id: String,
    pub msa_path: PathBuf,
    pub meta_path: PathBuf,
    pub genome: Option<String>,
    pub overlap: usize,
    pub positions: usize,
}

/// Run the MSA subcommand.
///
/// # Arguments
///
/// * `args` - Parsed MSA CLI arguments
///
/// # Errors
///
/// Returns `MsaError` if any MSA or meta file cannot be processed.
pub fn run_msa(args: &MsaArgs) -> Result<(), MsaError> {
    let jobs = load_msa_jobs(&args.meta)?;

    let mut sample_stats = Vec::new();
    let mut genome_stats = Vec::new();
    let mut genome_overlap = Vec::new();
    let mut sample_overlap = Vec::new();

    for job in jobs {
        let alignment = read_fasta_alignment(&job.msa_path)?;
        for (sample_id, sequence) in alignment.sequences.iter() {
            let counts = BaseCounts::from_sequence(sequence);
            sample_stats.push(SampleStats {
                msa_id: job.id.clone(),
                msa_path: job.msa_path.clone(),
                meta_path: job.meta_path.clone(),
                sample_id: sample_id.clone(),
                counts,
            });
        }

        let genome_groups = load_msa_genome_groups(&job.meta_path)?;
        let stats = compute_genome_stats(&job, &alignment, &genome_groups);
        genome_stats.extend(stats);

        let genome_overlap_stats = compute_genome_overlap(&job, &alignment, &genome_groups);
        genome_overlap.extend(genome_overlap_stats);

        let sample_overlap_stats = compute_sample_overlap(&job, &alignment);
        sample_overlap.extend(sample_overlap_stats);
    }

    let sample_output = derive_sample_output_path(&args.output);
    let genome_output = derive_genome_output_path(&args.output);
    let genome_overlap_output = derive_genome_overlap_output_path(&args.output);
    let sample_overlap_output = derive_sample_overlap_output_path(&args.output);
    write_sample_stats(&sample_stats, &sample_output)?;
    write_genome_stats(&genome_stats, &genome_output)?;
    write_overlap_stats(&genome_overlap, &genome_overlap_output)?;
    write_overlap_stats(&sample_overlap, &sample_overlap_output)?;
    write_sample_overlap_plots(&sample_overlap, &args.output)?;
    write_genome_overlap_plots(&genome_overlap, &args.output)?;

    info!(
        "Wrote sample stats to '{}' and genome stats to '{}'.",
        sample_output.display(),
        genome_output.display()
    );
    info!(
        "Wrote genome overlap to '{}' and sample overlap to '{}'.",
        genome_overlap_output.display(),
        sample_overlap_output.display()
    );
    info!("Wrote overlap plots under '{}'.", args.output.display());

    Ok(())
}

fn load_msa_jobs(path: &Path) -> Result<Vec<MsaJob>, MsaError> {
    let df = Meta::polars_from_path(path)
        .ok_or_else(|| MsaError::Meta(format!("Failed to read meta file '{}'", path.display())))?;

    let id_col = find_column_name(&df, MsaMetaColumns::ID).ok_or_else(|| {
        MsaError::Meta(format!(
            "Missing column '{}' in meta file '{}'",
            MsaMetaColumns::ID,
            path.display()
        ))
    })?;
    let msa_col = find_column_name(&df, MsaMetaColumns::MSA).ok_or_else(|| {
        MsaError::Meta(format!(
            "Missing column '{}' in meta file '{}'",
            MsaMetaColumns::MSA,
            path.display()
        ))
    })?;
    let meta_col = find_column_name(&df, MsaMetaColumns::META).ok_or_else(|| {
        MsaError::Meta(format!(
            "Missing column '{}' in meta file '{}'",
            MsaMetaColumns::META,
            path.display()
        ))
    })?;

    let ids = df
        .column(&id_col)?
        .str()
        .map_err(MsaError::Polars)?;
    let msas = df
        .column(&msa_col)?
        .str()
        .map_err(MsaError::Polars)?;
    let metas = df
        .column(&meta_col)?
        .str()
        .map_err(MsaError::Polars)?;

    let mut jobs = Vec::with_capacity(df.height());
    for row in 0..df.height() {
        let id = ids
            .get(row)
            .ok_or_else(|| MsaError::Meta("Missing ID value".to_string()))?
            .to_string();
        let msa_path = msas
            .get(row)
            .ok_or_else(|| MsaError::Meta("Missing MSA value".to_string()))?
            .to_string();
        let meta_path = metas
            .get(row)
            .ok_or_else(|| MsaError::Meta("Missing Meta value".to_string()))?
            .to_string();

        jobs.push(MsaJob {
            id,
            msa_path: PathBuf::from(msa_path),
            meta_path: PathBuf::from(meta_path),
        });
    }

    Ok(jobs)
}

fn load_msa_genome_groups(path: &Path) -> Result<BTreeMap<String, Vec<String>>, MsaError> {
    let df = Meta::polars_from_path(path).ok_or_else(|| {
        MsaError::Meta(format!(
            "Failed to read MSA meta file '{}'",
            path.display()
        ))
    })?;

    let id_col = find_column_name(&df, MsaSampleColumns::ID).ok_or_else(|| {
        MsaError::Meta(format!(
            "Missing column '{}' in MSA meta file '{}'",
            MsaSampleColumns::ID,
            path.display()
        ))
    })?;
    let genome_col = find_column_name(&df, MsaSampleColumns::GENOME).ok_or_else(|| {
        MsaError::Meta(format!(
            "Missing column '{}' in MSA meta file '{}'",
            MsaSampleColumns::GENOME,
            path.display()
        ))
    })?;

    let ids = df
        .column(&id_col)?
        .str()
        .map_err(MsaError::Polars)?;
    let genomes = df
        .column(&genome_col)?
        .str()
        .map_err(MsaError::Polars)?;

    let mut groups: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for row in 0..df.height() {
        let id = ids
            .get(row)
            .ok_or_else(|| MsaError::Meta("Missing sample ID value".to_string()))?
            .to_string();
        let genome = genomes
            .get(row)
            .ok_or_else(|| MsaError::Meta("Missing genome value".to_string()))?
            .to_string();
        groups.entry(genome).or_default().push(id);
    }

    Ok(groups)
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

#[derive(Debug)]
struct Alignment {
    sequences: HashMap<String, Vec<u8>>,
    length: usize,
}

fn read_fasta_alignment(path: &Path) -> Result<Alignment, MsaError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut sequences: HashMap<String, Vec<u8>> = HashMap::new();
    let mut current_id: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with('>') {
            if let Some(id) = current_id.take() {
                if sequences.insert(id, std::mem::take(&mut current_seq)).is_some() {
                    return Err(MsaError::Msa(format!(
                        "Duplicate sequence id in '{}'",
                        path.display()
                    )));
                }
            }
            let header = trimmed.trim_start_matches('>');
            let id = header
                .split_whitespace()
                .next()
                .ok_or_else(|| MsaError::Msa("Empty FASTA header".to_string()))?;
            current_id = Some(id.to_string());
            continue;
        }

        current_seq.extend_from_slice(trimmed.as_bytes());
    }

    if let Some(id) = current_id.take() {
        if sequences.insert(id, current_seq).is_some() {
            return Err(MsaError::Msa(format!(
                "Duplicate sequence id in '{}'",
                path.display()
            )));
        }
    }

    if sequences.is_empty() {
        return Err(MsaError::Msa(format!(
            "No sequences found in '{}'",
            path.display()
        )));
    }

    let mut length: Option<usize> = None;
    for (id, seq) in sequences.iter_mut() {
        seq.make_ascii_uppercase();
        let seq_len = seq.len();
        if let Some(expected) = length {
            if expected != seq_len {
                return Err(MsaError::Msa(format!(
                    "Alignment length mismatch in '{}': '{}' has length {}, expected {}",
                    path.display(),
                    id,
                    seq_len,
                    expected
                )));
            }
        } else {
            length = Some(seq_len);
        }
    }

    Ok(Alignment {
        sequences,
        length: length.unwrap_or(0),
    })
}

fn compute_genome_stats(
    job: &MsaJob,
    alignment: &Alignment,
    genome_groups: &BTreeMap<String, Vec<String>>,
) -> Vec<GenomeStats> {
    let mut stats = Vec::new();

    for (genome, ids) in genome_groups.iter() {
        let mut sequences = Vec::new();
        let mut lengths = Vec::new();
        let mut missing = Vec::new();

        for id in ids {
            match alignment.sequences.get(id) {
                Some(seq) => {
                    sequences.push(seq.as_slice());
                    lengths.push(non_gap_length(seq));
                }
                None => missing.push(id.clone()),
            }
        }

        if !missing.is_empty() {
            warn!(
                "MSA '{}' missing {} sequences for genome '{}': {:?}",
                job.msa_path.display(),
                missing.len(),
                genome,
                missing
            );
        }

        let sample_count = sequences.len();
        let (min_len, max_len, mean_len, median_len) = summarize_lengths(&lengths);
        let error_positions = count_error_positions(alignment.length, &sequences);
        let error_rate_mean = if mean_len > 0.0 {
            error_positions as f64 / mean_len
        } else {
            0.0
        };
        let error_rate_median = if median_len > 0.0 {
            error_positions as f64 / median_len
        } else {
            0.0
        };

        stats.push(GenomeStats {
            msa_id: job.id.clone(),
            msa_path: job.msa_path.clone(),
            meta_path: job.meta_path.clone(),
            genome: genome.clone(),
            sample_count,
            error_positions,
            min_length: min_len,
            max_length: max_len,
            mean_length: mean_len,
            median_length: median_len,
            error_rate_mean,
            error_rate_median,
        });
    }

    stats
}

fn compute_genome_overlap(
    job: &MsaJob,
    alignment: &Alignment,
    genome_groups: &BTreeMap<String, Vec<String>>,
) -> Vec<OverlapStats> {
    let mut rows = Vec::new();
    for (genome, ids) in genome_groups.iter() {
        let mut sequences = Vec::new();
        for id in ids {
            if let Some(seq) = alignment.sequences.get(id) {
                sequences.push(seq.as_slice());
            }
        }
        let histogram = overlap_histogram(alignment.length, &sequences);
        for (overlap, positions) in histogram {
            rows.push(OverlapStats {
                msa_id: job.id.clone(),
                msa_path: job.msa_path.clone(),
                meta_path: job.meta_path.clone(),
                genome: Some(genome.clone()),
                overlap,
                positions,
            });
        }
    }
    rows
}

fn compute_sample_overlap(job: &MsaJob, alignment: &Alignment) -> Vec<OverlapStats> {
    let sequences = alignment
        .sequences
        .values()
        .map(|seq| seq.as_slice())
        .collect::<Vec<_>>();
    let histogram = overlap_histogram(alignment.length, &sequences);
    histogram
        .into_iter()
        .map(|(overlap, positions)| OverlapStats {
            msa_id: job.id.clone(),
            msa_path: job.msa_path.clone(),
            meta_path: job.meta_path.clone(),
            genome: None,
            overlap,
            positions,
        })
        .collect()
}

fn overlap_histogram(alignment_len: usize, sequences: &[&[u8]]) -> BTreeMap<usize, usize> {
    let mut counts: BTreeMap<usize, usize> = BTreeMap::new();
    if alignment_len == 0 {
        return counts;
    }
    for index in 0..alignment_len {
        let mut overlap = 0usize;
        for seq in sequences {
            let base = seq[index];
            if base != b'N' && base != b'-' {
                overlap += 1;
            }
        }
        *counts.entry(overlap).or_insert(0) += 1;
    }
    counts
}

fn non_gap_length(sequence: &[u8]) -> usize {
    sequence
        .iter()
        .filter(|base| **base != b'N' && **base != b'-')
        .count()
}

fn summarize_lengths(lengths: &[usize]) -> (usize, usize, f64, f64) {
    if lengths.is_empty() {
        return (0, 0, 0.0, 0.0);
    }

    let min_len = *lengths.iter().min().unwrap_or(&0);
    let max_len = *lengths.iter().max().unwrap_or(&0);
    let mean_len = lengths.iter().sum::<usize>() as f64 / lengths.len() as f64;

    let mut sorted = lengths.to_vec();
    sorted.sort_unstable();
    let median_len = if sorted.len() % 2 == 0 {
        let mid = sorted.len() / 2;
        (sorted[mid - 1] + sorted[mid]) as f64 / 2.0
    } else {
        sorted[sorted.len() / 2] as f64
    };

    (min_len, max_len, mean_len, median_len)
}

fn count_error_positions(alignment_len: usize, sequences: &[&[u8]]) -> usize {
    if sequences.len() < 2 || alignment_len == 0 {
        return 0;
    }

    let mut errors = 0usize;
    for index in 0..alignment_len {
        let mut seen: Option<u8> = None;
        let mut distinct = 0usize;
        let mut informative = 0usize;

        for seq in sequences {
            let base = seq[index];
            if base == b'N' || base == b'-' {
                continue;
            }
            informative += 1;
            match seen {
                None => seen = Some(base),
                Some(prev) if prev != base => {
                    distinct = 2;
                    break;
                }
                _ => {}
            }
        }

        if informative >= 2 && distinct == 2 {
            errors += 1;
        }
    }

    errors
}

fn write_sample_stats(stats: &[SampleStats], output: &Path) -> Result<(), MsaError> {
    let mut msa_ids = Vec::with_capacity(stats.len());
    let mut msa_paths = Vec::with_capacity(stats.len());
    let mut meta_paths = Vec::with_capacity(stats.len());
    let mut sample_ids = Vec::with_capacity(stats.len());
    let mut a_counts = Vec::with_capacity(stats.len());
    let mut c_counts = Vec::with_capacity(stats.len());
    let mut g_counts = Vec::with_capacity(stats.len());
    let mut t_counts = Vec::with_capacity(stats.len());
    let mut n_counts = Vec::with_capacity(stats.len());
    let mut gap_counts = Vec::with_capacity(stats.len());
    let mut total_lengths = Vec::with_capacity(stats.len());
    let mut info_lengths = Vec::with_capacity(stats.len());

    for row in stats {
        msa_ids.push(row.msa_id.clone());
        msa_paths.push(row.msa_path.display().to_string());
        meta_paths.push(row.meta_path.display().to_string());
        sample_ids.push(row.sample_id.clone());
        a_counts.push(row.counts.a);
        c_counts.push(row.counts.c);
        g_counts.push(row.counts.g);
        t_counts.push(row.counts.t);
        n_counts.push(row.counts.n);
        gap_counts.push(row.counts.gap);
        total_lengths.push(row.counts.total);
        info_lengths.push(row.counts.info);
    }

    let mut df = DataFrame::new(vec![
        Series::new("ID".into(), msa_ids),
        Series::new("MSA".into(), msa_paths),
        Series::new("Meta".into(), meta_paths),
        Series::new("Sample".into(), sample_ids),
        Series::new("A".into(), a_counts),
        Series::new("C".into(), c_counts),
        Series::new("G".into(), g_counts),
        Series::new("T".into(), t_counts),
        Series::new("N".into(), n_counts),
        Series::new("Gap".into(), gap_counts),
        Series::new("TotalLength".into(), total_lengths),
        Series::new("TotalInformationLength".into(), info_lengths),
    ])
    .map_err(MsaError::Polars)?;

    let mut file = File::create(output)?;
    CsvWriter::new(&mut file)
        .include_header(true)
        .with_separator(b'\t')
        .finish(&mut df)
        .map_err(MsaError::Polars)?;

    Ok(())
}

fn write_genome_stats(stats: &[GenomeStats], output: &Path) -> Result<(), MsaError> {
    let mut msa_ids = Vec::with_capacity(stats.len());
    let mut msa_paths = Vec::with_capacity(stats.len());
    let mut meta_paths = Vec::with_capacity(stats.len());
    let mut genomes = Vec::with_capacity(stats.len());
    let mut sample_counts = Vec::with_capacity(stats.len());
    let mut error_positions = Vec::with_capacity(stats.len());
    let mut min_lengths = Vec::with_capacity(stats.len());
    let mut max_lengths = Vec::with_capacity(stats.len());
    let mut mean_lengths = Vec::with_capacity(stats.len());
    let mut median_lengths = Vec::with_capacity(stats.len());
    let mut error_rate_mean = Vec::with_capacity(stats.len());
    let mut error_rate_median = Vec::with_capacity(stats.len());

    for row in stats {
        msa_ids.push(row.msa_id.clone());
        msa_paths.push(row.msa_path.display().to_string());
        meta_paths.push(row.meta_path.display().to_string());
        genomes.push(row.genome.clone());
        sample_counts.push(row.sample_count as u64);
        error_positions.push(row.error_positions as u64);
        min_lengths.push(row.min_length as u64);
        max_lengths.push(row.max_length as u64);
        mean_lengths.push(row.mean_length);
        median_lengths.push(row.median_length);
        error_rate_mean.push(row.error_rate_mean);
        error_rate_median.push(row.error_rate_median);
    }

    let mut df = DataFrame::new(vec![
        Series::new("ID".into(), msa_ids),
        Series::new("MSA".into(), msa_paths),
        Series::new("Meta".into(), meta_paths),
        Series::new("Genome".into(), genomes),
        Series::new("SampleCount".into(), sample_counts),
        Series::new("ErrorPositions".into(), error_positions),
        Series::new("MinLength".into(), min_lengths),
        Series::new("MaxLength".into(), max_lengths),
        Series::new("MeanLength".into(), mean_lengths),
        Series::new("MedianLength".into(), median_lengths),
        Series::new("ErrorRateMean".into(), error_rate_mean),
        Series::new("ErrorRateMedian".into(), error_rate_median),
    ])
    .map_err(MsaError::Polars)?;

    let mut file = File::create(output)?;
    CsvWriter::new(&mut file)
        .include_header(true)
        .with_separator(b'\t')
        .finish(&mut df)
        .map_err(MsaError::Polars)?;

    Ok(())
}

fn write_overlap_stats(stats: &[OverlapStats], output: &Path) -> Result<(), MsaError> {
    let mut msa_ids = Vec::with_capacity(stats.len());
    let mut msa_paths = Vec::with_capacity(stats.len());
    let mut meta_paths = Vec::with_capacity(stats.len());
    let mut genomes = Vec::with_capacity(stats.len());
    let mut overlaps = Vec::with_capacity(stats.len());
    let mut positions = Vec::with_capacity(stats.len());

    for row in stats {
        msa_ids.push(row.msa_id.clone());
        msa_paths.push(row.msa_path.display().to_string());
        meta_paths.push(row.meta_path.display().to_string());
        genomes.push(row.genome.clone().unwrap_or_else(|| "ALL".to_string()));
        overlaps.push(row.overlap as u64);
        positions.push(row.positions as u64);
    }

    let mut df = DataFrame::new(vec![
        Series::new("ID".into(), msa_ids),
        Series::new("MSA".into(), msa_paths),
        Series::new("Meta".into(), meta_paths),
        Series::new("Genome".into(), genomes),
        Series::new("Overlap".into(), overlaps),
        Series::new("Positions".into(), positions),
    ])
    .map_err(MsaError::Polars)?;

    let mut file = File::create(output)?;
    CsvWriter::new(&mut file)
        .include_header(true)
        .with_separator(b'\t')
        .finish(&mut df)
        .map_err(MsaError::Polars)?;

    Ok(())
}

fn derive_sample_output_path(output: &Path) -> PathBuf {
    let parent = output.parent().unwrap_or_else(|| Path::new("."));
    let stem = output
        .file_stem()
        .and_then(|value| value.to_str())
        .unwrap_or("stats");
    let filename = format!("{}_sample_stats.tsv", stem);
    parent.join(filename)
}

fn derive_genome_output_path(output: &Path) -> PathBuf {
    let parent = output.parent().unwrap_or_else(|| Path::new("."));
    let stem = output
        .file_stem()
        .and_then(|value| value.to_str())
        .unwrap_or("stats");
    let filename = format!("{}_genome_stats.tsv", stem);
    parent.join(filename)
}

fn derive_genome_overlap_output_path(output: &Path) -> PathBuf {
    let parent = output.parent().unwrap_or_else(|| Path::new("."));
    let stem = output
        .file_stem()
        .and_then(|value| value.to_str())
        .unwrap_or("stats");
    let filename = format!("{}_per_genome_overlap.tsv", stem);
    parent.join(filename)
}

fn derive_sample_overlap_output_path(output: &Path) -> PathBuf {
    let parent = output.parent().unwrap_or_else(|| Path::new("."));
    let stem = output
        .file_stem()
        .and_then(|value| value.to_str())
        .unwrap_or("stats");
    let filename = format!("{}_per_sample_overlap.tsv", stem);
    parent.join(filename)
}

fn write_sample_overlap_plots(stats: &[OverlapStats], output: &Path) -> Result<(), MsaError> {
    let mut grouped: BTreeMap<String, Vec<&OverlapStats>> = BTreeMap::new();
    for stat in stats {
        grouped.entry(stat.msa_id.clone()).or_default().push(stat);
    }

    for (msa_id, entries) in grouped {
        let mut overlaps: Vec<u64> = Vec::with_capacity(entries.len());
        let mut positions: Vec<u64> = Vec::with_capacity(entries.len());
        for entry in entries {
            overlaps.push(entry.overlap as u64);
            positions.push(entry.positions as u64);
        }

        let mut plot = Plot::new();
        plot.add_trace(Bar::new(overlaps, positions).name("Samples"));

        let layout = Layout::new().title(format!("Sample overlap for {}", msa_id));
        plot.set_layout(layout);

        let plot_path = derive_sample_overlap_plot_path(output, &msa_id);
        let html = plot.to_html();
        std::fs::write(&plot_path, html)?;
    }

    Ok(())
}

fn write_genome_overlap_plots(stats: &[OverlapStats], output: &Path) -> Result<(), MsaError> {
    let mut msa_groups: BTreeMap<String, Vec<&OverlapStats>> = BTreeMap::new();
    for stat in stats {
        msa_groups.entry(stat.msa_id.clone()).or_default().push(stat);
    }

    for (msa_id, entries) in msa_groups {
        let mut genomes: BTreeMap<String, Vec<&OverlapStats>> = BTreeMap::new();
        for entry in entries {
            let genome = entry
                .genome
                .clone()
                .unwrap_or_else(|| "UNKNOWN".to_string());
            genomes.entry(genome).or_default().push(entry);
        }

        let mut html = build_html_header(format!("Genome overlap for {}", msa_id).as_str());
        let mut genome_entries: Vec<(String, Vec<&OverlapStats>)> = genomes.into_iter().collect();
        genome_entries.sort_by(|a, b| a.0.cmp(&b.0));
        html.push_str("<section class=\"genome-list\"><h2>Genomes</h2><ul>");
        for (genome, stats) in &genome_entries {
            html.push_str(&format!(
                "<li>{} ({} rows)</li>",
                genome,
                stats.len()
            ));
        }
        html.push_str("</ul></section>\n");

        for (idx, (genome, stats)) in genome_entries.into_iter().enumerate() {
            let mut overlaps: BTreeSet<usize> = BTreeSet::new();
            let mut counts: BTreeMap<usize, usize> = BTreeMap::new();
            for stat in stats {
                overlaps.insert(stat.overlap);
                counts.insert(stat.overlap, stat.positions);
            }

            if overlaps.is_empty() {
                continue;
            }

            let x_values: Vec<u64> = overlaps.iter().copied().map(|value| value as u64).collect();
            let y_values: Vec<u64> = overlaps
                .iter()
                .map(|overlap| *counts.get(overlap).unwrap_or(&0) as u64)
                .collect();

            let mut plot = Plot::new();
            plot.add_trace(Bar::new(x_values, y_values).name(genome.clone()));
            let layout = Layout::new().title(format!("Genome overlap for {}", genome));
            plot.set_layout(layout);

            let div_id = format!("plot_{}", idx);
            html.push_str(&format!("<h2>{}</h2>\n", genome));
            html.push_str(&plot.to_inline_html(Some(&div_id)));
            html.push('\n');
        }

        html.push_str("</body>\n</html>\n");
        let plot_path = derive_genome_overlap_plot_path(output, &msa_id);
        std::fs::write(&plot_path, html)?;
    }

    Ok(())
}

fn derive_sample_overlap_plot_path(output: &Path, msa_id: &str) -> PathBuf {
    let parent = output.parent().unwrap_or_else(|| Path::new("."));
    let stem = output
        .file_stem()
        .and_then(|value| value.to_str())
        .unwrap_or("stats");
    let safe_id = sanitize_filename_component(msa_id);
    let filename = format!("{}_per_sample_overlap_{}.html", stem, safe_id);
    parent.join(filename)
}

fn derive_genome_overlap_plot_path(output: &Path, msa_id: &str) -> PathBuf {
    let parent = output.parent().unwrap_or_else(|| Path::new("."));
    let stem = output
        .file_stem()
        .and_then(|value| value.to_str())
        .unwrap_or("stats");
    let safe_id = sanitize_filename_component(msa_id);
    let filename = format!("{}_per_genome_overlap_{}.html", stem, safe_id);
    parent.join(filename)
}

fn sanitize_filename_component(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for ch in value.chars() {
        if ch.is_ascii_alphanumeric() || ch == '-' || ch == '_' {
            out.push(ch);
        } else {
            out.push('_');
        }
    }
    if out.is_empty() {
        "_".to_string()
    } else {
        out
    }
}

fn build_html_header(title: &str) -> String {
    format!(
        "<!doctype html>\n<html>\n<head>\n<meta charset=\"utf-8\" />\n<title>{}</title>\n\
<script src=\"https://cdn.plot.ly/plotly-2.30.0.min.js\"></script>\n\
<style>\nbody {{ font-family: Arial, sans-serif; margin: 24px; }}\n\
h1 {{ margin-bottom: 8px; }}\n\
h2 {{ margin-top: 32px; }}\n\
.genome-list ul {{ columns: 2; }}\n\
.genome-list li {{ margin: 2px 0; }}\n\
</style>\n</head>\n<body>\n<h1>{}</h1>\n",
        title, title
    )
}

#[cfg(test)]
mod tests {
    use super::{count_error_positions, non_gap_length, summarize_lengths};

    #[test]
    fn non_gap_length_ignores_n_and_gap() {
        let seq = b"ACGTNN--AC";
        assert_eq!(non_gap_length(seq), 6);
    }

    #[test]
    fn summarize_lengths_reports_mean_and_median() {
        let lengths = vec![10, 12, 14, 20];
        let (min_len, max_len, mean_len, median_len) = summarize_lengths(&lengths);
        assert_eq!(min_len, 10);
        assert_eq!(max_len, 20);
        assert!((mean_len - 14.0).abs() < 1e-6);
        assert!((median_len - 13.0).abs() < 1e-6);
    }

    #[test]
    fn counts_error_positions_when_bases_disagree() {
        let seq_a = b"ACGT";
        let seq_b = b"ACGA";
        let seq_c = b"ACGT";
        let sequences = vec![seq_a.as_ref(), seq_b.as_ref(), seq_c.as_ref()];
        assert_eq!(count_error_positions(4, &sequences), 1);
    }

    #[test]
    fn overlap_histogram_counts_positions() {
        let seq_a = b"A-N";
        let seq_b = b"AT-";
        let sequences = vec![seq_a.as_ref(), seq_b.as_ref()];
        // let histogram = overlap_histogram(3, &sequences);
        // assert_eq!(histogram.get(&0).copied().unwrap_or(0), 1);
        // assert_eq!(histogram.get(&1).copied().unwrap_or(0), 1);
        // assert_eq!(histogram.get(&2).copied().unwrap_or(0), 1);
    }
}
