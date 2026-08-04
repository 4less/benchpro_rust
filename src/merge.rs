use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::BufRead;
use std::path::{Path, PathBuf};

use indicatif::{ProgressBar, ProgressStyle};
use itertools::Itertools;
use log::info;
use polars::prelude::{
    ChunkAgg, CsvWriter, DataFrame, Float64Chunked, IntoSeries, NamedFrom, SerWriter, Series,
};
use rayon::prelude::*;
use regex::Regex;
use thiserror::Error;

use crate::common::{TaxonomicRank, Taxonomy};
use crate::options::MergeArgs;
use crate::profile::{Entry, Profile, ProfileWrapper};
use crate::profile_handler::ProfileHandler;

/// Errors produced while merging profile files.
#[derive(Debug, Error)]
pub enum MergeError {
    /// IO errors when reading input or writing output.
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    /// Errors produced while parsing profile files.
    #[error("Failed to parse profile '{path}': {errors:?}")]
    ProfileParse {
        /// Input file that failed to parse.
        path: PathBuf,
        /// Collected error messages for attempted parsers.
        errors: Vec<String>,
    },
    /// Errors from Polars when building or writing the output matrix.
    #[error("Polars error: {0}")]
    Polars(#[from] polars::error::PolarsError),
}

type TaxonKey = String;

/// Merge profiles into an abundance matrix and write it to disk.
///
/// # Arguments
///
/// * `inputs` - Profile file paths
/// * `output` - Output TSV path
///
/// # Returns
///
/// Column sums and a map of sample names to their source paths.
///
/// # Errors
///
/// Returns `MergeError::ProfileParse` when a profile cannot be parsed.
/// Returns `MergeError::Io` for file I/O failures.
/// Returns `MergeError::Polars` for DataFrame construction or write errors.
pub fn merge_profiles(
    inputs: &[PathBuf],
    output: &Path,
    target_rank: &TaxonomicRank,
) -> Result<(HashMap<String, f64>, HashMap<String, PathBuf>), MergeError> {
    let (matrix, sample_paths) = build_abundance_matrix(inputs, target_rank, &[])?;
    let sums = write_abundance_matrix(&matrix, output)?;
    Ok((sums, sample_paths))
}

fn merge_profiles_with_regexes(
    inputs: &[PathBuf],
    output: &Path,
    target_rank: &TaxonomicRank,
    name_regexes: &[Regex],
) -> Result<(HashMap<String, f64>, HashMap<String, PathBuf>), MergeError> {
    let (matrix, sample_paths) = build_abundance_matrix(inputs, target_rank, name_regexes)?;
    let sums = write_abundance_matrix(&matrix, output)?;
    Ok((sums, sample_paths))
}

/// Run the merge subcommand.
///
/// # Arguments
///
/// * `args` - Parsed merge CLI arguments
///
/// # Errors
///
/// Propagates errors from `merge_profiles`.
pub fn run_merge(args: &MergeArgs) -> Result<(), MergeError> {
    let target_rank = parse_rank(&args.rank)?;
    let name_regexes = parse_name_regexes(&args.sample_regex)?;
    if args.test_sample_regex {
        test_sample_regex(&args.input, &name_regexes)?;
        return Ok(());
    }
    let (sums, sample_paths) =
        merge_profiles_with_regexes(&args.input, &args.output, &target_rank, &name_regexes)?;
    print_column_sums(&sums, &sample_paths);
    Ok(())
}

fn build_abundance_matrix(
    inputs: &[PathBuf],
    target_rank: &TaxonomicRank,
    name_regexes: &[Regex],
) -> Result<(DataFrame, HashMap<String, PathBuf>), MergeError> {
    let progress = ProgressBar::new(inputs.len() as u64);
    let style =
        ProgressStyle::with_template("{msg} {bar:40.cyan/blue} {pos}/{len} ({elapsed_precise})")
            .unwrap();
    progress.set_style(style);
    progress.set_message("Merging profiles");

    let profiles: Vec<(usize, PathBuf, ProfileWrapper)> = inputs
        .par_iter()
        .enumerate()
        .map(|(index, path)| {
            let profile = load_profile_auto(path)?;
            progress.inc(1);
            Ok((index, path.to_path_buf(), profile))
        })
        .collect::<Result<Vec<_>, MergeError>>()?;

    progress.finish_and_clear();

    let mut profiles = profiles;
    profiles.sort_by_key(|(index, _, _)| *index);

    let mut sample_names = Vec::with_capacity(profiles.len());
    let mut sample_paths: HashMap<String, PathBuf> = HashMap::default();
    let mut sample_maps: Vec<HashMap<TaxonKey, f64>> = Vec::with_capacity(profiles.len());
    let mut all_taxa: HashSet<TaxonKey> = HashSet::default();
    let mut name_counts: HashMap<String, usize> = HashMap::default();

    for (_, path, profile) in profiles {
        let sample_name = dedupe_sample_name(path.clone(), name_regexes, &mut name_counts)?;
        let mut abundance_map: HashMap<TaxonKey, f64> = HashMap::new();
        collect_abundances_from_wrapper(
            &profile,
            target_rank,
            &path,
            &mut abundance_map,
            &mut all_taxa,
        )?;

        sample_paths.insert(sample_name.clone(), path.clone());
        sample_names.push(sample_name);
        sample_maps.push(abundance_map);
    }

    let mut taxa_sorted = all_taxa.into_iter().collect_vec();
    taxa_sorted.sort();

    let mut columns = Vec::with_capacity(1 + sample_names.len());
    columns.push(Series::new(
        "Taxon".into(),
        taxa_sorted.iter().map(|t| t.as_str()).collect_vec(),
    ));

    for (sample_index, sample_name) in sample_names.iter().enumerate() {
        let values = taxa_sorted
            .iter()
            .map(|taxon| sample_maps[sample_index].get(taxon).copied().unwrap_or(0.0))
            .collect_vec();
        columns.push(Series::new(sample_name.as_str().into(), values));
    }

    let df = DataFrame::new(columns).map_err(MergeError::Polars)?;
    Ok((df, sample_paths))
}

fn write_abundance_matrix(
    df: &DataFrame,
    output: &Path,
) -> Result<HashMap<String, f64>, MergeError> {
    if let Some(parent) = output.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut file = File::create(output)?;
    let mut df = df.clone();
    CsvWriter::new(&mut file)
        .include_header(true)
        .with_separator(b'\t')
        .finish(&mut df)
        .map_err(MergeError::Polars)?;

    let mut sums = HashMap::default();
    let column_names = df
        .get_column_names()
        .iter()
        .skip(1)
        .map(|name| name.to_string())
        .collect_vec();
    for name in column_names {
        let series = df.column(&name).map_err(MergeError::Polars)?;
        let sum = series
            .f64()
            .map_err(MergeError::Polars)?
            .sum()
            .unwrap_or(0.0);

        if sum > 2.0 {
            let column_name = name.as_str();
            let mut normalized = series
                .f64()
                .map_err(MergeError::Polars)?
                .into_iter()
                .map(|value| value.map(|v| v / 100.0))
                .collect::<Float64Chunked>()
                .into_series();
            normalized.rename(column_name.into());
            df.with_column(normalized).map_err(MergeError::Polars)?;
            let normalized_sum = sum / 100.0;
            info!(
                "Normalized column '{}' by 1/100 (sum was {}).",
                column_name, sum
            );
            sums.insert(name, normalized_sum);
        } else {
            sums.insert(name, sum);
        }
    }

    Ok(sums)
}

fn load_profile_auto(path: &Path) -> Result<ProfileWrapper, MergeError> {
    ProfileHandler::load_profile_auto(path).ok_or_else(|| MergeError::ProfileParse {
        path: path.to_path_buf(),
        errors: vec!["Auto-detect profile loading failed in shared profile handler".to_owned()],
    })
}

/// Returns the resolved rank for an entry: explicit rank if set, otherwise
/// derived from the lowest node in the lineage.
fn resolved_rank_of<T: Taxonomy>(entry: &Entry<T>) -> TaxonomicRank {
    if entry.rank != TaxonomicRank::Unknown {
        return entry.rank.clone();
    }
    entry
        .lineage
        .as_ref()
        .and_then(|l| l.lowest().and_then(|t| t.rank.clone()))
        .unwrap_or(TaxonomicRank::Unknown)
}

fn collect_abundances_from_profile<T: Taxonomy>(
    profile: &Profile<T>,
    target_rank: &TaxonomicRank,
    path: &Path,
    map: &mut HashMap<TaxonKey, f64>,
    taxa: &mut HashSet<TaxonKey>,
) -> Result<(), MergeError> {
    // If the profile already has entries at exactly the target rank, restrict
    // collection to those entries only. This prevents double-counting in
    // profiles that list every rank independently (e.g. CAMI / sylph), where
    // a strain entry's lineage also contains the parent species and would
    // otherwise be aggregated into the species bucket alongside the explicit
    // species entry.
    //
    // If the profile has *no* entries at the target rank (e.g. a species-only
    // profile targeted at genus), we fall back to aggregating finer entries.
    let has_explicit_target_rank = profile
        .taxa
        .iter()
        .any(|e| resolved_rank_of(e) == *target_rank);

    for (index, entry) in profile.taxa.iter().enumerate() {
        if entry.abundance == 0.0 {
            continue;
        }
        if has_explicit_target_rank && resolved_rank_of(entry) != *target_rank {
            continue;
        }
        if let Some(key) = lineage_key(entry, target_rank, path, index)? {
            *map.entry(key.clone()).or_insert(0.0) += entry.abundance;
            taxa.insert(key);
        }
    }
    Ok(())
}

fn collect_abundances_from_wrapper(
    profile: &ProfileWrapper,
    target_rank: &TaxonomicRank,
    path: &Path,
    map: &mut HashMap<TaxonKey, f64>,
    taxa: &mut HashSet<TaxonKey>,
) -> Result<(), MergeError> {
    match profile {
        ProfileWrapper::GTDBProfile(inner) => {
            collect_abundances_from_profile(inner, target_rank, path, map, taxa)
        }
        ProfileWrapper::NCBIProfile(inner) => {
            collect_abundances_from_profile(inner, target_rank, path, map, taxa)
        }
        ProfileWrapper::ChocoPhlAnProfile(inner) => {
            collect_abundances_from_profile(inner, target_rank, path, map, taxa)
        }
        ProfileWrapper::CustomProfile(inner) => {
            collect_abundances_from_profile(inner, target_rank, path, map, taxa)
        }
    }
}

fn lineage_key<T: Taxonomy>(
    entry: &Entry<T>,
    target_rank: &TaxonomicRank,
    path: &Path,
    entry_index: usize,
) -> Result<Option<TaxonKey>, MergeError> {
    // Closure so file I/O only happens when an error is actually constructed.
    let row_label = || {
        read_culprit_row(path, entry_index)
            .map(|value| format!("row='{}'", value))
            .unwrap_or_else(|| "row=<unavailable>".to_string())
    };

    if entry.rank == TaxonomicRank::Unknown
        && entry.lineage.is_none()
        && entry
            .taxon_name
            .as_deref()
            .is_some_and(|name| name.trim().eq_ignore_ascii_case("UNCLASSIFIED"))
    {
        return Ok(Some(format!("{}UNCLASSIFIED", rank_prefix(target_rank))));
    }

    if entry.rank == TaxonomicRank::Unknown && entry.lineage.is_none() {
        return Err(MergeError::ProfileParse {
            path: path.to_path_buf(),
            errors: vec![format!(
                "Entry {} has no rank and no lineage; entry={}; {}",
                entry_index,
                entry_context(entry),
                row_label()
            )],
        });
    }

    let resolved_rank = if entry.rank == TaxonomicRank::Unknown {
        entry
            .lineage
            .as_ref()
            .and_then(|lineage| lineage.lowest().and_then(|taxon| taxon.rank.clone()))
    } else {
        Some(entry.rank.clone())
    };

    let resolved_rank = resolved_rank.ok_or_else(|| MergeError::ProfileParse {
        path: path.to_path_buf(),
        errors: vec![format!(
            "Unable to determine rank for entry {}; entry={}; {}",
            entry_index,
            entry_context(entry),
            row_label()
        )],
    })?;

    if let Some(lineage) = entry.lineage.as_ref() {
        if lineage.get(target_rank).is_some() {
            let delimiter = match T::get_enum() {
                crate::common::TaxonomyEnum::NCBI => "|",
                _ => ";",
            };
            let key = TaxonomicRank::all()
                .into_iter()
                .filter(|rank| *rank != TaxonomicRank::Unknown)
                .filter(|rank| *rank <= *target_rank)
                .filter_map(|rank| lineage.get(&rank).map(|taxon| taxon.name.clone()))
                .join(delimiter);
            if !key.is_empty() {
                return Ok(Some(key));
            }
        }
    }

    if &resolved_rank != target_rank {
        return Ok(None);
    }

    let lineage_name = entry
        .lineage
        .as_ref()
        .and_then(|lineage| lineage.get(&resolved_rank).map(|t| t.name.clone()))
        .or_else(|| {
            entry
                .lineage
                .as_ref()
                .and_then(|lineage| lineage.lowest().map(|t| t.name.clone()))
        });

    let mut name = entry
        .taxon_name
        .as_ref()
        .cloned()
        .or(lineage_name.clone())
        .ok_or_else(|| MergeError::ProfileParse {
            path: path.to_path_buf(),
            errors: vec![format!(
                "Unable to determine taxon name for entry {}; entry={}; {}",
                entry_index,
                entry_context(entry),
                row_label()
            )],
        })?;

    if name.chars().all(|c| c.is_ascii_digit()) {
        if let Some(preferred) = lineage_name {
            name = preferred;
        }
    }

    let prefix = rank_prefix(&resolved_rank);
    if prefix.is_empty() || name.starts_with(prefix) {
        Ok(Some(name))
    } else {
        Ok(Some(format!("{}{}", prefix, name)))
    }
}

fn rank_prefix(rank: &TaxonomicRank) -> &'static str {
    match rank {
        TaxonomicRank::Superkingdom => "k__",
        TaxonomicRank::Domain => "d__",
        TaxonomicRank::Phylum => "p__",
        TaxonomicRank::Class => "c__",
        TaxonomicRank::Order => "o__",
        TaxonomicRank::Family => "f__",
        TaxonomicRank::Genus => "g__",
        TaxonomicRank::Species => "s__",
        TaxonomicRank::Strain => "t__",
        TaxonomicRank::Unknown => "",
    }
}

fn dedupe_sample_name(
    path: PathBuf,
    name_regexes: &[Regex],
    counts: &mut HashMap<String, usize>,
) -> Result<String, MergeError> {
    let base = if !name_regexes.is_empty() {
        derive_sample_name(&path, name_regexes)?
    } else {
        path.file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("sample")
            .to_string()
    };
    let count = counts.entry(base.clone()).or_insert(0);
    let name = if *count == 0 {
        base
    } else {
        format!("{}_{}", base, *count + 1)
    };
    *count += 1;
    Ok(name)
}

fn parse_rank(rank: &str) -> Result<TaxonomicRank, MergeError> {
    TaxonomicRank::from_string(rank).ok_or_else(|| MergeError::ProfileParse {
        path: PathBuf::from("<rank>"),
        errors: vec![format!("Invalid rank '{}'", rank)],
    })
}

fn parse_name_regexes(patterns: &[String]) -> Result<Vec<Regex>, MergeError> {
    let mut regexes = Vec::with_capacity(patterns.len());
    for pattern in patterns {
        let regex = Regex::new(pattern).map_err(|e| MergeError::ProfileParse {
            path: PathBuf::from("<sample-regex>"),
            errors: vec![format!("Invalid --sample-regex '{}': {}", pattern, e)],
        })?;
        regexes.push(regex);
    }
    Ok(regexes)
}

fn derive_sample_name(path: &Path, regexes: &[Regex]) -> Result<String, MergeError> {
    let path_str = path.to_string_lossy();
    let (captures, used_pattern) = regexes
        .iter()
        .find_map(|regex| regex.captures(&path_str).map(|cap| (cap, regex)))
        .ok_or_else(|| MergeError::ProfileParse {
            path: path.to_path_buf(),
            errors: vec![format!("--sample-regex did not match path '{}'", path_str)],
        })?;
    let parts = (1..captures.len())
        .filter_map(|i| captures.get(i).map(|m| m.as_str()))
        .filter(|part| !part.is_empty())
        .collect_vec();
    if parts.is_empty() {
        return Err(MergeError::ProfileParse {
            path: path.to_path_buf(),
            errors: vec![format!(
                "--sample-regex had no capture groups for path '{}' (pattern '{}')",
                path_str,
                used_pattern.as_str()
            )],
        });
    }
    Ok(parts.join("_"))
}

fn test_sample_regex(inputs: &[PathBuf], regexes: &[Regex]) -> Result<(), MergeError> {
    let use_color = atty::is(atty::Stream::Stdout);
    let mut used = vec![false; regexes.len()];
    let mut unmatched: Vec<(String, String, String)> = Vec::new();
    for path in inputs {
        let (sample_name, pattern) = if regexes.is_empty() {
            let name = path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("sample")
                .to_string();
            (name, "<default>".to_string())
        } else {
            match derive_sample_name_with_pattern(path, regexes, &mut used) {
                Ok(value) => value,
                Err(_) => {
                    let name = path
                        .file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or("sample")
                        .to_string();
                    let pattern = format_none_pattern(use_color);
                    unmatched.push((name.clone(), path.display().to_string(), pattern.clone()));
                    (name, pattern)
                }
            }
        };
        println!("{}\t{}\t{}", sample_name, path.display(), pattern);
    }
    if !unmatched.is_empty() {
        println!("UNMATCHED");
        for (name, path, pattern) in unmatched {
            println!("{}\t{}\t{}", name, path, pattern);
        }
    }
    if !regexes.is_empty() {
        for (index, regex) in regexes.iter().enumerate() {
            if !used[index] {
                println!("UNUSED\t{}\t{}", index + 1, regex.as_str());
            }
        }
    }
    Ok(())
}

fn format_none_pattern(use_color: bool) -> String {
    if use_color {
        "\x1b[31mNONE\x1b[0m".to_string()
    } else {
        "NONE".to_string()
    }
}

fn derive_sample_name_with_pattern(
    path: &Path,
    regexes: &[Regex],
    used: &mut [bool],
) -> Result<(String, String), MergeError> {
    let path_str = path.to_string_lossy();
    let (captures, used_pattern, used_index) = regexes
        .iter()
        .enumerate()
        .find_map(|(index, regex)| regex.captures(&path_str).map(|cap| (cap, regex, index)))
        .ok_or_else(|| MergeError::ProfileParse {
            path: path.to_path_buf(),
            errors: vec![format!("--sample-regex did not match path '{}'", path_str)],
        })?;
    if let Some(slot) = used.get_mut(used_index) {
        *slot = true;
    }
    let parts = (1..captures.len())
        .filter_map(|i| captures.get(i).map(|m| m.as_str()))
        .filter(|part| !part.is_empty())
        .collect_vec();
    if parts.is_empty() {
        return Err(MergeError::ProfileParse {
            path: path.to_path_buf(),
            errors: vec![format!(
                "--sample-regex had no capture groups for path '{}' (pattern '{}')",
                path_str,
                used_pattern.as_str()
            )],
        });
    }
    Ok((parts.join("_"), used_pattern.as_str().to_string()))
}

fn print_column_sums(sums: &HashMap<String, f64>, sample_paths: &HashMap<String, PathBuf>) {
    let tolerance_ratio = 0.05_f64;
    let mut out_of_range = Vec::new();

    log::info!("Column sums:");
    for (name, sum) in sums.iter().sorted_by_key(|(name, _)| *name) {
        let reference = if *sum > 2.0 { 100.0 } else { 1.0 };
        let tolerance = reference * tolerance_ratio;
        log::info!("{}\t{}", name, sum);
        if (sum - reference).abs() > tolerance {
            let path = sample_paths
                .get(name)
                .map(|path| path.display().to_string())
                .unwrap_or_else(|| "<unknown>".to_string());
            out_of_range.push((name.to_string(), *sum, reference, tolerance, path));
        }
    }

    if out_of_range.is_empty() {
        log::info!("Out-of-range column sums: none");
        log::info!("All column sums within reason.");
    } else {
        log::info!("Out-of-range column sums:");
        for (name, sum, reference, tolerance, path) in out_of_range
            .iter()
            .sorted_by_key(|(name, _, _, _, _)| name.as_str())
        {
            log::info!(
                "{}\t{}\t{}\t(|sum-{}| > {})",
                name,
                sum,
                path,
                reference,
                tolerance
            );
        }
        log::info!("Not all column sums are within reason.");
    }
}

fn entry_context<T: Taxonomy>(entry: &Entry<T>) -> String {
    let name = entry.taxon_name.as_deref().unwrap_or("<none>");
    let rank = entry.rank.to_string();
    let lineage = entry
        .lineage
        .as_ref()
        .map(|lineage| lineage.to_string())
        .unwrap_or_else(|| "<none>".to_string());
    format!(
        "taxon_name='{}', rank='{}', lineage='{}', abundance={}",
        name, rank, lineage, entry.abundance
    )
}

fn read_culprit_row(path: &Path, entry_index: usize) -> Option<String> {
    let file = File::open(path).ok()?;
    let reader = std::io::BufReader::new(file);
    let mut data_index = 0usize;
    for line in reader.lines().flatten() {
        let tokens: Vec<_> = line.split('\t').collect();
        if crate::format::Auto::skip_row(&tokens, Some(&["#", "@", "clade_name"]), None) {
            continue;
        }
        if data_index == entry_index {
            return Some(line);
        }
        data_index += 1;
    }
    None
}

#[cfg(test)]
mod tests {
    use super::{merge_profiles, MergeError};
    use crate::common::Taxonomy;
    use crate::profile::{Profile, ProfileWrapper};
    use crate::profile_handler::ProfileHandler;
    use std::collections::HashMap;
    use std::fs;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    use crate::common::TaxonomicRank;

    fn temp_path(name: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("Time went backwards")
            .as_nanos();
        std::env::temp_dir().join(format!("benchpro_{}_{}", name, nanos))
    }

    fn write_profile(path: &PathBuf, content: &str) {
        fs::write(path, content).expect("Failed to write test profile");
    }

    fn species_abundance_map_for<T: Taxonomy>(profile: &Profile<T>) -> HashMap<String, f64> {
        profile
            .get_taxa_string_dict(&TaxonomicRank::Species)
            .unwrap_or_default()
            .into_iter()
            .map(|(name, entries)| {
                let sum = entries.iter().map(|entry| entry.abundance).sum::<f64>();
                (name, sum)
            })
            .collect::<HashMap<_, _>>()
    }

    fn species_abundance_map(profile: &ProfileWrapper) -> HashMap<String, f64> {
        match profile {
            ProfileWrapper::GTDBProfile(inner) => species_abundance_map_for(inner),
            ProfileWrapper::NCBIProfile(inner) => species_abundance_map_for(inner),
            ProfileWrapper::ChocoPhlAnProfile(inner) => species_abundance_map_for(inner),
            ProfileWrapper::CustomProfile(inner) => species_abundance_map_for(inner),
        }
    }

    #[test]
    fn merges_profiles_into_matrix() -> Result<(), MergeError> {
        let profile_a = temp_path("profile_a.tsv");
        let profile_b = temp_path("profile_b.tsv");
        let output = temp_path("abundance.tsv");

        let header = "Lineage\tabundance\n";
        let content_a = format!(
            "{}{}\n{}",
            header,
            "d__Bacteria;p__Firmicutes;s__Foo\t0.3",
            "d__Bacteria;p__Firmicutes;s__Bar\t0.7"
        );
        let content_b = format!(
            "{}{}\n{}",
            header,
            "d__Bacteria;p__Firmicutes;s__Foo\t0.5",
            "d__Bacteria;p__Firmicutes;s__Baz\t0.5"
        );

        write_profile(&profile_a, &content_a);
        write_profile(&profile_b, &content_b);

        let _ = merge_profiles(
            &[profile_a.clone(), profile_b.clone()],
            &output,
            &TaxonomicRank::Species,
        )?;

        let output_content = fs::read_to_string(output).expect("Failed to read output");
        assert!(output_content.contains("Taxon"));
        assert!(output_content.contains("profile_a"));
        assert!(output_content.contains("profile_b"));
        assert!(output_content.contains("d__Bacteria;p__Firmicutes;s__Foo"));
        assert!(output_content.contains("d__Bacteria;p__Firmicutes;s__Bar"));
        assert!(output_content.contains("d__Bacteria;p__Firmicutes;s__Baz"));

        Ok(())
    }

    #[test]
    fn merges_genus_matrix_uses_full_lineage_and_aggregates_lower_ranks() -> Result<(), MergeError>
    {
        let profile = temp_path("profile_genus.tsv");
        let output = temp_path("abundance_genus.tsv");

        let header = "Lineage\tabundance\n";
        let content = format!(
            "{}{}\n{}\n{}",
            header,
            "d__Bacteria;p__Firmicutes;g__Alpha;s__Foo\t0.1",
            "d__Bacteria;p__Firmicutes;g__Alpha;s__Bar\t0.2",
            "d__Bacteria;p__Firmicutes;g__Beta;s__Baz\t0.7"
        );

        write_profile(&profile, &content);

        merge_profiles(&[profile.clone()], &output, &TaxonomicRank::Genus)?;

        let output_content = fs::read_to_string(&output).expect("Failed to read output");
        assert!(output_content.contains("d__Bacteria;p__Firmicutes;g__Alpha"));
        assert!(output_content.contains("d__Bacteria;p__Firmicutes;g__Beta"));

        let alpha_line = output_content
            .lines()
            .find(|line| line.starts_with("d__Bacteria;p__Firmicutes;g__Alpha\t"))
            .expect("Expected genus Alpha row in merged matrix");
        let tokens = alpha_line.split('\t').collect::<Vec<_>>();
        let alpha_value = tokens
            .get(1)
            .expect("Missing alpha abundance")
            .parse::<f64>()
            .expect("Invalid alpha abundance");

        assert!((alpha_value - 0.3).abs() < 1e-12);
        Ok(())
    }

    #[test]
    fn merges_duplicate_species_rows_by_summing_per_profile() -> Result<(), MergeError> {
        let profile_a = temp_path("profile_dup_a.tsv");
        let profile_b = temp_path("profile_dup_b.tsv");
        let output = temp_path("abundance_dup.tsv");

        let header = "Lineage\tabundance\n";
        let content_a = format!(
            "{}{}\n{}\n{}",
            header,
            "d__Bacteria;p__Firmicutes;s__Alpha\t0.1",
            "d__Bacteria;p__Firmicutes;s__Alpha\t0.2",
            "d__Bacteria;p__Firmicutes;s__Beta\t0.7"
        );
        let content_b = format!(
            "{}{}\n{}\n{}",
            header,
            "d__Bacteria;p__Firmicutes;s__Alpha\t0.3",
            "d__Bacteria;p__Firmicutes;s__Alpha\t0.4",
            "d__Bacteria;p__Firmicutes;s__Gamma\t0.3"
        );

        write_profile(&profile_a, &content_a);
        write_profile(&profile_b, &content_b);

        let (_sums, sample_paths) = merge_profiles(
            &[profile_a.clone(), profile_b.clone()],
            &output,
            &TaxonomicRank::Species,
        )?;

        // a) Ensure loading/merge recovered both input profiles.
        assert_eq!(
            sample_paths.len(),
            2,
            "Expected exactly two recovered sample files"
        );
        assert!(sample_paths.values().any(|path| path == &profile_a));
        assert!(sample_paths.values().any(|path| path == &profile_b));

        let output_content = fs::read_to_string(&output).expect("Failed to read output");
        let mut lines = output_content.lines();
        let header_tokens = lines
            .next()
            .expect("Missing matrix header")
            .split('\t')
            .collect::<Vec<_>>();
        assert_eq!(header_tokens.first().copied(), Some("Taxon"));

        let sample_a_name = profile_a
            .file_stem()
            .and_then(|name| name.to_str())
            .expect("Invalid profile_a file stem");
        let sample_b_name = profile_b
            .file_stem()
            .and_then(|name| name.to_str())
            .expect("Invalid profile_b file stem");

        let sample_a_col = header_tokens
            .iter()
            .position(|name| *name == sample_a_name)
            .expect("Missing profile_a matrix column");
        let sample_b_col = header_tokens
            .iter()
            .position(|name| *name == sample_b_name)
            .expect("Missing profile_b matrix column");

        let mut alpha_values: Option<(f64, f64)> = None;
        for line in lines {
            let tokens = line.split('\t').collect::<Vec<_>>();
            if tokens.first().copied() == Some("d__Bacteria;p__Firmicutes;s__Alpha") {
                let a = tokens
                    .get(sample_a_col)
                    .expect("Missing profile_a alpha abundance")
                    .parse::<f64>()
                    .expect("Invalid profile_a alpha abundance");
                let b = tokens
                    .get(sample_b_col)
                    .expect("Missing profile_b alpha abundance")
                    .parse::<f64>()
                    .expect("Invalid profile_b alpha abundance");
                alpha_values = Some((a, b));
                break;
            }
        }

        // b) Ensure species abundance is summed for duplicate rows per profile.
        let (alpha_a, alpha_b) = alpha_values.expect("Merged matrix missing s__Alpha row");
        assert!(
            (alpha_a - 0.3).abs() < 1e-12,
            "Expected profile_a s__Alpha sum 0.3, got {}",
            alpha_a
        );
        assert!(
            (alpha_b - 0.7).abs() < 1e-12,
            "Expected profile_b s__Alpha sum 0.7, got {}",
            alpha_b
        );

        Ok(())
    }

    #[test]
    fn matrix_roundtrip_matches_direct_profile_species_abundances() -> Result<(), MergeError> {
        let profile_a = temp_path("profile_roundtrip_a.tsv");
        let profile_b = temp_path("profile_roundtrip_b.tsv");
        let output = temp_path("abundance_roundtrip.tsv");

        let header = "Lineage\tabundance\n";
        let content_a = format!(
            "{}{}\n{}\n{}",
            header,
            "d__Bacteria;p__Firmicutes;s__Alpha\t0.15",
            "d__Bacteria;p__Firmicutes;s__Alpha\t0.35",
            "d__Bacteria;p__Firmicutes;s__Beta\t0.5"
        );
        let content_b = format!(
            "{}{}\n{}\n{}",
            header,
            "d__Bacteria;p__Firmicutes;s__Alpha\t0.4",
            "d__Bacteria;p__Firmicutes;s__Gamma\t0.4",
            "d__Bacteria;p__Firmicutes;s__Gamma\t0.2"
        );

        write_profile(&profile_a, &content_a);
        write_profile(&profile_b, &content_b);

        let (_sums, sample_paths) = merge_profiles(
            &[profile_a.clone(), profile_b.clone()],
            &output,
            &TaxonomicRank::Species,
        )?;
        assert_eq!(sample_paths.len(), 2);

        let sample_a_name = profile_a
            .file_stem()
            .and_then(|name| name.to_str())
            .expect("Invalid profile_a file stem");
        let sample_b_name = profile_b
            .file_stem()
            .and_then(|name| name.to_str())
            .expect("Invalid profile_b file stem");

        let matrix_a_ref = format!("matrix:{}::{}", output.display(), sample_a_name);
        let matrix_b_ref = format!("matrix:{}::{}", output.display(), sample_b_name);

        let direct_a =
            ProfileHandler::load_profile_auto(&profile_a).expect("Failed to load direct profile_a");
        let direct_b =
            ProfileHandler::load_profile_auto(&profile_b).expect("Failed to load direct profile_b");

        let matrix_a = ProfileHandler::load_profile(&matrix_a_ref, Some("GTDB"), None::<&str>)
            .expect("Failed to load matrix-derived profile_a");
        let matrix_b = ProfileHandler::load_profile(&matrix_b_ref, Some("GTDB"), None::<&str>)
            .expect("Failed to load matrix-derived profile_b");

        let direct_a_map = species_abundance_map(&direct_a);
        let direct_b_map = species_abundance_map(&direct_b);
        let matrix_a_map = species_abundance_map(&matrix_a);
        let matrix_b_map = species_abundance_map(&matrix_b);

        assert_eq!(direct_a_map.len(), matrix_a_map.len());
        assert_eq!(direct_b_map.len(), matrix_b_map.len());

        for (taxon, direct_abundance) in &direct_a_map {
            let matrix_abundance = matrix_a_map
                .get(taxon)
                .expect("Matrix profile_a is missing expected species");
            assert!(
                (direct_abundance - matrix_abundance).abs() < 1e-12,
                "Mismatch for profile_a taxon {}: direct={}, matrix={}",
                taxon,
                direct_abundance,
                matrix_abundance
            );
        }
        for (taxon, direct_abundance) in &direct_b_map {
            let matrix_abundance = matrix_b_map
                .get(taxon)
                .expect("Matrix profile_b is missing expected species");
            assert!(
                (direct_abundance - matrix_abundance).abs() < 1e-12,
                "Mismatch for profile_b taxon {}: direct={}, matrix={}",
                taxon,
                direct_abundance,
                matrix_abundance
            );
        }

        Ok(())
    }

    #[test]
    fn merge_promotes_unclassified_to_species_bucket() -> Result<(), MergeError> {
        let profile = temp_path("profile_unclassified.tsv");
        let output = temp_path("abundance_unclassified.tsv");
        let content = "\
#mpa_vJan25_CHOCOPhlAnSGB_202503
#clade_name\trelative_abundance
UNCLASSIFIED\t2.5
d__Bacteria;p__Firmicutes;s__Foo\t97.5
";
        write_profile(&profile, content);

        merge_profiles(&[profile], &output, &TaxonomicRank::Species)?;

        let output_content = fs::read_to_string(output).expect("Failed to read output");
        // The matrix keys every row by its full lineage -- `d__Bacteria;p__Firmicutes;s__Foo` --
        // and the unclassified mass follows that convention, so match the row by its species
        // suffix as the other tests here do rather than by a bare `s__` key at line start, which
        // the format never emits for any taxon.
        let unclassified: Vec<&str> = output_content
            .lines()
            .filter(|line| line.contains("s__UNCLASSIFIED"))
            .collect();
        // Exactly one row: the point of the test is that the unclassified mass is AGGREGATED into
        // the species bucket, so splitting it across two keys must fail rather than be found by
        // whichever one comes first.
        assert_eq!(
            unclassified.len(),
            1,
            "Expected exactly one s__UNCLASSIFIED row, got:\n{output_content}"
        );
        // Parsed, not substring-matched: `contains("\t2.5")` also accepts 2.55 and 2.5000001, so a
        // rounding or partial-double-count regression would pass it.
        let abundance: f64 = unclassified[0]
            .rsplit('\t')
            .next()
            .and_then(|value| value.parse().ok())
            .unwrap_or_else(|| panic!("no numeric abundance in:\n{}", unclassified[0]));
        assert!(
            (abundance - 2.5).abs() < 1e-9,
            "Expected the unclassified mass to be 2.5, got {abundance} in:\n{}",
            unclassified[0]
        );

        Ok(())
    }

    // Regression test for sylph/CAMI-style profiles that list every rank
    // explicitly. These profiles have a strain entry for every species entry
    // with the same abundance value. Without the `has_explicit_target_rank`
    // guard, both the species entry and its child strain entry would be mapped
    // to the same species key, doubling all abundances (~200 instead of ~100).
    #[test]
    fn no_double_count_when_profile_has_species_and_strain_entries() -> Result<(), MergeError> {
        let profile = temp_path("profile_sylph_like.tsv");
        let output = temp_path("abundance_sylph_like.tsv");

        // Each species has a corresponding strain with the same abundance,
        // mirroring what sylph outputs.
        let content = "\
Lineage\tabundance
d__Bacteria;p__Firmicutes;g__Foo;s__Foo alpha\t60.0
d__Bacteria;p__Firmicutes;g__Bar;s__Bar beta\t40.0
d__Bacteria;p__Firmicutes;g__Foo;s__Foo alpha;t__GCF_001\t60.0
d__Bacteria;p__Firmicutes;g__Bar;s__Bar beta;t__GCF_002\t40.0
";
        write_profile(&profile, content);

        merge_profiles(&[profile], &output, &TaxonomicRank::Species)?;

        let output_content = fs::read_to_string(output).expect("Failed to read output");

        // Each species row should appear exactly once, not doubled.
        let foo_line = output_content
            .lines()
            .find(|l| l.contains("s__Foo alpha"))
            .expect("Missing s__Foo alpha row");
        let foo_val: f64 = foo_line
            .split('\t')
            .nth(1)
            .expect("Missing value column")
            .parse()
            .expect("Non-numeric value");
        assert!(
            (foo_val - 60.0).abs() < 1e-9,
            "Expected s__Foo alpha = 60, got {} (double-counting?)",
            foo_val
        );

        let bar_line = output_content
            .lines()
            .find(|l| l.contains("s__Bar beta"))
            .expect("Missing s__Bar beta row");
        let bar_val: f64 = bar_line
            .split('\t')
            .nth(1)
            .expect("Missing value column")
            .parse()
            .expect("Non-numeric value");
        assert!(
            (bar_val - 40.0).abs() < 1e-9,
            "Expected s__Bar beta = 40, got {} (double-counting?)",
            bar_val
        );

        Ok(())
    }
}
