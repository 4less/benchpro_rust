use std::{
    cmp::max,
    collections::{HashMap, HashSet},
    default,
    fmt::Display,
    io::{BufRead, BufReader, Read, Seek},
    process::exit,
    str::FromStr,
};

use itertools::Itertools;
use log::{debug, error, warn};
use polars::{error::PolarsResult, frame::DataFrame, prelude::NamedFrom, series::Series};

use crate::{
    common::{
        Custom, Lineage, LineageFromString, Taxon, TaxonomicRank, Taxonomy, TaxonomyEnum, GTDB,
        NCBI,
    },
    format::{Columns, Format},
    utils::add_string_columns,
};

pub type StringList = Vec<String>;
pub type StringStringMap = HashMap<String, String>;
pub type Entries<T: Taxonomy> = Vec<Entry<T>>;
pub type EntriesRef<'a, T: Taxonomy> = Vec<&'a Entry<T>>;

/// Single taxon abundance entry parsed from a profile.
///
/// # Fields
///
/// - `taxon_name` - Optional raw taxon name as provided in the input
/// - `lineage` - Optional lineage parsed into taxonomy-specific ranks
/// - `alternative_names` - Optional alternative names for matching
/// - `abundance` - Relative abundance for this taxon
/// - `rank` - Parsed taxonomic rank for this entry
#[derive(Debug, Default, Clone, PartialEq)]
pub struct Entry<T: Taxonomy> {
    pub taxon_name: Option<String>,
    pub lineage: Option<Lineage<T>>,
    pub alternative_names: Option<Vec<String>>,
    pub abundance: f64,
    pub rank: TaxonomicRank,
}

impl<T: Taxonomy> Entry<T> {
    /// Returns the lineage for this entry if available.
    ///
    /// # Returns
    ///
    /// Optional reference to the lineage.
    pub fn lineage(&self) -> Option<&Lineage<T>> {
        self.lineage.as_ref()
    }
}

#[derive(Default, Clone, Debug)]
struct BCVectors {
    pred_type: Vec<String>,
    taxon_names: Vec<String>,
    taxon_ids: Vec<usize>,
    prediction_abundances: Vec<f64>,
    gold_std_abundances: Vec<f64>,
    prediction_count: Vec<usize>,
    gold_std_count: Vec<usize>,
}

/// Binary classification label used in benchmarking.
#[derive(Debug, Clone, PartialEq)]
pub enum BC {
    TP,
    FP,
    FN,
    TN,
    FFP,
    FFN,
    Unknown,
}

impl Default for BC {
    fn default() -> Self {
        BC::Unknown
    }
}

impl BC {
    /// Returns the short string representation used in reports.
    ///
    /// # Returns
    ///
    /// Short label such as `TP`, `FP`, or `FN`.
    pub fn to_string(&self) -> String {
        match self {
            BC::TP => "TP".to_string(),
            BC::FP => "FP".to_string(),
            BC::FN => "FN".to_string(),
            BC::TN => "TN".to_string(),
            BC::FFP => "FFP".to_string(),
            BC::FFN => "FFN".to_string(),
            BC::Unknown => "Unknown".to_string(),
        }
    }
}

impl TryFrom<&str> for BC {
    type Error = String;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            value if value == "TP" => Ok(BC::TP),
            value if value == "FP" => Ok(BC::FP),
            value if value == "FN" => Ok(BC::FN),
            value if value == "TN" => Ok(BC::TN),
            value if value == "FFP" => Ok(BC::FFP),
            value if value == "FFN" => Ok(BC::FFN),
            value if value == "Unknown" => Ok(BC::Unknown),
            _ => Err(format!("BC Type does not exist: {}", value)),
        }
    }
}

impl FromStr for BC {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        Self::try_from(value)
    }
}

impl BCVectors {
    pub fn add_entry(
        &mut self,
        pred_type: &BC,
        taxon_name: &str,
        taxon_id: usize,
        prediction_abundance: f64,
        gold_std_abundance: f64,
        prediction_count: usize,
        gold_std_count: usize,
    ) {
        self.pred_type.push(pred_type.to_string());
        self.taxon_names.push(taxon_name.to_string());
        self.taxon_ids.push(taxon_id);
        self.prediction_abundances.push(prediction_abundance);
        self.gold_std_abundances.push(gold_std_abundance);
        self.prediction_count.push(prediction_count);
        self.gold_std_count.push(gold_std_count);
    }

    pub fn into_polars_df(self) -> PolarsResult<DataFrame> {
        let df = DataFrame::new(vec![
            Series::new("Name".into(), self.taxon_names),
            Series::new("Type".into(), self.pred_type),
            Series::new("PredictionAbundance".into(), self.prediction_abundances),
            Series::new("GoldStdAbundance".into(), self.gold_std_abundances),
            Series::new(
                "PredictionCount".into(),
                self.prediction_count
                    .iter()
                    .map(|&x| x as i64)
                    .collect::<Vec<_>>(),
            ),
            Series::new(
                "GoldStdCount".into(),
                self.gold_std_count
                    .iter()
                    .map(|&x| x as i64)
                    .collect::<Vec<_>>(),
            ),
        ]);
        df
    }

    pub fn polars_df(&self) -> PolarsResult<DataFrame> {
        let df = DataFrame::new(vec![
            Series::new("Name".into(), self.taxon_names.clone()),
            Series::new("Type".into(), self.pred_type.clone()),
            Series::new(
                "PredictionAbundance".into(),
                self.prediction_abundances.clone(),
            ),
            Series::new("GoldStdAbundance".into(), self.gold_std_abundances.clone()),
            Series::new(
                "PredictionCount".into(),
                self.prediction_count
                    .clone()
                    .iter()
                    .map(|&x| x as i64)
                    .collect::<Vec<_>>(),
            ),
            Series::new(
                "GoldStdCount".into(),
                self.gold_std_count
                    .clone()
                    .iter()
                    .map(|&x| x as i64)
                    .collect::<Vec<_>>(),
            ),
        ]);
        df
    }
}

type TaxonEntryMap<'a, T: Taxonomy> = HashMap<&'a Taxon, EntriesRef<'a, T>>;
type NamesEntryMap<'a, T: Taxonomy> = HashMap<String, EntriesRef<'a, T>>;

fn is_prefix_only_taxon_name(name: &str) -> bool {
    let trimmed = name.trim();
    let Some((_, rest)) = trimmed.split_once("__") else {
        return false;
    };

    rest.trim().is_empty()
}

fn filter_prefix_only_taxa_names<'a, T: Taxonomy>(
    taxa: NamesEntryMap<'a, T>,
) -> NamesEntryMap<'a, T> {
    taxa.into_iter()
        .filter(|(name, _)| !is_prefix_only_taxon_name(name))
        .collect()
}

fn filter_prefix_only_taxa_refs<'a, T: Taxonomy>(
    taxa: TaxonEntryMap<'a, T>,
) -> TaxonEntryMap<'a, T> {
    taxa.into_iter()
        .filter(|(taxon, _)| !is_prefix_only_taxon_name(&taxon.name))
        .collect()
}

/// Computes TP/FP/FN counts based on exact taxon name matches.
///
/// # Arguments
///
/// * `prediction` - Predicted taxa grouped by name
/// * `gold_std` - Gold standard taxa grouped by name
///
/// # Returns
///
/// DataFrame containing TP/FP/FN rows with abundances and counts.
///
/// # Errors
///
/// Returns a Polars error when DataFrame creation fails.
pub fn binary_classification_df<T: Taxonomy>(
    prediction: &NamesEntryMap<T>,
    gold_std: &NamesEntryMap<T>,
) -> PolarsResult<DataFrame> {
    let mut df_vectors = BCVectors::default();

    let pred_set = prediction.keys().collect::<HashSet<_>>();
    let gold_set = gold_std.keys().collect::<HashSet<_>>();

    let tps = pred_set
        .intersection(&gold_set)
        .map(|&x| x)
        .collect::<HashSet<_>>();
    let fps = pred_set
        .difference(&gold_set)
        .map(|&x| x)
        .collect::<HashSet<_>>();
    let fns = gold_set
        .difference(&pred_set)
        .map(|&x| x)
        .collect::<HashSet<_>>();

    tps.iter()
        .map(|&t| (BC::TP, t))
        .chain(fps.iter().map(|&t| (BC::FP, t)))
        .chain(fns.iter().map(|&t| (BC::FN, t)))
        .for_each(|(mut bp_type, taxon)| {
            let (prediction_abundance, prediction_entry_count) =
                if matches!(bp_type, BC::TP) || matches!(bp_type, BC::FP) {
                    let entries = prediction.get(taxon).expect(&format!(
                        "Gold set does not contain TP taxon. Unrecoverable error. '{}', {}",
                        taxon,
                        bp_type.to_string()
                    ));
                    (entries.iter().map(|e| e.abundance).sum(), entries.len())
                } else {
                    (0f64, 0usize)
                };

            let (gold_std_abundance, gold_std_entry_count) =
                if matches!(bp_type, BC::TP) || matches!(bp_type, BC::FN) {
                    let entries = gold_std.get(taxon).expect(&format!(
                        "Gold set does not contain FN taxon. Unrecoverable error. '{}', {}",
                        taxon,
                        bp_type.to_string()
                    ));
                    (entries.iter().map(|e| e.abundance).sum(), entries.len())
                } else {
                    (0f64, 0usize)
                };

            if taxon == "" && matches!(bp_type, BC::FP) {
                debug!("FIX BC TYPE");
                bp_type = BC::Unknown;
            }

            df_vectors.add_entry(
                &bp_type,
                taxon,
                0,
                prediction_abundance,
                gold_std_abundance,
                prediction_entry_count,
                gold_std_entry_count,
            );
        });
    df_vectors.into_polars_df()
}

/// Computes binary classification counts when alternative names are allowed.
///
/// # Arguments
///
/// * `prediction` - Predicted taxa mapped by taxon object
/// * `gold_std` - Gold standard taxa mapped by taxon object
///
/// # Returns
///
/// DataFrame containing TP/FP/FN entries with abundance summaries.
///
/// # Errors
///
/// Returns a Polars error when the output DataFrame cannot be created.
pub fn binary_classification_alternatives_df<T: Taxonomy>(
    prediction: &TaxonEntryMap<T>,
    gold_std: &TaxonEntryMap<T>,
) -> PolarsResult<DataFrame> {
    let mut df_vectors = BCVectors::default();

    let prediction_taxa = prediction.keys().collect::<HashSet<_>>();
    let gold_std_taxa = gold_std.keys().collect::<HashSet<_>>();

    let BCData {
        tps,
        fps,
        fns,
        prediction_counts,
    } = binary_classification_helper(&prediction_taxa, &gold_std_taxa);

    if tps
        .iter()
        .any(|(_, v)| v.secondary.is_some() && v.primary.is_none())
    {
        for (taxon, value) in tps.iter() {
            debug!("Gold {:?}\n\t\tPred {:?}", taxon, value);
        }
    }

    tps.iter().for_each(|(&taxon, predictions)| {
        let (prediction_abundance, prediction_entry_count) = {
            // if primary is set, primary abundance always counts.
            // if secondary is set, only secondary entries count that only match ONE gold std entry
            // Ambiguous matches with more than one gold_std are not counted.

            let primary_abundance_list = predictions.primary.as_ref().map_or(vec![], |x| {
                prediction
                    .get(taxon)
                    .unwrap()
                    .iter()
                    .map(|t| t.abundance)
                    .collect_vec()
            });
            let primary_abundance: f64 = primary_abundance_list.iter().sum();
            let primary_entry_count = primary_abundance_list.len();

            let unique_secondary_abundance_list =
                predictions.secondary.as_ref().map_or(vec![], |x| {
                    x.iter()
                        .filter(|&&x| prediction_counts.get(x).expect("asdads") == &1)
                        .flat_map(|t| prediction.get(t).unwrap())
                        .map(|&e| e.abundance)
                        .collect_vec()
                });

            let non_unique_secondary_abundance_list =
                predictions.secondary.as_ref().map_or(vec![], |x| {
                    x.iter()
                        .filter(|&&x| prediction_counts.get(x).expect("asdads") > &1)
                        .flat_map(|t| prediction.get(t).unwrap())
                        .map(|&e| e.abundance)
                        .collect_vec()
                });
            let unique_secondary_abundance: f64 = unique_secondary_abundance_list.iter().sum();
            let non_unique_secondary_abundance: f64 =
                non_unique_secondary_abundance_list.iter().sum();
            let unique_secondary_count = unique_secondary_abundance_list.len();
            let non_unique_secondary_count = non_unique_secondary_abundance_list.len();

            if non_unique_secondary_abundance != 0.0f64 {
                debug!("---- {:?}", taxon);
                debug!(
                    "Abundance: Primary {}  Secondary (U) {} Secondary (NU) {}",
                    primary_abundance, unique_secondary_abundance, non_unique_secondary_abundance
                );
                debug!(
                    "Counts:    Primary {}  Secondary (U) {} Secondary (NU) {}",
                    primary_entry_count, unique_secondary_count, non_unique_secondary_count
                );

                for g in gold_std {
                    debug!("GOLD: {:?}", g.0);
                }
                for p in prediction {
                    debug!("PRED: {:?}", p.0);
                }
            }

            (
                primary_abundance + unique_secondary_abundance,
                primary_entry_count + unique_secondary_count,
            )
        };

        let (gold_std_abundance, gold_std_entry_count) = {
            let entries = match gold_std.get(taxon) {
                Some(entries) => entries,
                None => {
                    debug!("--------------------------------");
                    gold_std.iter().for_each(|(t, _)| {
                        debug!("Gold: {} {}", t.match_exact(taxon), t.match_any(taxon));
                    });
                    panic!(
                        "Gold set does not contain FN taxon. Unrecoverable error. '{:?}', TP",
                        taxon
                    )
                }
            };
            (entries.iter().map(|e| e.abundance).sum(), entries.len())
        };

        df_vectors.add_entry(
            &BC::TP,
            &taxon.name,
            0,
            prediction_abundance,
            gold_std_abundance,
            prediction_entry_count,
            gold_std_entry_count,
        );
    });

    fps.iter().map(|&t| (BC::FP, t))
        .chain(fns.iter().map(|&t| (BC::FN, t)))
        .for_each(|(mut bp_type, taxon)|
    {
        debug!("{}: {:?}", bp_type.to_string(), taxon);

        let (prediction_abundance, prediction_entry_count) = if matches!(bp_type, BC::TP) || matches!(bp_type, BC::FP) {

            let entries = match prediction.get(taxon) {
                Some(entries) => entries,
                None => {
                    debug!("--------------------------------");
                    gold_std.iter().for_each(|(t, _)| {
                        debug!("Gold: {} {}", t.match_exact(taxon), t.match_any(taxon));
                    });
                    panic!("Prediction set does not contain P taxon. Unrecoverable error. '{:?}', {}", taxon, bp_type.to_string())
                },
            };
            let entries = prediction.get(taxon)
            .expect(&format!("Gold set does not contain TP taxon. Unrecoverable error. '{:?}', {}", taxon, bp_type.to_string()));
            (
                 entries
                     .iter().map(|e| e.abundance)
                     .sum(),
                 entries.len()
            )
        } else { (0f64, 0usize) };



        let (gold_std_abundance, gold_std_entry_count) = if matches!(bp_type, BC::TP) || matches!(bp_type, BC::FN) {
            let entries = match gold_std.get(taxon) {
                Some(entries) => entries,
                None => {
                    debug!("--------------------------------");
                    gold_std.iter().for_each(|(t, _)| {
                        debug!("Gold: {} {}", t.match_exact(taxon), t.match_any(taxon));
                    });
                    panic!("Gold set does not contain FN taxon. Unrecoverable error. '{:?}', {}", taxon, bp_type.to_string())
                },
            };
            (
                entries
                    .iter().map(|e| e.abundance)
                    .sum(),
                entries.len()
            )
        } else { (0f64, 0usize) };

        // Later check
        if taxon.name == "" && matches!(bp_type, BC::FP) {
            debug!("FIX BC TYPE");
            bp_type = BC::Unknown;
        }

        df_vectors.add_entry(
            &bp_type,
            &taxon.name,
            0,
            prediction_abundance,
            gold_std_abundance,
            prediction_entry_count,
            gold_std_entry_count);
    });
    df_vectors.into_polars_df()
}

#[derive(Default, Debug)]
struct TPMapValue<'a> {
    primary: Option<&'a Taxon>,
    secondary: Option<Vec<&'a Taxon>>,
}

pub type TPMap<'a> = HashMap<&'a Taxon, TPMapValue<'a>>;
pub type OtherSet<'a> = HashSet<&'a Taxon>;
pub type MatchCountMap<'a> = HashMap<&'a Taxon, usize>;
/// Aggregated TP/FP/FN sets and match counts for alternative matching.
#[derive(Default)]
pub struct BCData<'a> {
    tps: TPMap<'a>,
    fps: OtherSet<'a>,
    fns: OtherSet<'a>,
    prediction_counts: MatchCountMap<'a>,
}

/// Computes TP/FP/FN sets and match counts for alternative-name matching.
///
/// # Arguments
///
/// * `prediction` - Predicted taxa set
/// * `gold_std` - Gold standard taxa set
///
/// # Returns
///
/// Aggregated TP/FP/FN sets with prediction match counts.
pub fn binary_classification_helper<'a, TREF: AsRef<Taxon>>(
    prediction: &'a HashSet<TREF>,
    gold_std: &'a HashSet<TREF>,
) -> BCData<'a> {
    let mut fn_return = BCData::default();

    let mut matched_gold_stds = HashSet::new();
    let mut matched_predictions = HashSet::new();

    let all_gold_no_alt = gold_std
        .iter()
        .all(|x| x.as_ref().alternative_names.is_none());

    if !all_gold_no_alt {
        for gs in gold_std.iter() {
            debug!("{:?}", gs.as_ref());
        }

        error!("Gold standard may not have any alternative names.");
        panic!("Gold standard may not have any alternative names.")
    }

    // iterator over all pairs of gold std and prediction
    let pair_iter = gold_std
        .iter()
        .flat_map(|gold| prediction.iter().map(|pred| (gold.as_ref(), pred.as_ref())));

    // iterator over exact match pairs.
    let exact_match_pairs = pair_iter.clone().filter(|(gold, prediction)| {
        assert!(gold.alternative_names.is_none());
        gold.match_exact(prediction)
    });

    // iterator over alternative name match pairs.
    let alternative_match_pairs = pair_iter.filter(|(gold, prediction)| {
        assert!(gold.alternative_names.is_none());
        gold.match_any(prediction) && !gold.match_exact(prediction)
    });

    // Extract exact matches first
    exact_match_pairs.clone().for_each(|(gold, prediction)| {
        assert!(gold.alternative_names.is_none());
        if matched_gold_stds.contains(gold) {
            warn!("Entry: {:?}", matched_gold_stds.get(gold).unwrap());
            warn!("GoldEntry: {:?}", fn_return.tps.get(gold).unwrap());
            warn!("Pred: {:?}", prediction);

            exact_match_pairs.clone().for_each(|(g, p)| {
                debug!("----MATCH\nG: {:?}\nP: {:?}", g, p);
            });
        }
        assert!(!matched_gold_stds.contains(gold));
        assert!(!fn_return.tps.contains_key(gold));
        matched_gold_stds.insert(gold);

        matched_predictions.insert(prediction);
        fn_return.tps.insert(
            gold,
            TPMapValue {
                primary: Some(gold),
                secondary: None,
            },
        );
    });

    alternative_match_pairs.for_each(|(gold, prediction)| {
        let entry = fn_return.tps.entry(gold).or_insert(TPMapValue::default());

        matched_predictions.insert(prediction);
        if let Some(list) = &mut entry.secondary {
            list.push(prediction);
        } else {
            entry.secondary = Some(vec![prediction]);
        }

        fn_return
            .prediction_counts
            .entry(prediction)
            .and_modify(|v| *v += 1)
            .or_insert(1);
    });

    fn_return.fps = prediction
        .iter()
        .filter(|p| !matched_predictions.contains(p.as_ref()))
        .map(|x| x.as_ref())
        .collect::<HashSet<_>>();

    fn_return.fns = gold_std
        .iter()
        .filter(|g| !matched_gold_stds.contains(g.as_ref()))
        .map(|x| x.as_ref())
        .collect::<HashSet<_>>();

    fn_return
}

/// Parsed taxonomic profile plus optional metadata.
#[derive(Debug, Default, Clone, PartialEq)]
pub struct Profile<T: Taxonomy> {
    pub unstructured_meta: StringList,
    pub meta: StringStringMap,
    pub taxa: Entries<T>,
}

/// Report of abundance normalization decisions and per-rank sums.
#[derive(Debug, Clone, PartialEq)]
pub struct NormalizeReport {
    /// Whether normalization was applied.
    pub applied: bool,
    /// Per-rank sums before normalization.
    pub pre_sums: HashMap<TaxonomicRank, f64>,
    /// Per-rank sums after normalization.
    pub post_sums: HashMap<TaxonomicRank, f64>,
}

/// Allowed relative tolerance when validating abundance sums.
pub const ABUNDANCE_SUM_TOLERANCE: f64 = 0.20;

impl Profile<GTDB> {
    /// Wraps a GTDB profile for type-erased handling.
    ///
    /// # Returns
    ///
    /// `ProfileWrapper::GTDBProfile`.
    pub fn wrap(self) -> ProfileWrapper {
        ProfileWrapper::GTDBProfile(self)
    }
}

impl Profile<NCBI> {
    /// Wraps an NCBI profile for type-erased handling.
    ///
    /// # Returns
    ///
    /// `ProfileWrapper::NCBIProfile`.
    pub fn wrap(self) -> ProfileWrapper {
        ProfileWrapper::NCBIProfile(self)
    }
}

impl Profile<Custom> {
    /// Wraps a custom profile for type-erased handling.
    ///
    /// # Returns
    ///
    /// `ProfileWrapper::CustomProfile`.
    pub fn wrap(self) -> ProfileWrapper {
        ProfileWrapper::CustomProfile(self)
    }
}

impl<T: Taxonomy> Profile<T> {
    /// Returns the set of unique ranks present in the profile.
    ///
    /// # Returns
    ///
    /// `None` when ranks are unknown; otherwise the set of ranks.
    pub fn unique_ranks(&self) -> Option<HashSet<TaxonomicRank>> {
        if self
            .taxa
            .iter()
            .all(|entry| matches!(entry.rank, TaxonomicRank::Unknown))
        {
            return None;
        }
        assert!(self
            .taxa
            .iter()
            .all(|entry| !matches!(entry.rank, TaxonomicRank::Unknown)));

        Some(
            self.taxa
                .iter()
                .map(|e| e.rank.clone())
                .collect::<HashSet<_>>(),
        )
    }

    fn rank_sums(&self) -> HashMap<TaxonomicRank, f64> {
        let mut rank_sums: HashMap<TaxonomicRank, f64> = HashMap::new();
        for entry in &self.taxa {
            let sum = rank_sums.entry(entry.rank.clone()).or_insert(0.0);
            *sum += entry.abundance;
        }
        rank_sums
    }

    /// Returns rank abundance sums that exceed a maximum threshold.
    ///
    /// # Arguments
    ///
    /// * `max_sum` - Maximum allowed total abundance per rank
    ///
    /// # Returns
    ///
    /// List of `(rank, sum)` pairs that exceed `max_sum`.
    pub fn rank_sum_violations(&self, max_sum: f64) -> Vec<(TaxonomicRank, f64)> {
        self.rank_sums()
            .into_iter()
            .filter(|(_, sum)| *sum > max_sum)
            .collect()
    }

    /// Normalizes abundances to fractions when sums indicate percent-scale data.
    ///
    /// Normalization is decided per sample, not per rank. When multiple ranks are
    /// explicitly present, either all ranks must be within tolerance of 1.0 or
    /// all ranks must be within tolerance of 100.0, otherwise an error is raised.
    ///
    /// # Arguments
    ///
    /// * `tolerance` - Relative tolerance for deciding whether sums are near 1.0 or 100.0
    ///
    /// # Returns
    ///
    /// Report containing pre/post rank sums and whether normalization was applied.
    ///
    /// # Errors
    ///
    /// Returns `ProfileError::NormalizationError` when sums are inconsistent or
    /// remain out of tolerance after normalization.
    pub fn normalize_abundances_with_checks(
        &mut self,
        tolerance: f64,
    ) -> Result<NormalizeReport, ProfileError> {
        let pre_sums = self.rank_sums();
        if pre_sums.is_empty() {
            return Ok(NormalizeReport {
                applied: false,
                pre_sums,
                post_sums: HashMap::new(),
            });
        }

        let within = |value: f64, target: f64| (value - target).abs() <= target * tolerance;
        let active_threshold = tolerance * 1.0;
        let active_pre_sums = pre_sums
            .iter()
            .filter(|(_, sum)| **sum >= active_threshold)
            .map(|(rank, sum)| (rank.clone(), *sum))
            .collect::<HashMap<_, _>>();

        if active_pre_sums.is_empty() {
            return Err(ProfileError::NormalizationError(
                "No active rank sums above threshold; cannot normalize".to_owned(),
            ));
        }

        let all_within_one = active_pre_sums.values().all(|sum| within(*sum, 1.0));
        let all_within_hundred = active_pre_sums.values().all(|sum| within(*sum, 100.0));

        if all_within_one {
            return Ok(NormalizeReport {
                applied: false,
                pre_sums: pre_sums.clone(),
                post_sums: pre_sums,
            });
        }

        if all_within_hundred {
            for entry in &mut self.taxa {
                entry.abundance /= 100.0;
            }

            let post_sums = self.rank_sums();
            let active_post_sums = post_sums
                .iter()
                .filter(|(_, sum)| **sum >= active_threshold)
                .map(|(rank, sum)| (rank.clone(), *sum))
                .collect::<HashMap<_, _>>();
            let all_post_within_one = active_post_sums.values().all(|sum| within(*sum, 1.0));
            if !all_post_within_one {
                return Err(ProfileError::NormalizationError(format!(
                    "Abundance normalization failed; rank sums not within tolerance of 1.0 after normalization: {:?}",
                    post_sums
                )));
            }

            return Ok(NormalizeReport {
                applied: true,
                pre_sums,
                post_sums,
            });
        }

        Err(ProfileError::NormalizationError(format!(
            "Abundance sums are not within tolerance of 1.0 or 100.0; pre-normalization sums: {:?}",
            pre_sums
        )))
    }

    /// If ranks are defined and there is more than one rank,
    /// I deduce that ranks are defined separately on purpose.
    /// If ranks are undefined, see if the lineage contains the target ranks.
    /// - This only works for the GTDB taxonomy.
    /// - For the NCBI taxonomy, something more elaborate needss to be figured out
    ///
    /// # Arguments
    ///
    /// * `rank` - Target rank to extract
    ///
    /// # Returns
    ///
    /// `Some` set of taxon names for the rank or `None` if unavailable.
    pub fn get_taxa_with_rank(&self, rank: &TaxonomicRank) -> Option<HashSet<String>> {
        match self.unique_ranks() {
            Some(ranks) => {
                if (ranks.len() > 1 && !ranks.contains(rank))
                    || rank > ranks.iter().max().expect("Has no max")
                {
                    return None;
                }
                if ranks.len() == 1 && rank <= ranks.iter().max().expect("Has no max") {
                    let set = self
                        .taxa
                        .iter()
                        .map(|entry| entry.lineage().unwrap().get(rank))
                        .filter(|name| name.is_some())
                        .map(|taxon| taxon.unwrap())
                        .map(|taxon| taxon.name.clone())
                        .collect::<HashSet<_>>();
                    if set.is_empty() {
                        return None;
                    } else {
                        return Some(set);
                    }
                }

                let set = self
                    .taxa
                    .iter()
                    .filter(|entry| matches!(&entry.rank, rank))
                    .map(|entry| entry.lineage().unwrap().get(rank))
                    .filter(|name| name.is_some())
                    .map(|taxon| taxon.unwrap())
                    .map(|taxon| taxon.name.clone())
                    .collect::<HashSet<_>>();

                if set.is_empty() {
                    return None;
                } else {
                    return Some(set);
                }
            }
            None => {
                let set = self
                    .taxa
                    .iter()
                    .map(|entry| entry.lineage().unwrap().get(rank))
                    .filter(|name| name.is_some())
                    .map(|taxon| taxon.unwrap())
                    .map(|taxon| taxon.name.clone())
                    .collect::<HashSet<_>>();

                if set.is_empty() {
                    return None;
                } else {
                    return Some(set);
                }
            }
        }
    }

    /// Groups entries by taxon name at the requested rank.
    ///
    /// # Arguments
    ///
    /// * `rank` - Target rank to extract
    ///
    /// # Returns
    ///
    /// Map from taxon name to matching entries, or `None` if the rank is absent.
    pub fn get_taxa_string_dict(
        &self,
        rank: &TaxonomicRank,
    ) -> Option<HashMap<String, Vec<&Entry<T>>>> {
        match self.unique_ranks() {
            Some(ranks) => {
                if (ranks.len() > 1 && !ranks.contains(rank))
                    || rank > ranks.iter().max().expect("Has no max")
                {
                    return None;
                }
                if ranks.len() == 1 && rank <= ranks.iter().max().expect("Has no max") {
                    debug!("Option 1: There is only one rank defined {:?}", ranks);
                    let set = self
                        .taxa
                        .iter()
                        .map(|entry| (entry.lineage().unwrap().get(rank), entry))
                        .filter(|(name, entry)| name.is_some())
                        .map(|(taxon, entry)| (taxon.unwrap(), entry))
                        .map(|(taxon, entry)| (taxon.name.clone(), entry))
                        .fold(HashMap::new(), |mut acc, (key, value)| {
                            acc.entry(key).or_insert_with(Vec::new).push(value);
                            acc
                        });

                    if set.is_empty() {
                        return None;
                    } else {
                        return Some(set);
                    }
                }
                debug!(
                    "Option 2: There are multiple ranks defined {:?} (Target: {:?})",
                    ranks, rank
                );
                let set = self
                    .taxa
                    .iter()
                    .filter(|&entry| rank == &entry.rank)
                    .map(|entry| (entry.lineage().unwrap().get(rank), entry))
                    .filter(|(name, _)| name.is_some())
                    .map(|(taxon, entry)| (taxon.unwrap(), entry))
                    .map(|(taxon, entry)| (taxon.name.clone(), entry))
                    .fold(HashMap::new(), |mut acc, (key, value)| {
                        acc.entry(key).or_insert_with(Vec::new).push(value);
                        acc
                    });

                if set.is_empty() {
                    return None;
                } else {
                    return Some(set);
                }
            }
            None => {
                debug!("Option3");
                let set = self
                    .taxa
                    .iter()
                    .map(|entry| (entry.lineage().unwrap().get(rank), entry))
                    .filter(|(name, entry)| name.is_some())
                    .map(|(taxon, entry)| (taxon.unwrap(), entry))
                    .map(|(taxon, entry)| (taxon.name.clone(), entry))
                    .fold(HashMap::new(), |mut acc, (key, value)| {
                        acc.entry(key).or_insert_with(Vec::new).push(value);
                        acc
                    });

                if set.is_empty() {
                    return None;
                } else {
                    return Some(set);
                }
            }
        }
    }

    /// Groups entries by taxon object at the requested rank.
    ///
    /// # Arguments
    ///
    /// * `rank` - Target rank to extract
    ///
    /// # Returns
    ///
    /// Map from taxon reference to matching entries, or `None` if the rank is absent.
    pub fn get_taxa_dict<'a>(
        &'a self,
        rank: &TaxonomicRank,
    ) -> Option<HashMap<&'a Taxon, Vec<&Entry<T>>>> {
        match self.unique_ranks() {
            Some(ranks) => {
                if (ranks.len() > 1 && !ranks.contains(rank))
                    || rank > ranks.iter().max().expect("Has no max")
                {
                    return None;
                }
                if ranks.len() == 1 && rank <= ranks.iter().max().expect("Has no max") {
                    let set = self
                        .taxa
                        .iter()
                        .map(|entry| (entry.lineage().unwrap().get(rank), entry))
                        .filter(|(name, _)| name.is_some())
                        .map(|(taxon, entry)| (taxon.unwrap(), entry))
                        .fold(HashMap::new(), |mut acc, (key, value)| {
                            acc.entry(key).or_insert_with(Vec::new).push(value);
                            acc
                        });

                    if set.is_empty() {
                        return None;
                    } else {
                        return Some(set);
                    }
                }
                let set = self
                    .taxa
                    .iter()
                    .filter(|&entry| &entry.rank == rank)
                    .map(|entry| (entry.lineage().unwrap().get(rank), entry))
                    .filter(|(name, _)| name.is_some())
                    .map(|(taxon, entry)| (taxon.unwrap(), entry))
                    .fold(HashMap::new(), |mut acc, (key, value)| {
                        acc.entry(key).or_insert_with(Vec::new).push(value);
                        acc
                    });

                if set.is_empty() {
                    return None;
                } else {
                    return Some(set);
                }
            }
            None => {
                let set = self
                    .taxa
                    .iter()
                    .map(|entry| (entry.lineage().unwrap().get(rank), entry))
                    .filter(|(name, _)| name.is_some())
                    .map(|(taxon, entry)| (taxon.unwrap(), entry))
                    .fold(HashMap::new(), |mut acc, (key, value)| {
                        acc.entry(key).or_insert_with(Vec::new).push(value);
                        acc
                    });

                if set.is_empty() {
                    return None;
                } else {
                    return Some(set);
                }
            }
        }
    }

    /// Runs binary classification for a list of ranks.
    ///
    /// # Arguments
    ///
    /// * `gold_std` - Gold standard profile for comparison
    /// * `ranks` - Ranks to evaluate
    /// * `allow_ambiguity` - Whether to allow alternative name matching
    ///
    /// # Returns
    ///
    /// DataFrame combining TP/FP/FN rows across ranks.
    ///
    /// # Errors
    ///
    /// Returns a Polars error when DataFrame operations fail.
    pub fn binary_classification(
        &self,
        gold_std: &Self,
        ranks: &[TaxonomicRank],
        allow_ambiguity: bool,
    ) -> PolarsResult<DataFrame> {
        let mut dfs = Vec::default();
        let mut dfs_alt = Vec::default();

        for rank in ranks {
            debug!("----------{}-----------", rank.to_string());
            let prediction_names = self.get_taxa_string_dict(rank);
            let gold_std_names = gold_std
                .get_taxa_string_dict(rank)
                .map(filter_prefix_only_taxa_names);

            if let (Some(prediction), Some(gold_std)) = (prediction_names, gold_std_names) {
                match binary_classification_df(&prediction, &gold_std) {
                    Ok(mut df) => {
                        let _ = df.with_column(Series::new(
                            "Rank".into(),
                            std::iter::repeat_n(rank.to_owned(), df.height())
                                .map(|x| x.to_string())
                                .collect::<Vec<_>>(),
                        ));
                        dfs.push(df)
                    }
                    Err(_) => {
                        debug!("Nothing");
                    }
                }
            }

            let prediction_taxa = self.get_taxa_dict(rank);
            let gold_std_taxa = gold_std
                .get_taxa_dict(rank)
                .map(filter_prefix_only_taxa_refs);

            if allow_ambiguity {
                if let (Some(prediction), Some(gold_std)) = (prediction_taxa, gold_std_taxa) {
                    match binary_classification_alternatives_df(&prediction, &gold_std) {
                        Ok(mut df) => {
                            let _ = add_string_columns(
                                &mut df,
                                &[("Rank".to_string(), rank.to_string())],
                            );
                            dfs_alt.push(df);
                        }
                        Err(_) => {
                            debug!("Nothing");
                        }
                    }
                }
            }
        }

        let mut joined_df = dfs
            .into_iter()
            .reduce(|df1, df2| df1.vstack(&df2).unwrap())
            .unwrap();
        joined_df.with_column(Series::new(
            "AllowAlternatives".into(),
            std::iter::repeat(false)
                .take(joined_df.height())
                .collect_vec(),
        ))?;

        if allow_ambiguity {
            let mut joined_df_alt = dfs_alt
                .into_iter()
                .reduce(|df1, df2| df1.vstack(&df2).unwrap())
                .unwrap();
            joined_df_alt.with_column(Series::new(
                "AllowAlternatives".into(),
                std::iter::repeat(true)
                    .take(joined_df.height())
                    .collect_vec(),
            ))?;
            joined_df.vstack(&joined_df_alt)
        } else {
            Ok(joined_df)
        }
    }
}

/// Type-erased wrapper around profiles of different taxonomies.
#[derive(Debug, Clone)]
pub enum ProfileWrapper {
    GTDBProfile(Profile<GTDB>),
    NCBIProfile(Profile<NCBI>),
    CustomProfile(Profile<Custom>),
}

impl ProfileWrapper {
    /// Returns the taxonomy identifier for the wrapped profile.
    ///
    /// # Returns
    ///
    /// `TaxonomyEnum` describing the profile taxonomy.
    pub fn taxonomy(&self) -> TaxonomyEnum {
        match self {
            ProfileWrapper::GTDBProfile(_) => TaxonomyEnum::GTDB,
            ProfileWrapper::NCBIProfile(_) => TaxonomyEnum::NCBI,
            ProfileWrapper::CustomProfile(_) => TaxonomyEnum::Custom,
        }
    }

    /// Returns taxa names at a given rank for the wrapped profile.
    ///
    /// # Arguments
    ///
    /// * `rank` - Target rank to extract
    ///
    /// # Returns
    ///
    /// Set of taxon names, or `None` if unavailable.
    pub fn taxa(&self, rank: &TaxonomicRank) -> Option<HashSet<String>> {
        match self {
            ProfileWrapper::GTDBProfile(profile) => profile.get_taxa_with_rank(rank),
            ProfileWrapper::NCBIProfile(profile) => profile.get_taxa_with_rank(rank),
            ProfileWrapper::CustomProfile(profile) => profile.get_taxa_with_rank(rank),
        }
    }

    /// Returns the set of ranks present in the wrapped profile.
    ///
    /// # Returns
    ///
    /// `None` if ranks are unknown; otherwise a set of ranks.
    pub fn unique_ranks(&self) -> Option<HashSet<TaxonomicRank>> {
        match self {
            ProfileWrapper::GTDBProfile(profile) => profile.unique_ranks(),
            ProfileWrapper::NCBIProfile(profile) => profile.unique_ranks(),
            ProfileWrapper::CustomProfile(profile) => profile.unique_ranks(),
        }
    }

    /// Runs binary classification between matching taxonomy profiles.
    ///
    /// # Arguments
    ///
    /// * `gold_std` - Gold standard profile with the same taxonomy
    /// * `allow_alternatives` - Whether to consider alternative names
    ///
    /// # Returns
    ///
    /// DataFrame containing TP/FP/FN rows for each rank.
    ///
    /// # Errors
    ///
    /// Returns a Polars error when DataFrame operations fail.
    pub fn binary_classification(
        &self,
        gold_std: &ProfileWrapper,
        allow_alternatives: bool,
    ) -> PolarsResult<DataFrame> {
        match (self, gold_std) {
            (ProfileWrapper::GTDBProfile(pred), ProfileWrapper::GTDBProfile(gold)) => {
                pred.binary_classification(gold, &TaxonomicRank::all(), allow_alternatives)
            }
            (ProfileWrapper::NCBIProfile(pred), ProfileWrapper::NCBIProfile(gold)) => {
                pred.binary_classification(gold, &TaxonomicRank::all(), allow_alternatives)
            }
            (ProfileWrapper::CustomProfile(pred), ProfileWrapper::CustomProfile(gold)) => {
                pred.binary_classification(gold, &TaxonomicRank::all(), allow_alternatives)
            }
            (_, _) => panic!("Invalid types"),
        }
    }
}

impl ProfileWrapper {
    /// Normalizes abundances with sample-wide consistency checks.
    ///
    /// # Arguments
    ///
    /// * `tolerance` - Relative tolerance for deciding scale
    ///
    /// # Returns
    ///
    /// Report containing pre/post rank sums and whether normalization was applied.
    ///
    /// # Errors
    ///
    /// Returns `ProfileError::NormalizationError` when sums are inconsistent or invalid.
    pub fn normalize_abundances_with_checks(
        &mut self,
        tolerance: f64,
    ) -> Result<NormalizeReport, ProfileError> {
        match self {
            ProfileWrapper::GTDBProfile(profile) => {
                profile.normalize_abundances_with_checks(tolerance)
            }
            ProfileWrapper::NCBIProfile(profile) => {
                profile.normalize_abundances_with_checks(tolerance)
            }
            ProfileWrapper::CustomProfile(profile) => {
                profile.normalize_abundances_with_checks(tolerance)
            }
        }
    }

    /// Returns the number of unique ranks present in the profile.
    ///
    /// # Returns
    ///
    /// Number of ranks or 0 when ranks are unknown.
    pub fn unique_rank_count(&self) -> usize {
        match self {
            ProfileWrapper::GTDBProfile(profile) => {
                profile.unique_ranks().map(|r| r.len()).unwrap_or(0)
            }
            ProfileWrapper::NCBIProfile(profile) => {
                profile.unique_ranks().map(|r| r.len()).unwrap_or(0)
            }
            ProfileWrapper::CustomProfile(profile) => {
                profile.unique_ranks().map(|r| r.len()).unwrap_or(0)
            }
        }
    }
}

/// Errors returned while loading or validating profiles.
#[derive(thiserror::Error, Debug, Clone)]
pub enum ProfileError {
    #[error("{0}")]
    GenericError(String),
    #[error("{0}")]
    CamiFormatError(String),
    #[error("{0}")]
    FormatError(String),
    #[error("{0}")]
    TaxonomyError(String),
    #[error("{0}")]
    NormalizationError(String),
}

pub type ProfileResult<T: Taxonomy> = Result<Profile<T>, ProfileError>;

/// Loads a profile from a reader using a format-specific parser.
pub trait LoadProfile<T: Taxonomy + Default> {
    /// Parse the provided reader using the given format and columns.
    ///
    /// # Arguments
    ///
    /// * `input` - Reader positioned at the start of the profile
    /// * `columns` - Optional column definition override
    ///
    /// # Returns
    ///
    /// Parsed profile with taxonomy-aware entries.
    ///
    /// # Errors
    ///
    /// Returns `ProfileError` when parsing fails.
    fn load<F: Format, R: Read + Seek>(input: &mut R, columns: Option<Columns>)
        -> ProfileResult<T>;
}

// Implement LoadProfile for Profile
impl<T: Taxonomy + Default> LoadProfile<T> for Profile<T> {
    fn load<F: Format, R: Read + Seek>(
        input: &mut R,
        columns: Option<Columns>,
    ) -> ProfileResult<T> {
        F::load_profile(input, columns)
    }
}

#[cfg(test)]
mod tests {
    use std::{collections::HashMap, io::Cursor};

    use crate::{
        common::{LineageFromString, Taxon, TaxonomicRank, GTDB},
        format::{Auto, Columns, CAMI},
    };

    use super::*;

    #[test]
    fn test_prefix_only_taxon_name() {
        assert!(is_prefix_only_taxon_name("s__"));
        assert!(is_prefix_only_taxon_name("c__ "));
        assert!(is_prefix_only_taxon_name("d__\t"));
        assert!(!is_prefix_only_taxon_name("s__Lactobacillus"));
        assert!(!is_prefix_only_taxon_name("Bacteria"));
    }

    #[test]
    fn test_filter_prefix_only_taxa_refs() {
        let entry = Entry::<GTDB> {
            taxon_name: None,
            lineage: None,
            alternative_names: None,
            abundance: 0.0,
            rank: TaxonomicRank::Species,
        };

        let prefix_taxon = Taxon::with_name_and_rank("s__", &TaxonomicRank::Species);
        let named_taxon = Taxon::with_name_and_rank("s__Lactobacillus", &TaxonomicRank::Species);

        let taxa = HashMap::from([(&prefix_taxon, vec![&entry]), (&named_taxon, vec![&entry])]);

        let filtered = filter_prefix_only_taxa_refs(taxa);

        assert!(!filtered.contains_key(&prefix_taxon));
        assert!(filtered.contains_key(&named_taxon));
    }

    #[test]
    fn test_filter_prefix_only_taxa_names() {
        let entry = Entry::<GTDB> {
            taxon_name: None,
            lineage: None,
            alternative_names: None,
            abundance: 0.0,
            rank: TaxonomicRank::Species,
        };

        let taxa = HashMap::from([
            ("s__".to_string(), vec![&entry]),
            ("s__Lactobacillus".to_string(), vec![&entry]),
        ]);

        let filtered = filter_prefix_only_taxa_names(taxa);

        assert!(!filtered.contains_key("s__"));
        assert!(filtered.contains_key("s__Lactobacillus"));
    }

    #[test]
    fn test_normalize_abundances_if_percent_single_rank() {
        let mut profile = Profile::<GTDB>::default();
        profile.taxa = vec![
            Entry {
                abundance: 50.0,
                rank: TaxonomicRank::Species,
                ..Entry::default()
            },
            Entry {
                abundance: 50.0,
                rank: TaxonomicRank::Species,
                ..Entry::default()
            },
        ];

        let report = profile
            .normalize_abundances_with_checks(ABUNDANCE_SUM_TOLERANCE)
            .expect("Normalization should succeed");

        assert!(report.applied);
        assert!((profile.taxa[0].abundance - 0.5).abs() < 1e-6);
        assert!((profile.taxa[1].abundance - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_normalize_abundances_if_percent_multi_rank_total_over_100() {
        let mut profile = Profile::<GTDB>::default();
        profile.taxa = vec![
            Entry {
                abundance: 100.0,
                rank: TaxonomicRank::Domain,
                ..Entry::default()
            },
            Entry {
                abundance: 60.0,
                rank: TaxonomicRank::Species,
                ..Entry::default()
            },
            Entry {
                abundance: 40.0,
                rank: TaxonomicRank::Species,
                ..Entry::default()
            },
        ];

        let report = profile
            .normalize_abundances_with_checks(ABUNDANCE_SUM_TOLERANCE)
            .expect("Normalization should succeed");

        assert!(report.applied);
        assert!((profile.taxa[0].abundance - 1.0).abs() < 1e-6);
        assert!((profile.taxa[1].abundance - 0.6).abs() < 1e-6);
        assert!((profile.taxa[2].abundance - 0.4).abs() < 1e-6);
    }

    #[test]
    fn test_normalize_abundances_if_percent_no_change_for_fractional() {
        let mut profile = Profile::<GTDB>::default();
        profile.taxa = vec![
            Entry {
                abundance: 0.5,
                rank: TaxonomicRank::Species,
                ..Entry::default()
            },
            Entry {
                abundance: 0.5,
                rank: TaxonomicRank::Species,
                ..Entry::default()
            },
        ];

        let report = profile
            .normalize_abundances_with_checks(ABUNDANCE_SUM_TOLERANCE)
            .expect("Normalization should succeed");

        assert!(!report.applied);
        assert!((profile.taxa[0].abundance - 0.5).abs() < 1e-6);
        assert!((profile.taxa[1].abundance - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_normalize_abundances_with_checks_mixed_scales_error() {
        let mut profile = Profile::<GTDB>::default();
        profile.taxa = vec![
            Entry {
                abundance: 100.0,
                rank: TaxonomicRank::Domain,
                ..Entry::default()
            },
            Entry {
                abundance: 0.9,
                rank: TaxonomicRank::Species,
                ..Entry::default()
            },
        ];

        let result = profile.normalize_abundances_with_checks(ABUNDANCE_SUM_TOLERANCE);

        assert!(result.is_err());
    }

    #[test]
    fn test_normalize_abundances_ignores_small_rank_sums() {
        let mut profile = Profile::<GTDB>::default();
        profile.taxa = vec![
            Entry {
                abundance: 0.99,
                rank: TaxonomicRank::Species,
                ..Entry::default()
            },
            Entry {
                abundance: 0.01,
                rank: TaxonomicRank::Order,
                ..Entry::default()
            },
        ];

        let report = profile
            .normalize_abundances_with_checks(ABUNDANCE_SUM_TOLERANCE)
            .expect("Normalization should succeed");

        assert!(!report.applied);
    }

    #[test]
    fn test_rank_sum_violations_detects_excess() {
        let mut profile = Profile::<GTDB>::default();
        profile.taxa = vec![
            Entry {
                abundance: 0.6,
                rank: TaxonomicRank::Species,
                ..Entry::default()
            },
            Entry {
                abundance: 0.6,
                rank: TaxonomicRank::Species,
                ..Entry::default()
            },
        ];

        let violations = profile.rank_sum_violations(1.0001);

        assert_eq!(violations.len(), 1);
        assert_eq!(violations[0].0, TaxonomicRank::Species);
        assert!((violations[0].1 - 1.2).abs() < 1e-6);
    }

    #[test]
    fn test_rank_sum_violations_allows_within_threshold() {
        let mut profile = Profile::<GTDB>::default();
        profile.taxa = vec![
            Entry {
                abundance: 0.4,
                rank: TaxonomicRank::Species,
                ..Entry::default()
            },
            Entry {
                abundance: 0.6,
                rank: TaxonomicRank::Species,
                ..Entry::default()
            },
        ];

        let violations = profile.rank_sum_violations(1.0001);

        assert!(violations.is_empty());
    }

    #[test]
    fn test_get_taxa_with_rank_from_species_only_profile() {
        let lineage_str = "d__Bacteria;p__Verrucomicrobiota;c__Verrucomicrobiae;\
o__Verrucomicrobiales;f__Akkermansiaceae;g__Akkermansia;\
s__Akkermansia muciniphila_A";
        let lineage = GTDB::lineage_from_string(lineage_str, None);

        let entry = Entry::<GTDB> {
            taxon_name: None,
            lineage: Some(lineage),
            alternative_names: None,
            abundance: 0.1,
            rank: TaxonomicRank::Species,
        };

        let profile = Profile::<GTDB> {
            unstructured_meta: Vec::new(),
            meta: HashMap::new(),
            taxa: vec![entry],
        };

        let genus = profile.get_taxa_with_rank(&TaxonomicRank::Genus);
        assert!(genus.is_some());
        assert!(genus.unwrap().contains("g__Akkermansia"));
    }

    #[test]
    fn test_get_taxa_string_dict_from_species_only_profile() {
        let lineage_str = "d__Bacteria;p__Verrucomicrobiota;c__Verrucomicrobiae;\
o__Verrucomicrobiales;f__Akkermansiaceae;g__Akkermansia;\
s__Akkermansia muciniphila_A";
        let lineage = GTDB::lineage_from_string(lineage_str, None);

        let entry = Entry::<GTDB> {
            taxon_name: None,
            lineage: Some(lineage),
            alternative_names: None,
            abundance: 0.1,
            rank: TaxonomicRank::Species,
        };

        let profile = Profile::<GTDB> {
            unstructured_meta: Vec::new(),
            meta: HashMap::new(),
            taxa: vec![entry],
        };

        let genus_map = profile.get_taxa_string_dict(&TaxonomicRank::Genus);
        assert!(genus_map.is_some());
        assert!(genus_map.unwrap().contains_key("g__Akkermansia"));
    }

    #[test]
    fn test_load_profile_cami() {
        let mut test_data = Cursor::new("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n");

        let profile = Profile::<GTDB>::default();
        let result = Profile::<GTDB>::load::<CAMI, _>(&mut test_data, None);

        assert!(result.is_ok());
    }

    #[test]
    fn test_auto_columns() {
        let mut test_data1 = Cursor::new("RS_GCF_001436455.1	d__Bacteria|p__Bacillota|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae|g__Lactobacillus|s__Lactobacillus jensenii	0.00113805");
        let mut test_data2 = Cursor::new("d__Bacteria|p__Pseudomonadota|c__Gammaproteobacteria|o__Pseudomonadales	0.7119	1.8877000000000002	NA	NA");

        let columns1 = Auto::derive_columns(&mut test_data1);
        let columns2 = Auto::derive_columns(&mut test_data2);
        debug!("Columns1: {:?}", columns1);
        debug!("Columns2: {:?}", columns2);

        assert!(columns1.is_some());
        assert!(columns2.is_none());
    }

    #[test]
    fn test_is_gtdb() {
        let test1 = "d__Bacteria";
        let test2 = "d__Bacteria|p__Pseudomonadota|c__Gammaproteobacteria|o__Pseudomonadales";
        let test3 = "12332|123|213";

        let result1 = Columns::is_gtdb_lineage(test1);
        let result2 = Columns::is_gtdb_lineage(test2);
        let result3 = Columns::is_gtdb_lineage(test3);

        assert_eq!(result1, true);
        assert_eq!(result2, true);
        assert_eq!(result3, false);
    }
}
