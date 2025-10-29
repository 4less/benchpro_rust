use std::{cmp::max, collections::{HashMap, HashSet}, default, fmt::Display, io::{BufRead, BufReader, Read, Seek}, process::exit, str::FromStr};

use itertools::Itertools;
use polars::{error::PolarsResult, frame::DataFrame, prelude::NamedFrom, series::Series};

use crate::{common::{Custom, Lineage, LineageFromString, Taxon, TaxonomicRank, Taxonomy, TaxonomyEnum, GTDB, NCBI}, format::{Columns, Format}, utils::add_string_columns};



pub type StringList = Vec<String>;
pub type StringStringMap = HashMap<String, String>;
pub type Entries<T: Taxonomy> = Vec<Entry<T>>;
pub type EntriesRef<'a, T: Taxonomy> = Vec<&'a Entry<T>>;


#[derive(Debug, Default, Clone, PartialEq)]
pub struct Entry<T: Taxonomy> {
    pub taxon_name: Option<String>,
    pub lineage: Option<Lineage<T>>,
    pub alternative_names: Option<Vec<String>>,
    pub abundance: f64,
    pub rank: TaxonomicRank,
}

impl<T: Taxonomy> Entry<T> {
    pub fn lineage(&self) -> Option<&Lineage<T>> {
        self.lineage.as_ref()
    }
}

#[derive(Default, Clone, Debug)]
struct BCVectors {
    pred_type: Vec::<String>,
    taxon_names: Vec::<String>,
    taxon_ids: Vec::<usize>,
    prediction_abundances: Vec::<f64>,
    gold_std_abundances: Vec::<f64>,
    prediction_count: Vec::<usize>,
    gold_std_count: Vec::<usize>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum BC {
    TP, FP, FN, TN, FFP, FFN, Unknown
}

impl Default for BC {
    fn default() -> Self {
        BC::Unknown
    }
}

impl BC {
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
            _ => Err(format!("BC Type does not exist: {}", value))
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
            gold_std_count: usize) {

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
            // Series::new("Rank".into(), rank_vec.iter().map(|x| x.to_string()).collect::<Vec<_>>()),
            Series::new("Name".into(), self.taxon_names),
            Series::new("Type".into(), self.pred_type),
            Series::new("PredictionAbundance".into(), self.prediction_abundances),
            Series::new("GoldStdAbundance".into(), self.gold_std_abundances),
            Series::new("PredictionCount".into(), self.prediction_count.iter().map(|&x| x as i64).collect::<Vec<_>>()),
            Series::new("GoldStdCount".into(), self.gold_std_count.iter().map(|&x| x as i64).collect::<Vec<_>>()),
        ]);
        df
    }

    pub fn polars_df(&self) -> PolarsResult<DataFrame> {
        let df = DataFrame::new(vec![
            // Series::new("Rank".into(), rank_vec.iter().map(|x| x.to_string()).collect::<Vec<_>>()),
            Series::new("Name".into(), self.taxon_names.clone()),
            Series::new("Type".into(), self.pred_type.clone()),
            Series::new("PredictionAbundance".into(), self.prediction_abundances.clone()),
            Series::new("GoldStdAbundance".into(), self.gold_std_abundances.clone()),
            Series::new("PredictionCount".into(), self.prediction_count.clone().iter().map(|&x| x as i64).collect::<Vec<_>>()),
            Series::new("GoldStdCount".into(), self.gold_std_count.clone().iter().map(|&x| x as i64).collect::<Vec<_>>()),
        ]);
        df
    }
}

type TaxonEntryMap<'a, T: Taxonomy> = HashMap<&'a Taxon, EntriesRef<'a, T>>;
type NamesEntryMap<'a, T: Taxonomy> = HashMap<String, EntriesRef<'a, T>>;

// fn get_binary_classification_alternative_names<T: Taxonomy>(
//         vectors: &mut BCVectors, 
//         prediction: TaxonEntryMap<T>,
//         gold_std: TaxonEntryMap<T>,
//     ) -> PolarsResult<DataFrame> 
// {

// }

pub fn binary_classification_df<T: Taxonomy>(
    prediction: &NamesEntryMap<T>,
    gold_std: &NamesEntryMap<T>,
) -> PolarsResult<DataFrame> 
{
    let mut df_vectors = BCVectors::default();

    let pred_set = prediction.keys().collect::<HashSet<_>>();
    let gold_set = gold_std.keys().collect::<HashSet<_>>();

    let tps = pred_set.intersection(&gold_set).map(|&x| x).collect::<HashSet<_>>();
    let fps = pred_set.difference(&gold_set).map(|&x| x).collect::<HashSet<_>>();
    let fns = gold_set.difference(&pred_set).map(|&x| x).collect::<HashSet<_>>();

    // eprintln!("Prediction set:");

    // for (name, entries) in prediction {
    //     println!("---{}", name);
    //     for entry in entries {
    //         println!("{} - {:?}", entry.abundance, entry.lineage().unwrap().to_string());
    //     }
    // }

    tps.iter().map(|&t| (BC::TP, t))
        .chain(fps.iter().map(|&t| (BC::FP, t)))
        .chain(fns.iter().map(|&t| (BC::FN, t)))
        .for_each(|(mut bp_type, taxon)| 
    {
        // match (prediction.get(taxon), gold_std.get(taxon) {

        // }

        let (prediction_abundance, prediction_entry_count) = if matches!(bp_type, BC::TP) || matches!(bp_type, BC::FP) {
            let entries = prediction.get(taxon)
                 .expect(&format!("Gold set does not contain TP taxon. Unrecoverable error. '{}', {}", taxon, bp_type.to_string()));
            (
                 entries
                     .iter().map(|e| e.abundance)
                     .sum(),
                 entries.len()
            )
        } else { (0f64, 0usize) };

        let (gold_std_abundance, gold_std_entry_count) = if matches!(bp_type, BC::TP) || matches!(bp_type, BC::FN) {
           let entries = gold_std.get(taxon)
           .expect(&format!("Gold set does not contain FN taxon. Unrecoverable error. '{}', {}", taxon, bp_type.to_string()));
           (
                entries
                    .iter().map(|e| e.abundance)
                    .sum(),
                entries.len()
           )
        } else { (0f64, 0usize) };

        // println!(">> {} {}", prediction_abundance, gold_std_abundance);



        if taxon == "" && matches!(bp_type, BC::FP) {
            eprintln!("FIX BC TYPE");
            bp_type = BC::Unknown;
        }

        df_vectors.add_entry(
            &bp_type,
            taxon, 
            0, 
            prediction_abundance, 
            gold_std_abundance, 
            prediction_entry_count, 
            gold_std_entry_count);
    });
    df_vectors.into_polars_df()
}


pub fn binary_classification_alternatives_df<T: Taxonomy>(
    prediction: &TaxonEntryMap<T>,
    gold_std: &TaxonEntryMap<T>,
) -> PolarsResult<DataFrame> 
{
    // Differences to Name based classification:
    // - some Taxa have alternative names e.g. (SpeciesA/SpeciesB/SpeciesC)
    // New information needed:
    // - Does an alternative names entry match with multiple entries in the Gold STd?

    let mut df_vectors = BCVectors::default();

    // let pred_set = prediction.keys().collect::<HashSet<_>>();
    // let gold_set = gold_std.keys().collect::<HashSet<_>>();

    // let tps = pred_set.intersection(&gold_set).map(|&x| x).collect::<HashSet<_>>();
    // let fps = pred_set.difference(&gold_set).map(|&x| x).collect::<HashSet<_>>();
    // let fns = gold_set.difference(&pred_set).map(|&x| x).collect::<HashSet<_>>();

    let prediction_taxa = prediction.keys().collect::<HashSet<_>>();
    let gold_std_taxa = gold_std.keys().collect::<HashSet<_>>();
    
    let BCData {
        tps, 
        fps,
        fns,
        prediction_counts
    } = binary_classification_helper(&prediction_taxa, &gold_std_taxa);

    if tps.iter().any(|(_,v)| {
        v.secondary.is_some() && v.primary.is_none()
    }) {
        for (taxon,value) in tps.iter() {
            eprintln!("Gold {:?}\n\t\tPred {:?}", taxon, value);
        }

        // exit(9);
    }

    tps.iter().for_each(|(&taxon, predictions)| {
        // println!("TP: {:?}", taxon);

        let (prediction_abundance, prediction_entry_count) = {
            // if primary is set, primary abundance always counts.
            // if secondary is set, only secondary entries count that only match ONE gold std entry
            // Ambiguous matches with more than one gold_std are not counted.

            let primary_abundance_list = predictions.primary.as_ref().map_or(vec![], |x| 
                prediction.get(taxon).unwrap().iter().map(|t| t.abundance).collect_vec());
            let primary_abundance: f64 = primary_abundance_list.iter().sum();
            let primary_entry_count = primary_abundance_list.len();

            let unique_secondary_abundance_list = predictions.secondary.as_ref().map_or(vec![], |x| 
                x.iter()
                    .filter(|&&x| prediction_counts.get(x).expect("asdads") == &1)
                    .flat_map(|t| prediction.get(t).unwrap())
                    .map(|&e| e.abundance)
                    .collect_vec()
            );

            let non_unique_secondary_abundance_list = predictions.secondary.as_ref().map_or(vec![], |x| 
                x.iter()
                    .filter(|&&x| prediction_counts.get(x).expect("asdads") > &1)
                    .flat_map(|t| prediction.get(t).unwrap())
                    .map(|&e| e.abundance)
                    .collect_vec()
            );
            let unique_secondary_abundance: f64 = unique_secondary_abundance_list.iter().sum();
            let non_unique_secondary_abundance: f64 = non_unique_secondary_abundance_list.iter().sum();
            let unique_secondary_count = unique_secondary_abundance_list.len();
            let non_unique_secondary_count = non_unique_secondary_abundance_list.len();

            if non_unique_secondary_abundance != 0.0f64 {

                println!("---- {:?}", taxon);
                println!("Abundance: Primary {}  Secondary (U) {} Secondary (NU) {}", primary_abundance, unique_secondary_abundance, non_unique_secondary_abundance);
                println!("Counts:    Primary {}  Secondary (U) {} Secondary (NU) {}", primary_entry_count, unique_secondary_count, non_unique_secondary_count);

                for g in gold_std {
                    eprintln!("GOLD: {:?}", g.0);
                }
                for p in prediction {
                    eprintln!("PRED: {:?}", p.0);
                }
            }

            (primary_abundance + unique_secondary_abundance,
            primary_entry_count + unique_secondary_count)
        };
        
        let (gold_std_abundance, gold_std_entry_count) = {
            let entries = match gold_std.get(taxon) {
                Some(entries) => entries,
                None => {
                    println!("--------------------------------");
                    gold_std.iter().for_each(|(t, _)| {
                        println!("Gold: {} {}", t.match_exact(taxon), t.match_any(taxon));
                    });
                    panic!("Gold set does not contain FN taxon. Unrecoverable error. '{:?}', TP", taxon)
                },
            };
            (
                entries
                    .iter().map(|e| e.abundance)
                    .sum(),
                entries.len()
            )
        };


        df_vectors.add_entry(
            &BC::TP,
            &taxon.name, 
            0, 
            prediction_abundance, 
            gold_std_abundance, 
            prediction_entry_count, 
            gold_std_entry_count);
    });

    fps.iter().map(|&t| (BC::FP, t))
        .chain(fns.iter().map(|&t| (BC::FN, t)))
        .for_each(|(mut bp_type, taxon)| 
    {
        println!("{}: {:?}", bp_type.to_string(), taxon);

        let (prediction_abundance, prediction_entry_count) = if matches!(bp_type, BC::TP) || matches!(bp_type, BC::FP) {
            
            let entries = match prediction.get(taxon) {
                Some(entries) => entries,
                None => {
                    println!("--------------------------------");
                    gold_std.iter().for_each(|(t, _)| {
                        println!("Gold: {} {}", t.match_exact(taxon), t.match_any(taxon));
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
                    println!("--------------------------------");
                    gold_std.iter().for_each(|(t, _)| {
                        println!("Gold: {} {}", t.match_exact(taxon), t.match_any(taxon));
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
            eprintln!("FIX BC TYPE");
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
#[derive(Default)]
pub struct BCData<'a> {
    tps: TPMap<'a>,
    fps: OtherSet<'a>,
    fns: OtherSet<'a>,
    prediction_counts: MatchCountMap<'a>,
}

pub fn binary_classification_helper<'a, TREF: AsRef<Taxon>>(
    prediction: &'a HashSet<TREF>,
    gold_std: &'a HashSet<TREF>,
) -> BCData<'a>
{
    let mut fn_return = BCData::default();

    let mut matched_gold_stds = HashSet::new();
    let mut matched_predictions = HashSet::new();


    let all_gold_no_alt = gold_std.iter().all(|x| x.as_ref().alternative_names.is_none());
    
    if !all_gold_no_alt {
        for gs in gold_std.iter() {
            println!("{:?}", gs.as_ref());
        }
        
        panic!("Gold standard may not have any alternative names.")
    }

    // iterator over all pairs of gold std and prediction
    let pair_iter = 
        gold_std.iter()
        .flat_map(|gold| {
            prediction.iter().map(|pred| (gold.as_ref(), pred.as_ref()))
        });

    // iterator over exact match pairs.
    let exact_match_pairs = pair_iter.clone()
        .filter(|(gold, prediction)| {
            assert!(gold.alternative_names.is_none());
            gold.match_exact(prediction)
        });
    
    // iterator over alternative name match pairs.
    let alternative_match_pairs = pair_iter
        .filter(|(gold, prediction)| {
            assert!(gold.alternative_names.is_none());
            gold.match_any(prediction) && !gold.match_exact(prediction)
        });

    // Extract exact matches first
    exact_match_pairs.clone()
        .for_each(|(gold, prediction)| {
            assert!(gold.alternative_names.is_none());
            if matched_gold_stds.contains(gold) {
                eprintln!("Entry: {:?}", matched_gold_stds.get(gold).unwrap());
                eprintln!("GoldEntry: {:?}", fn_return.tps.get(gold).unwrap());
                eprintln!("Pred: {:?}", prediction);

                exact_match_pairs.clone().for_each(|(g,p)| {
                    println!("----MATCH\nG: {:?}\nP: {:?}", g, p);
                });
            }
            assert!(!matched_gold_stds.contains(gold));
            assert!(!fn_return.tps.contains_key(gold));
            matched_gold_stds.insert(gold);

            matched_predictions.insert(prediction);
            fn_return.tps.insert(gold, TPMapValue {
                primary: Some(gold),
                secondary: None,
            });
        });

    alternative_match_pairs
        .for_each(|(gold, prediction)| {
            let entry = fn_return.tps
                .entry(gold).or_insert(TPMapValue::default());

            matched_predictions.insert(prediction);
            if let Some(list) = &mut entry.secondary {
                list.push(prediction);
            } else {
                entry.secondary = Some(vec![prediction]);
            }

            fn_return.prediction_counts
                .entry(prediction).and_modify(|v| *v += 1).or_insert(1); 
        });

    fn_return.fps = prediction.iter()
        .filter(|p| !matched_predictions.contains(p.as_ref()))
        .map(|x| x.as_ref())
        .collect::<HashSet<_>>();


    fn_return.fns = gold_std.iter()
        .filter(|g| !matched_gold_stds.contains(g.as_ref()))
        .map(|x| x.as_ref())
        .collect::<HashSet<_>>();
    
    fn_return
}


#[derive(Debug, Default, Clone, PartialEq)]
pub struct Profile<T: Taxonomy> {
        pub unstructured_meta: StringList,
        pub meta: StringStringMap,
        pub taxa: Entries<T>,
}

impl Profile<GTDB> {
    pub fn wrap(self) -> ProfileWrapper {
        ProfileWrapper::GTDBProfile(self)
    }
}

impl Profile<NCBI> {
    pub fn wrap(self) -> ProfileWrapper {
        ProfileWrapper::NCBIProfile(self)
    }
}

impl Profile<Custom> {
    pub fn wrap(self) -> ProfileWrapper {
        ProfileWrapper::CustomProfile(self)
    }
}

impl<T: Taxonomy> Profile<T> {
    pub fn unique_ranks(&self) -> Option<HashSet::<TaxonomicRank>> {



        if self.taxa.iter().all(|entry| matches!(entry.rank, TaxonomicRank::Unknown)) {
            return None
        }
        assert!(self.taxa.iter().all(|entry| !matches!(entry.rank, TaxonomicRank::Unknown)));

        Some(self.taxa.iter().map(|e| e.rank.clone()).collect::<HashSet<_>>())
    }

    /// If ranks are defined and there is more than one rank, 
    /// I deduce that ranks are defined separately on purpose.
    /// If ranks are undefined, see if the lineage contains the target ranks. 
    /// - This only works for the GTDB taxonomy. 
    /// - For the NCBI taxonomy, something more elaborate needss to be figured out
    pub fn get_taxa_with_rank(&self, rank: &TaxonomicRank) -> Option<HashSet<String>> {
        match self.unique_ranks() {
            Some(ranks) => {
                if (ranks.len() > 1 && !ranks.contains(rank)) || rank > ranks.iter().max().expect("Has no max") {
                    return None
                }
                if ranks.len() == 1 && rank <= ranks.iter().max().expect("Has no max") {
                    let set = self.taxa.iter()
                        .map(|entry| entry.lineage().unwrap().get(rank))
                        .filter(|name| name.is_some())
                        .map(|taxon| taxon.unwrap())
                        .map(|taxon| taxon.name.clone())
                        .collect::<HashSet<_>>();
                }

                let set = self.taxa.iter()
                    .filter(|entry| matches!(&entry.rank, rank))
                    .map(|entry| entry.lineage().unwrap().get(rank))
                    .filter(|name| name.is_some())
                    .map(|taxon| taxon.unwrap())
                    .map(|taxon| taxon.name.clone())
                    .collect::<HashSet<_>>();

                if set.is_empty() {
                    return None
                } else {
                    return Some(set)
                }
            },
            None => {
                let set = self.taxa.iter()
                    .map(|entry| entry.lineage().unwrap().get(rank))
                    .filter(|name| name.is_some())
                    .map(|taxon| taxon.unwrap())
                    .map(|taxon| taxon.name.clone())
                    .collect::<HashSet<_>>();
                
                if set.is_empty() {
                    return None
                } else {
                    return Some(set)
                }
            },
        }
    }

    pub fn get_taxa_string_dict(&self, rank: &TaxonomicRank) -> Option<HashMap<String, Vec<&Entry<T>>>> {
        match self.unique_ranks() {
            Some(ranks) => {
                if (ranks.len() > 1 && !ranks.contains(rank)) || rank > ranks.iter().max().expect("Has no max") {
                    return None
                }
                if ranks.len() == 1 && rank <= ranks.iter().max().expect("Has no max") {
                    println!("Option 1: There is only one rank defined {:?}", ranks);
                    let set = self.taxa.iter()
                        .map(|entry| (entry.lineage().unwrap().get(rank), entry))
                        .filter(|(name, entry)| name.is_some())
                        .map(|(taxon, entry)| (taxon.unwrap(), entry))
                        .map(|(taxon, entry)| (taxon.name.clone(), entry))
                        .fold(HashMap::new(), |mut acc, (key, value)| {
                            acc.entry(key).or_insert_with(Vec::new).push(value);
                            acc
                        });

                    if set.is_empty() {
                        return None
                    } else {
                        return Some(set).clone()
                    }
                }
                println!("Option 2: There are multiple ranks defined {:?} (Target: {:?})", ranks, rank);
                let set = self.taxa.iter()
                    .filter(|&entry| rank == &entry.rank)
                    .map(|entry| (entry.lineage().unwrap().get(rank), entry))
                    .filter(|(name, _)| name.is_some())
                    .map(|(taxon, entry)| (taxon.unwrap(), entry))
                    .map(|(taxon, entry)| (taxon.name.clone(), entry))
                    .fold(HashMap::new(), |mut acc, (key, value)| {
                        acc.entry(key).or_insert_with(Vec::new).push(value);
                        acc
                    });
                    
                // set.iter().flat_map(|(name, entries)| entries)//.all(|entry|)
                //     .for_each(|entry| println!("{}", entry.rank.to_string()));

                if set.is_empty() {
                    return None
                } else {
                    return Some(set)
                }
            },
            None => {
                println!("Option3");
                let set = self.taxa.iter()
                    .map(|entry| (entry.lineage().unwrap().get(rank), entry))
                    .filter(|(name, entry)| name.is_some())
                    .map(|(taxon, entry)| (taxon.unwrap(), entry))
                    .map(|(taxon, entry)| (taxon.name.clone(), entry))
                    .fold(HashMap::new(), |mut acc, (key, value)| {
                        acc.entry(key).or_insert_with(Vec::new).push(value);
                        acc
                    });
                
                if set.is_empty() {
                    return None
                } else {
                    return Some(set)
                }
            },
        }
    }


    pub fn get_taxa_dict<'a>(&'a self, rank: &TaxonomicRank) -> Option<HashMap<&'a Taxon, Vec<&Entry<T>>>> {
        match self.unique_ranks() {
            Some(ranks) => {
                if (ranks.len() > 1 && !ranks.contains(rank)) || rank > ranks.iter().max().expect("Has no max") {
                    return None
                }
                if ranks.len() == 1 && rank <= ranks.iter().max().expect("Has no max") {
                    // println!("Option 1: There is only one rank defined {:?}", ranks);
                    let set = self.taxa.iter()
                        .map(|entry| (entry.lineage().unwrap().get(rank), entry))
                        .filter(|(name, _)| name.is_some())
                        .map(|(taxon, entry)| (taxon.unwrap(), entry))
                        .fold(HashMap::new(), |mut acc, (key, value)| {
                            acc.entry(key).or_insert_with(Vec::new).push(value);
                            acc
                        });

                    if set.is_empty() {
                        return None
                    } else {
                        return Some(set)
                    }
                }
                // println!("Option 2: There are multiple ranks defined {:?}", ranks);
                let set = self.taxa.iter()
                    .filter(|&entry| &entry.rank == rank)
                    .map(|entry| (entry.lineage().unwrap().get(rank), entry))
                    .filter(|(name, _)| name.is_some())
                    .map(|(taxon, entry)| (taxon.unwrap(), entry))
                    .fold(HashMap::new(), |mut acc, (key, value)| {
                        acc.entry(key).or_insert_with(Vec::new).push(value);
                        acc
                    });
                    

                if set.is_empty() {
                    return None
                } else {
                    return Some(set)
                }
            },
            None => {
                // println!("Option3");
                let set = self.taxa.iter()
                    .map(|entry| (entry.lineage().unwrap().get(rank), entry))
                    .filter(|(name, _)| name.is_some())
                    .map(|(taxon, entry)| (taxon.unwrap(), entry))
                    .fold(HashMap::new(), |mut acc, (key, value)| {
                        acc.entry(key).or_insert_with(Vec::new).push(value);
                        acc
                    });
                
                if set.is_empty() {
                    return None
                } else {
                    return Some(set)
                }
            },
        }
    }


    

    pub fn binary_classification(&self, gold_std: &Self, ranks: &[TaxonomicRank], allow_ambiguity: bool) -> PolarsResult<DataFrame> {
        let mut dfs = Vec::default();
        let mut dfs_alt = Vec::default();
        
        for rank in ranks {
            println!("----------{}-----------", rank.to_string());
            let prediction_names = self.get_taxa_string_dict(rank);
            let gold_std_names = gold_std.get_taxa_string_dict(rank);


            if let (Some(prediction), Some(gold_std)) = (prediction_names, gold_std_names) {
                match binary_classification_df(&prediction, &gold_std) {
                    Ok(mut df) => {

                        let _ = df.with_column(Series::new("Rank".into(), std::iter::repeat_n(rank.to_owned(), df.height()).map(|x| x.to_string()).collect::<Vec<_>>()));
                        dfs.push(df)
                    },
                    Err(_) => println!("Nothing"),
                }
            }


            let prediction_taxa = self.get_taxa_dict(rank);
            let gold_std_taxa = gold_std.get_taxa_dict(rank);

            if let (Some(prediction), Some(gold_std)) = (prediction_taxa, gold_std_taxa) {
                match binary_classification_alternatives_df(&prediction, &gold_std) {
                    Ok(mut df) => {
                        let _ = add_string_columns(&mut df, 
                            &[("Rank".to_string(), rank.to_string())]);
                        dfs_alt.push(df);
                    },
                    Err(_) => println!("Nothing"),
                }
            }
        }

        let mut joined_df = dfs.into_iter().reduce(|df1, df2| df1.vstack(&df2).unwrap()).unwrap();
        joined_df.with_column(Series::new("AllowAlternatives".into(), std::iter::repeat(false).take(joined_df.height()).collect_vec()))?;

        let mut joined_df_alt= dfs_alt.into_iter().reduce(|df1, df2| df1.vstack(&df2).unwrap()).unwrap();
        joined_df_alt.with_column(Series::new("AllowAlternatives".into(), std::iter::repeat(true).take(joined_df.height()).collect_vec()))?;


        joined_df.vstack(&joined_df_alt)
    }
}

#[derive(Debug)]
pub enum ProfileWrapper {
    GTDBProfile(Profile<GTDB>),
    NCBIProfile(Profile<NCBI>),
    CustomProfile(Profile<Custom>)
}

impl ProfileWrapper {
    pub fn taxonomy(&self) -> TaxonomyEnum {
        match self {
            ProfileWrapper::GTDBProfile(_) => TaxonomyEnum::GTDB,
            ProfileWrapper::NCBIProfile(_) => TaxonomyEnum::NCBI,
            ProfileWrapper::CustomProfile(_) => TaxonomyEnum::Custom,
        }
    }

    pub fn taxa(&self, rank: &TaxonomicRank) -> Option<HashSet<String>> {
        match self {
            ProfileWrapper::GTDBProfile(profile) => profile.get_taxa_with_rank(rank),
            ProfileWrapper::NCBIProfile(profile) => profile.get_taxa_with_rank(rank),
            ProfileWrapper::CustomProfile(profile) => profile.get_taxa_with_rank(rank),
        }
    }

    pub fn unique_ranks(&self) -> Option<HashSet<TaxonomicRank>> {
        match self {
            ProfileWrapper::GTDBProfile(profile) => profile.unique_ranks(),
            ProfileWrapper::NCBIProfile(profile) => profile.unique_ranks(),
            ProfileWrapper::CustomProfile(profile) => profile.unique_ranks(),
        }
    }

    pub fn binary_classification(&self, gold_std: &ProfileWrapper) -> PolarsResult<DataFrame> {
        match (self, gold_std) {
            (ProfileWrapper::GTDBProfile(pred), ProfileWrapper::GTDBProfile(gold)) => pred.binary_classification(gold, &TaxonomicRank::all(), false),
            (ProfileWrapper::NCBIProfile(pred), ProfileWrapper::NCBIProfile(gold)) => pred.binary_classification(gold, &TaxonomicRank::all(), false),
            (ProfileWrapper::CustomProfile(pred), ProfileWrapper::CustomProfile(gold)) => pred.binary_classification(gold, &TaxonomicRank::all(), false),
            (_, _) => panic!("Invalid types"),
        }
    }
}

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
}

pub type ProfileResult<T: Taxonomy> = Result<Profile<T>, ProfileError>;


// LoadProfile trait for dynamic loading
pub trait LoadProfile<T: Taxonomy + Default> {
    fn load<F: Format, R: Read + Seek>(input: &mut R, columns: Option<Columns>) -> ProfileResult<T>;
}

// Implement LoadProfile for Profile
impl<T: Taxonomy + Default> LoadProfile<T> for Profile<T> {
    fn load<F: Format, R: Read + Seek>(input: &mut R, columns: Option<Columns>) -> ProfileResult<T> {
        F::load_profile(input, columns)
    }
}



#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use crate::{common::{self, Taxonomy, GTDB}, format::{Auto, Columns, CAMI}, profile};

    use super::*;

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
        
        eprintln!("Columns1: {:?}", columns1);
        eprintln!("Columns2: {:?}", columns2);


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
