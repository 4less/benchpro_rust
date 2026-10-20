use core::f64;
use std::{
    collections::{HashMap, HashSet},
    str::FromStr,
    sync::Mutex,
};

use itertools::Itertools;
use log::debug;
use polars::{
    error::PolarsResult,
    frame::DataFrame,
    prelude::{col, lit, DataFrameJoinOps, IntoLazy, NamedFrom as _},
    series::Series,
};

use crate::{
    meta::{Meta, MetaColumn},
    profile::BC,
    profile_handler::ProfileHandler,
    taxonomic_profiling::detailed_data::DetailedData,
    tree_handler::{TaxaSet, TreeHandler},
    utils::{closest_neighbor, sample_apply, time, wrap_names},
};

/// Reclassifies TP/FP/FN using phylogenetic proximity.
///
/// # Arguments
///
/// * `data` - Detailed benchmark DataFrame to adjust
/// * `meta` - Parsed meta table for tree references
/// * `tree_handler` - Cached tree handler for subtree extraction
/// * `profile_handler` - Loaded profiles and meta context
///
/// # Returns
///
/// DataFrame with adjusted classification labels.
///
/// # Errors
///
/// Returns a Polars error when DataFrame operations fail.
pub fn get_adjusted_benchmarks(
    data: &DataFrame,
    _meta: &Meta,
    tree_handler: &Mutex<TreeHandler>,
    profile_handler: &ProfileHandler,
) -> PolarsResult<DataFrame> {
    let threshold = 0.04_f64;

    debug!("--FIND--");
    debug!("{:?}", data);

    // Filter for only species and remove allow alternatives
    let mut_species_df: DataFrame = data
        .clone()
        .lazy()
        .filter(
            col("Rank")
                .eq(lit("Species"))
                .and(col("AllowAlternatives").eq(lit(false))),
        )
        .collect()?;

    // TREEEEEE
    let bc_tree_df = mut_species_df
        .left_join(
            &profile_handler
                .meta
                .raw
                .select(["ID", "GoldStdTree"])
                .expect("Cannot subset df"),
            ["ID"],
            ["ID"],
        )
        .expect("Cannot join with meta-data");

    let mut taxon_to_type: HashMap<String, BC> = HashMap::new();
    let mut taxon_to_detectable: HashMap<String, bool> = HashMap::new();

    // TODO: With tree and dataset, now relabel FP and FN
    // Distance between nodes, sort by distance, resolve pairings..
    let mut all_dfs = sample_apply(&bc_tree_df, |df| {
        let mut rows = DetailedData::from_polars_df(&df).unwrap();

        let tree_path = df
            .column(MetaColumn::GoldStdTree.to_str())
            .unwrap()
            .str()
            .unwrap()
            .first()
            .expect("Must have GoldStdTree set");

        let names = df
            .column("Name")
            .unwrap()
            .str()
            .unwrap()
            .iter()
            .map(|x| x.unwrap().to_owned())
            .collect::<TaxaSet>();

        let th_locked = tree_handler.lock().unwrap();

        let (duration, result) = time(|| th_locked.get_subtree(tree_path, &names));

        debug!("Get Subtree took {:?}", duration);

        if let Some(tree) = result {
            if tree.n_leaves() < 2 {
                return Ok(df);
            };

            let mut print_tree = tree.clone();
            let _ = wrap_names(&mut print_tree);
            taxon_to_type.clear();
            taxon_to_detectable.clear();

            // Taxon to Type hashmap
            df.column("Name")
                .unwrap()
                .str()
                .unwrap()
                .iter()
                .zip(df.column("Type").unwrap().str().unwrap().iter())
                .map(|(a, b)| (a.unwrap().to_owned(), BC::from_str(b.unwrap()).unwrap()))
                .collect_into(&mut taxon_to_type);

            // Taxon to Detectable hashmap
            df.column("Name")
                .unwrap()
                .str()
                .unwrap()
                .iter()
                .zip(df.column("DetectableTaxon").unwrap().str().unwrap().iter())
                .map(|(a, b)| (a.unwrap().to_owned(), b.unwrap() == "True"))
                .collect_into(&mut taxon_to_detectable);

            let taxon_to_abundance: HashMap<String, f64> = rows
                .data()
                .iter()
                .map(|r| {
                    // FP rows carry a gold abundance of 0.0, so pick by type, not presence.
                    let abundance = if r.bc_type == BC::FP {
                        r.prediction_abundance
                    } else {
                        r.goldstd_abundance
                    };
                    (r.name.clone(), abundance.unwrap_or(0f64))
                })
                .collect::<HashMap<_, _>>();

            for row in rows.data_mut().iter_mut() {
                let node = tree.get_by_name(&row.name);

                if let Some(node) = node {
                    let dist = closest_neighbor(&tree, &node.id).unwrap();
                    row.closest_neigbor = Some(dist.name.clone());
                    row.closest_neighbor_distance = Some(dist.distance);
                    row.closest_neighbor_type = taxon_to_type.get(&dist.name).map(|x| x.to_owned());
                    row.closest_neighbor_detectable =
                        taxon_to_detectable.get(&dist.name).map(|x| x.to_owned());
                    row.closest_neighbor_abundance = taxon_to_abundance
                        .get(dist.name.as_str())
                        .map(|x| x.to_owned());
                }
            }

            rows.sort_by(
                |a, b| match (a.closest_neighbor_distance, b.closest_neighbor_distance) {
                    (Some(a_val), Some(b_val)) => a_val
                        .partial_cmp(&b_val)
                        .unwrap_or(std::cmp::Ordering::Equal),
                    (None, None) => std::cmp::Ordering::Equal,
                    (None, Some(_)) => std::cmp::Ordering::Greater,
                    (Some(_), None) => std::cmp::Ordering::Less,
                },
            );

            let mut used_fns = HashSet::new();
            let mut fn_to_fp = Vec::new();

            // Collect necessary changes without mutating rows directly
            rows.data_mut().iter_mut().for_each(|row| {
                if let (Some(name), Some(dist), Some(type_), Some(detectable)) = (
                    &row.closest_neigbor,
                    row.closest_neighbor_distance,
                    &row.closest_neighbor_type,
                    row.closest_neighbor_detectable,
                ) {
                    if row.bc_type == BC::FP && dist <= threshold && type_ == &BC::FN && !detectable
                    {
                        row.bc_type = BC::FFP;
                        if !used_fns.contains(name) {
                            fn_to_fp.push((name.clone(), row.name.clone())); // Clone values for later use
                            used_fns.insert(name.clone());
                        }
                    }
                }
            });

            // Apply the collected changes
            for (fn_, fp) in fn_to_fp {
                if let Some(row) = rows.get_mut(&fn_) {
                    row.bc_type = BC::TP;
                    row.prediction_abundance = Some(*taxon_to_abundance.get(fp.as_str()).unwrap());
                }
            }
        };

        rows.to_polars_df()
    });
    if let Ok(df) = all_dfs.as_mut() {
        let new_col = std::iter::repeat(false).take(df.height()).collect_vec();
        let new_series = Series::new("AllowAlternatives".into(), new_col);
        df.with_column(new_series)?;
    }

    all_dfs
}

#[cfg(test)]
mod tests {
    use std::{collections::HashMap, fs, sync::Mutex};

    use phylotree::tree::Tree;
    use polars::prelude::{DataFrame, NamedFrom, Series};

    use super::get_adjusted_benchmarks;
    use crate::{
        meta::Meta,
        profile_handler::ProfileHandler,
        tree_handler::{Name2Id, TreeHandler},
    };

    /// Pairs at distance 0.02: (fn1, fp1) with fn1 undetectable, (fn3, fp3) with fn3
    /// detectable. fp2 and fn2 have no neighbor within the 0.04 threshold.
    const NEWICK: &str = "((s__fn1:0.01,s__fp1:0.01):0.5,(s__tp:0.3,s__fp2:0.3):0.5,\
(s__fn3:0.01,s__fp3:0.01):0.5,s__fn2:1.0);";

    fn tree_handler_for(path: &str) -> TreeHandler {
        let tree = Tree::from_newick(NEWICK).unwrap();
        let name2id = tree
            .get_leaves()
            .iter()
            .map(|id| (tree.get(id).unwrap().name.clone().unwrap(), *id))
            .collect::<Name2Id>();
        let mut tree_map = HashMap::new();
        tree_map.insert(path.to_string(), (name2id, tree));
        TreeHandler { tree_map }
    }

    fn profile_handler_for(path: &str) -> ProfileHandler {
        let raw = DataFrame::new(vec![
            Series::new("ID".into(), &["S1"]),
            Series::new("GoldStdTree".into(), &[path]),
        ])
        .unwrap();
        ProfileHandler {
            meta: Meta {
                raw,
                entries: Vec::new(),
            },
            prediction_map: HashMap::new(),
            gold_std_map: HashMap::new(),
            available_taxa: HashMap::new(),
        }
    }

    fn detailed_df() -> DataFrame {
        let rows: &[(&str, &str, f64, f64, &str)] = &[
            ("s__tp", "TP", 0.2, 0.2, "True"),
            ("s__fp1", "FP", 0.3, 0.0, "False"),
            ("s__fn1", "FN", 0.0, 0.25, "False"),
            ("s__fp2", "FP", 0.1, 0.0, "False"),
            ("s__fn2", "FN", 0.0, 0.15, "False"),
            ("s__fp3", "FP", 0.4, 0.0, "False"),
            ("s__fn3", "FN", 0.0, 0.4, "True"),
        ];
        let n = rows.len();
        DataFrame::new(vec![
            Series::new("Name".into(), rows.iter().map(|r| r.0).collect::<Vec<_>>()),
            Series::new("Type".into(), rows.iter().map(|r| r.1).collect::<Vec<_>>()),
            Series::new(
                "PredictionAbundance".into(),
                rows.iter().map(|r| r.2).collect::<Vec<_>>(),
            ),
            Series::new(
                "GoldStdAbundance".into(),
                rows.iter().map(|r| r.3).collect::<Vec<_>>(),
            ),
            Series::new(
                "PredictionCount".into(),
                rows.iter().map(|r| (r.2 > 0.0) as i64).collect::<Vec<_>>(),
            ),
            Series::new(
                "GoldStdCount".into(),
                rows.iter().map(|r| (r.3 > 0.0) as i64).collect::<Vec<_>>(),
            ),
            Series::new(
                "DetectableTaxon".into(),
                rows.iter().map(|r| r.4).collect::<Vec<_>>(),
            ),
            Series::new("ValidTaxon".into(), vec!["false"; n]),
            Series::new("Rank".into(), vec!["Species"; n]),
            Series::new("ID".into(), vec!["S1"; n]),
            Series::new("AllowAlternatives".into(), vec![false; n]),
        ])
        .unwrap()
    }

    #[test]
    fn test_adjusted_reclassifies_fp_close_to_undetectable_fn() {
        let dir = std::env::temp_dir().join(format!("benchpro_adjusted_{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let tree_path = dir.join("tree.nwk");
        fs::write(&tree_path, NEWICK).unwrap();
        let tree_path = tree_path.to_str().unwrap();

        let profile_handler = profile_handler_for(tree_path);
        let tree_handler = Mutex::new(tree_handler_for(tree_path));

        let df = get_adjusted_benchmarks(
            &detailed_df(),
            &profile_handler.meta,
            &tree_handler,
            &profile_handler,
        )
        .unwrap();
        fs::remove_dir_all(&dir).unwrap();

        let names = df.column("Name").unwrap().str().unwrap();
        let types = df.column("Type").unwrap().str().unwrap();
        let pred = df.column("PredictionAbundance").unwrap().f64().unwrap();
        let gold = df.column("GoldStdAbundance").unwrap().f64().unwrap();
        let rows: HashMap<&str, (&str, f64, f64)> = (0..df.height())
            .map(|i| {
                let value = (
                    types.get(i).unwrap(),
                    pred.get(i).unwrap(),
                    gold.get(i).unwrap(),
                );
                (names.get(i).unwrap(), value)
            })
            .collect();

        assert_eq!(rows.len(), 7);
        assert_eq!(rows["s__tp"].0, "TP");

        // FP within threshold of an undetectable FN: the FP becomes FFP and the FN becomes
        // a TP carrying the FP's abundance.
        assert_eq!(rows["s__fp1"].0, "FFP");
        assert_eq!(rows["s__fn1"], ("TP", 0.3, 0.25));

        // No neighbor within the threshold: unchanged.
        assert_eq!(rows["s__fp2"].0, "FP");
        assert_eq!(rows["s__fn2"].0, "FN");

        // FN is detectable by the tool, so a nearby FP is a genuine error.
        assert_eq!(rows["s__fp3"].0, "FP");
        assert_eq!(rows["s__fn3"].0, "FN");

        let alt = df.column("AllowAlternatives").unwrap().bool().unwrap();
        assert!(alt.into_iter().all(|v| v == Some(false)));
    }
}
