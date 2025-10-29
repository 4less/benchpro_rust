use core::f64;
use std::{collections::{HashMap, HashSet}, process::exit, str::FromStr, sync::{Arc, Mutex}};

use itertools::Itertools;
use polars::{
    error::PolarsResult,
    frame::DataFrame,
    prelude::{all, col, lit, DataFrameJoinOps, IntoLazy, NamedFrom as _, NamedFromOwned},
    series::Series,
};

use crate::{
    meta::{Meta, MetaColumn},
    profile::BC,
    profile_handler::ProfileHandler,
    taxonomic_profiling::detailed_data::DetailedData,
    tree_handler::{TaxaSet, TreeHandler},
    utils::{closest_neighbor, sample_apply, time, wrap_names, NeighborDist},
};

pub fn get_adjusted_benchmarks(
    data: &DataFrame,
    meta: &Meta,
    tree_handler: &Mutex<TreeHandler>,
    profile_handler: &ProfileHandler,
) -> PolarsResult<DataFrame> {
    let threshold =  0.04_f64;

    println!("--FIND--");
    println!("{:?}", data);

    // Filter for only species and remove allow alternatives
    let mut_species_df: DataFrame = data
        .clone()
        .lazy()
        .filter(col("Rank").eq(lit("Species")).and(col("AllowAlternatives").eq(lit(false))))
        .collect()?;


    // TREEEEEE
    let mut bc_tree_df = mut_species_df
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

        println!("Get Subtree took {:?}", duration);

        if let Some(tree) = result {
            if tree.n_leaves() < 2 {
                return Ok(df);
            };

            let mut print_tree = tree.clone();
            let _ = wrap_names(&mut print_tree);

            // println!("Tree:\n: {}", print_tree.to_newick().unwrap());

            // for name in &names {
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


            let taxon_to_abundance: HashMap<String, f64> = rows.data().iter().map(|r| {
                let abundance: f64 = match (r.goldstd_abundance, r.prediction_abundance) {
                    (None, None) => 0f64,
                    (None, Some(p)) => p,
                    (Some(g), None) => g,
                    (Some(g), Some(_)) => g,
                };
                (r.name.clone(), abundance)
            }).collect::<HashMap<_, _>>();

            for row in rows.data_mut().iter_mut() {
                let node = tree.get_by_name(&row.name);


                if let Some(node) = node {
                    let dist = closest_neighbor(&tree, &node.id).unwrap();
                    row.closest_neigbor = Some(dist.name.clone());
                    row.closest_neighbor_distance = Some(dist.distance);
                    row.closest_neighbor_type = taxon_to_type.get(&dist.name).map(|x| x.to_owned());
                    row.closest_neighbor_detectable =
                        taxon_to_detectable.get(&dist.name).map(|x| x.to_owned());
                    row.closest_neighbor_abundance = taxon_to_abundance.get(dist.name.as_str()).map(|x| x.to_owned());
                }
            }

            rows.sort_by(|a, b| {
                match (a.closest_neighbor_distance, b.closest_neighbor_distance) {
                    (Some(a_val), Some(b_val)) => a_val.partial_cmp(&b_val).unwrap_or(std::cmp::Ordering::Equal),
                    (None, None) => std::cmp::Ordering::Equal,
                    (None, Some(_)) => std::cmp::Ordering::Greater,
                    (Some(_), None) => std::cmp::Ordering::Less,
                }
            });
            
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
                    if row.bc_type == BC::FP && dist <= threshold && type_ == &BC::FN && !detectable {
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
