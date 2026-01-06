use std::{
    collections::HashSet,
    fs::File,
    io::{BufReader, Cursor},
    iter::repeat,
    mem::swap,
    path::Path,
    process::exit,
    sync::{Arc, Mutex, RwLock},
};

use itertools::{izip, Itertools};
use log::{debug, info, warn};
use polars::{
    error::PolarsResult,
    frame::DataFrame,
    io::SerWriter,
    prelude::{
        col, cols, lit, when, CsvWriter, DataFrameJoinOps, DataType, IntoLazy, NamedFrom,
        PlSmallStr, QuoteStyle, SortMultipleOptions,
    },
    series::Series,
};

use crate::{
    common::{Detectable, TaxonomicRank},
    meta::{self, MetaColumn},
    options::Args,
    profile::{LoadProfile, Profile, ProfileWrapper},
    profile_handler::ProfileHandler,
    tree_adjusted_benchmarks::get_adjusted_benchmarks,
    tree_handler::{TaxaSet, TreeHandler},
    utils::{
        add_string_columns, f1_score, get_lca, get_subtree, get_subtree_with_leaves, precision,
        sample_apply, sensitivity, tree_collapse_edges, wrap_names, write_df,
    },
};

/// Run Benchpro using the provided CLI arguments.
///
/// # Arguments
///
/// * `args` - Parsed CLI arguments
pub fn run(args: &Args) {
    match &args.meta {
        Some(_) => meta_based_workflow(args),
        None => meta_free_workflow(args),
    }
}

/// Aggregate TP/FP/FN counts per group and compute summary metrics.
///
/// # Arguments
///
/// * `df` - Detailed per-taxon DataFrame
///
/// # Returns
///
/// DataFrame with aggregated counts and summary metrics.
///
/// # Errors
///
/// Returns a Polars error when aggregation fails.
pub fn add_binary_classification(df: DataFrame) -> PolarsResult<DataFrame> {
    let mut newdf = df
        .lazy()
        .group_by(["Rank", "ID", "AllowAlternatives", "Adjusted"])
        .agg([
            col("Type")
                .filter(col("Type").eq(lit("TP")))
                .count()
                .alias("TP"),
            col("Type")
                .filter(col("Type").eq(lit("FP")))
                .count()
                .alias("FP"),
            col("Type")
                .filter(col("Type").eq(lit("FN")))
                .count()
                .alias("FN"),
        ])
        .collect()
        .expect("Unsuccessful");

    let series = newdf.columns(["TP", "FP", "FN"]).unwrap();

    let iter = izip!(
        series[0].u32().unwrap().iter(),
        series[1].u32().unwrap().iter(),
        series[2].u32().unwrap().iter()
    )
    .map(|(tp, fp, fn_)| (tp.unwrap(), fp.unwrap(), fn_.unwrap()))
    .collect_vec();

    let f1_series = Series::new(
        "F1".into(),
        iter.iter()
            .map(|(tp, fp, fn_)| f1_score(*tp as usize, *fp as usize, *fn_ as usize))
            .collect_vec(),
    );
    let sensitivity_series = Series::new(
        "Sensitivity".into(),
        &iter
            .iter()
            .map(|(tp, fp, fn_)| sensitivity(*tp as usize, *fp as usize, *fn_ as usize))
            .collect_vec(),
    );
    let precision_series = Series::new(
        "Precision".into(),
        iter.iter()
            .map(|(tp, fp, fn_)| precision(*tp as usize, *fp as usize, *fn_ as usize))
            .collect_vec(),
    );

    let _ = newdf.with_column(f1_series)?;
    let _ = newdf.with_column(sensitivity_series)?;
    let _ = newdf.with_column(precision_series)?;

    Ok(newdf)
}

/// Build the per-taxon binary classification DataFrame for all meta entries.
///
/// # Arguments
///
/// * `handler` - Profile handler with loaded profiles and meta
/// * `allow_alternatives` - Whether to allow alternative name matching
///
/// # Returns
///
/// DataFrame with per-taxon classification rows.
pub fn get_taxon_df(handler: &ProfileHandler, allow_alternatives: bool) -> DataFrame {
    debug!("Number of profiles: {}", handler.prediction_map.len());
    debug!("Number of gold profiles: {}", handler.gold_std_map.len());

    let mut dfs: Vec<DataFrame> = Vec::default();
    let mut dfs_with_taxa: Vec<DataFrame> = Vec::default();

    for row in handler.meta.entries.iter() {
        let prediction = handler.prediction_map.get(row.profile.to_str().unwrap());
        let goldstd = handler.gold_std_map.get(row.goldstd.to_str().unwrap());

        let taxa = row.taxa_list.as_ref().map(|x| {
            handler
                .available_taxa
                .get(x.to_str().unwrap())
                .expect("Should be loaded...")
        });

        if prediction.is_some() && !goldstd.is_some() {
            warn!("Missing gold standard for row: {:?}", row);
        }

        if let (Some(prediction), Some(goldstd)) = (prediction, goldstd) {
            debug!("Dataset: {} ... {:?}", row.id, row.goldstd);
            let df = prediction.binary_classification(goldstd, allow_alternatives);

            if let Ok(mut df) = df {
                let _ = df.with_column(Series::new(
                    "ID".into(),
                    repeat(row.id.clone()).take(df.height()).collect::<Vec<_>>(),
                ));
                let _ = df.with_column(Series::new(
                    "Tool".into(),
                    repeat(row.tool.clone())
                        .take(df.height())
                        .collect::<Vec<_>>(),
                ));
                let _ = df.with_column(Series::new(
                    "Dataset".into(),
                    repeat(row.dataset.clone())
                        .take(df.height())
                        .collect::<Vec<_>>(),
                ));
                let _ = df.with_column(Series::new(
                    "ValidTaxon".into(),
                    repeat("false").take(df.height()).collect::<Vec<_>>(),
                ));
                let _ = df.with_column(Series::new(
                    "DetectableTaxon".into(),
                    repeat(Detectable::Unknown.to_string())
                        .take(df.height())
                        .collect::<Vec<_>>(),
                ));

                if let Some(available_taxa) = taxa.as_ref() {
                    let detectable_series = Series::new(
                        "detectable".into(),
                        available_taxa
                            .iter()
                            .map(|x| x.to_owned())
                            .collect::<Vec<_>>(),
                    );

                    // Modify the Detectable column
                    df = df
                        .lazy()
                        .with_column(
                            when(col("Rank").eq(lit(TaxonomicRank::Species.to_string())))
                                .then(
                                    when(col("Name").is_in(lit(detectable_series)))
                                        .then(lit("True"))
                                        .otherwise(lit("False")),
                                )
                                .otherwise(col("DetectableTaxon")) // Keep original value
                                .alias("DetectableTaxon"),
                        )
                        .collect()
                        .unwrap();
                }

                dfs.push(df);
            }
        }
    }
    let complete_df = if dfs.len() > 1 {
        dfs.into_iter()
            .reduce(|df1, df2| df1.vstack(&df2).unwrap())
            .unwrap()
    } else {
        dfs.first().unwrap().clone()
    };

    return complete_df;
}

/// Execute the meta-driven workflow and write detailed and summary outputs.
/// Execute the meta-driven benchmark workflow and write outputs.
///
/// # Arguments
///
/// * `args` - Parsed CLI arguments
pub fn meta_based_workflow(args: &Args) {
    type C = MetaColumn;
    let path = Path::new(args.meta.as_ref().unwrap());

    let handler = match ProfileHandler::from_meta(path) {
        Ok(handler) => handler,
        Err(e) => panic!("Failed to load profiles: {}", e),
    };

    let tree_handler = TreeHandler::from_meta(&handler.meta);

    let tree_handler = match tree_handler {
        Ok(th) => th,
        Err(e) => panic!("{}", e),
    };

    ///////////////////////////////////////////////////////////////////////
    // Binary classification DF

    let mut complete_df = get_taxon_df(&handler, args.allow_alternatives);

    

    let mut new_complete_df = complete_df
        .left_join(
            &handler
                .meta
                .raw
                .select(["ID", "Taxonomy"])
                .expect("Cannot subset df"),
            ["ID"],
            ["ID"],
        )
        .expect("Cannot join with meta-data");

    let _ = add_string_columns(
        &mut new_complete_df,
        &[("Adjusted".to_string(), "False".to_string())],
    );

    let safe_th = Mutex::new(tree_handler);

    if args.adjusted {
        let df_adjusted = get_adjusted_benchmarks(&complete_df, &handler.meta, &safe_th, &handler);
        debug!("DF Adjusted {:?}", df_adjusted);
        debug!("Height {}", new_complete_df.height());
        match df_adjusted {
            Ok(mut dfa) => {
                debug!("Height DFA {}", dfa.height());
                dfa = handler
                    .meta
                    .left_join_to(&dfa, &[MetaColumn::Dataset, MetaColumn::Taxonomy, MetaColumn::Tool], true)
                    .expect("Cannot join dfs?");

                debug!("Height DFA {}", dfa.height());

                add_string_columns(&mut dfa, &[("Adjusted".to_string(), "True".to_string())]);

                debug!("Height DFA {}", dfa.height());
                let _ = add_string_columns(
                    &mut new_complete_df,
                    &[
                        ("ClosestNeighbor", ""),
                        ("ClosestNeighborType", ""),
                    ]
                    .map(|(a, b)| (a.to_string(), b.to_string())),
                );
                new_complete_df.with_column(Series::new("ClosestNeighborDistance".into(), std::iter::repeat(0f64).take(new_complete_df.height()).collect_vec())).unwrap();
                new_complete_df.with_column(Series::new("ClosestNeighborAbundance".into(), std::iter::repeat(0f64).take(new_complete_df.height()).collect_vec())).unwrap();


                let reorder = new_complete_df.get_column_names().into_iter().map(|x| x.to_owned());
                let dfa = dfa.select(reorder).unwrap();

                new_complete_df = new_complete_df.vstack(&dfa).unwrap();

        }
            Err(e) => {
                warn!("Adjusted benchmarks failed: {:?}", e);
            },
        }
    }
    

    new_complete_df = new_complete_df
        .sort(
            ["ID", "Rank", "Name", "Type", "AllowAlternatives", "Adjusted"],
            SortMultipleOptions::default(),
        )
        .expect("Cannot sort detailed output");

    // Specify the path where you want to save the CSV file

    let file_name_detailed = format!("{}_detailed.tsv", args.outprefix);
    let file_path = Path::new(&file_name_detailed);

    // Open the file for writing
    let mut file = File::create(file_path).expect("Could not create file");

    // Write the DataFrame to the CSV file using CsvWriter
    let _ = CsvWriter::new(&mut file)
        .include_header(true)
        .with_quote_style(QuoteStyle::Never)
        .with_separator(b'\t')
        .finish(&mut new_complete_df);

    ///////////////////////////////////////////////////////////////////////
    // Summary Binary classification DF

    let mut newdf =
        add_binary_classification(new_complete_df.clone()).expect("Could not derive binclas");

    let mut newdf = newdf
        .left_join(&handler.meta.raw, ["ID"], ["ID"])
        .expect("Cannot join with meta-data");

    newdf = newdf
        .sort(
            ["ID", "Rank", "AllowAlternatives", "Adjusted"],
            SortMultipleOptions::default(),
        )
        .expect("Cannot sort summary output");

    // Open the file for writing
    let file_name_bc = format!("{}.tsv", args.outprefix);
    let file_path: &Path = Path::new(&file_name_bc);
    let mut file = File::create(file_path).expect("Could not create file");

    // Write the DataFrame to the CSV file using CsvWriter
    let _ = CsvWriter::new(&mut file)
        .include_header(true)
        .with_separator(b'\t')
        .with_quote_style(QuoteStyle::Never)
        .finish(&mut newdf);

    debug!("NEWDF {}", newdf);
    info!("File written to: {}", file_name_detailed);
    info!("File written to: {}", file_name_bc);
}

/// Placeholder for a meta-free workflow.
/// Execute the meta-free workflow (not implemented).
///
/// # Arguments
///
/// * `args` - Parsed CLI arguments
pub fn meta_free_workflow(args: &Args) {}
