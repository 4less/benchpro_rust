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
use polars::{
    error::PolarsResult,
    frame::DataFrame,
    io::SerWriter,
    prelude::{
        col, cols, lit, when, CsvWriter, DataFrameJoinOps, DataType, IntoLazy, NamedFrom,
        PlSmallStr, QuoteStyle,
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

pub fn run(args: &Args) {
    match &args.meta {
        Some(_) => meta_based_workflow(args),
        None => meta_free_workflow(args),
    }
}

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
        .expect("Unsuccesful");

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

pub fn get_taxon_df(handler: &ProfileHandler) -> DataFrame {
    println!("Number of profiles: {}", handler.prediction_map.len());
    println!("Number of gold profiles: {}", handler.gold_std_map.len());

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
            println!(">>>>>>>>>>>> {:?}", row);
        }

        if let (Some(prediction), Some(goldstd)) = (prediction, goldstd) {
            // println!("START---{}---{}-\n{}\n{}\n---------------------------", row.id, row.taxonomy, row.profile.to_str().unwrap(), row.goldstd.to_str().unwrap());
            
            println!("Dataset: {} ... {:?}", row.id, row.goldstd);
            let df = prediction.binary_classification(goldstd);

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
                    // let string_vec: Vec<String> = ;
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

    let mut complete_df = get_taxon_df(&handler);

    

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

    let safe_th = Mutex::new(tree_handler);

    let df_adjusted = get_adjusted_benchmarks(&complete_df, &handler.meta, &safe_th, &handler);
    eprintln!("DF Adjusted {:?}", df_adjusted);

    println!("Height {}", new_complete_df.height());
    match df_adjusted {
        Ok(mut dfa) => {
            println!("Height DFA {}", dfa.height());
            dfa = handler
                .meta
                .left_join_to(&dfa, &[MetaColumn::Dataset, MetaColumn::Taxonomy, MetaColumn::Tool], true)
                .expect("Cannot join dfs?");

            println!("Height DFA {}", dfa.height());

            add_string_columns(&mut dfa, &[("Adjusted".to_string(), "True".to_string())]);

            println!("Height DFA {}", dfa.height());
            let _ = add_string_columns(
                &mut new_complete_df,
                &[
                    ("ClosestNeighbor", ""),
                    ("ClosestNeighborType", ""),
                    ("Adjusted", "False"),
                ]
                .map(|(a, b)| (a.to_string(), b.to_string())),
            );
            new_complete_df.with_column(Series::new("ClosestNeighborDistance".into(), std::iter::repeat(0f64).take(new_complete_df.height()).collect_vec())).unwrap();
            new_complete_df.with_column(Series::new("ClosestNeighborAbundance".into(), std::iter::repeat(0f64).take(new_complete_df.height()).collect_vec())).unwrap();


            let reorder = new_complete_df.get_column_names().into_iter().map(|x| x.to_owned());
            let dfa = dfa.select(reorder).unwrap();

            new_complete_df = new_complete_df.vstack(&dfa).unwrap();

            // let file_name_detailed_adj = format!("{}_detailed_adj.tsv", args.outprefix);
            // let file_path = Path::new(&file_name_detailed_adj);

            // // Open the file for writing
            // let mut file = File::create(file_path).expect("Could not create file");

            // // Write the DataFrame to the CSV file using CsvWriter
            // let _ = CsvWriter::new(&mut file)
            //     .include_header(true)
            //     .with_separator(b'\t')
            //     .with_quote_style(polars::prelude::QuoteStyle::Never)
            //     .finish(&mut dfa);
            // println!("File written to: {}", file_name_detailed_adj);
        }
        Err(e) => println!("Disappointing {:?}", e),
    }
    

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

    println!("NEWDF {}", newdf);
    println!("File written to: {}", file_name_detailed);
    println!("File written to: {}", file_name_bc);
}

pub fn meta_free_workflow(args: &Args) {}
