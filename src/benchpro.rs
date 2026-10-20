use std::{
    collections::{HashMap, HashSet},
    fs::File,
    iter::repeat,
    path::Path,
    process::exit,
    sync::Mutex,
};

use itertools::{izip, Itertools};
use log::{debug, error, info, warn};
use polars::{
    error::PolarsResult,
    frame::DataFrame,
    io::SerWriter,
    prelude::{
        col, lit, when, CsvWriter, DataFrameJoinOps, IntoLazy, NamedFrom, QuoteStyle,
        SortMultipleOptions,
    },
    series::Series,
};

use crate::{
    common::{Detectable, TaxonomicRank},
    meta::{Meta, MetaColumn},
    options::CommonArgs,
    profile::{ProfileWrapper, ABUNDANCE_SUM_TOLERANCE},
    profile_handler::{ProfileHandler, ProfileHandlerError},
    taxonomy_inference::{
        describe_versions, resolve_meta, FileInference, InferenceOptions, RowTaxonomy,
    },
    tree_adjusted_benchmarks::get_adjusted_benchmarks,
    tree_handler::TreeHandler,
    utils::{
        add_string_columns, bray_curtis_similarity, f1_score, l2_similarity, pearson_correlation,
        precision, sensitivity, spearman_correlation,
    },
};

/// Run Benchpro using the provided CLI arguments.
///
/// # Arguments
///
/// * `args` - Parsed CLI arguments shared by subcommands
pub fn run(args: &CommonArgs) {
    if args.validate_meta {
        validate_meta_only(args);
        return;
    }
    let Some(outprefix) = args.outprefix.as_deref() else {
        error!("--outprefix is required unless --validate-meta is used");
        exit(2);
    };
    if args.meta.is_none() {
        error!("--meta is required");
        exit(2);
    }
    if let Err(message) = check_output_prefix(outprefix) {
        error!("{}", message);
        exit(2);
    }
    meta_based_workflow(args);
}

/// Checks that the folder for `outprefix` exists, so a run does not fail after computing.
///
/// # Arguments
///
/// * `outprefix` - Output prefix such as `results/run1`
///
/// # Returns
///
/// `Ok(())` when output files can be created next to the prefix.
///
/// # Errors
///
/// Returns a user-facing message when the output folder is missing.
///
/// # Examples
///
/// ```ignore
/// assert!(check_output_prefix("run1").is_ok());
/// assert!(check_output_prefix("missing_dir/run1").is_err());
/// ```
fn check_output_prefix(outprefix: &str) -> Result<(), String> {
    let folder = match Path::new(outprefix).parent() {
        Some(parent) if !parent.as_os_str().is_empty() => parent,
        _ => return Ok(()),
    };
    if folder.is_dir() {
        Ok(())
    } else {
        Err(format!(
            "Output folder '{}' does not exist. Create it with: mkdir -p {}",
            folder.display(),
            folder.display()
        ))
    }
}

fn validate_meta_only(args: &CommonArgs) {
    let meta_path = match &args.meta {
        Some(path) => path,
        None => {
            error!("--validate-meta requires --meta");
            exit(2);
        }
    };

    let mut meta = match Meta::from_path(meta_path) {
        Ok(meta) => meta,
        Err(err) => {
            error!("Meta validation failed: {}", err);
            exit(1);
        }
    };

    if let Err(err) = meta.validate_paths() {
        error!("Meta validation failed: {}", err);
        exit(1);
    }

    match resolve_meta(&mut meta, &inference_options(args)) {
        Ok(rows) => {
            for row in &rows {
                info!(
                    "Row '{}': Taxonomy {} ({})",
                    row.id, row.taxonomy, row.source
                );
            }
            log_taxonomy_warnings(&rows, &meta);
        }
        Err(details) => {
            error!("Meta validation failed: taxonomy problems in the meta file:");
            details.lines().for_each(|line| error!("  {}", line));
            exit(1);
        }
    }

    info!("Meta validation succeeded");
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
    let abundance_metrics = compute_abundance_metrics(&df)?;
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

    newdf.left_join(
        &abundance_metrics,
        ["Rank", "ID", "AllowAlternatives", "Adjusted"],
        ["Rank", "ID", "AllowAlternatives", "Adjusted"],
    )
}

#[derive(Debug)]
struct AbundanceMetrics {
    prediction: Vec<f64>,
    gold_std: Vec<f64>,
    prediction_tp: Vec<f64>,
    gold_std_tp: Vec<f64>,
}

/// Normalizes abundances to fractions; returns whether a percent-scale profile was rescaled.
fn normalize_profile_or_panic(
    profile: &mut ProfileWrapper,
    path: &Path,
    label: &str,
    ignore_error: bool,
) -> bool {
    let report = match profile.normalize_abundances_with_checks(ABUNDANCE_SUM_TOLERANCE) {
        Ok(report) => report,
        Err(e) => {
            if ignore_error {
                warn!(
                    "Ignoring abundance normalization error for {} profile {}: {}",
                    label,
                    path.display(),
                    e
                );
                return false;
            }
            panic!(
                "Abundance normalization failed for {} profile {}: {}",
                label,
                path.display(),
                e
            )
        }
    };

    if report.applied {
        debug!(
            "Rescaled {} profile {} from percentages: pre={:?} post={:?}",
            label,
            path.display(),
            report.pre_sums,
            report.post_sums
        );
    }

    // Percent-scale sums are expected before rescaling; only sums still above 1 are a problem.
    let mut over = report
        .post_sums
        .iter()
        .filter(|(_, sum)| **sum > 1.0 + ABUNDANCE_SUM_TOLERANCE)
        .map(|(rank, sum)| format!("{} {:.3}", rank, sum))
        .collect::<Vec<_>>();
    if !over.is_empty() {
        over.sort();
        warn!(
            "Abundances of {} profile {} sum to more than 1 at some ranks: {}",
            label,
            path.display(),
            over.join(", ")
        );
    }
    report.applied
}

fn compute_abundance_metrics(df: &DataFrame) -> PolarsResult<DataFrame> {
    let rank_col = df.column("Rank")?.str()?;
    let id_col = df.column("ID")?.str()?;
    let allow_alternatives_col = df.column("AllowAlternatives")?.bool()?;
    let adjusted_col = df.column("Adjusted")?.str()?;
    let type_col = df.column("Type")?.str()?;
    let pred_abundance_col = df.column("PredictionAbundance")?.f64()?;
    let gold_abundance_col = df.column("GoldStdAbundance")?.f64()?;

    let mut groups: HashMap<(String, String, bool, String), AbundanceMetrics> = HashMap::new();

    for row_index in 0..df.height() {
        let rank = rank_col.get(row_index).unwrap_or("");
        let id = id_col.get(row_index).unwrap_or("");
        let allow_alternatives = allow_alternatives_col.get(row_index).unwrap_or(false);
        let adjusted = adjusted_col.get(row_index).unwrap_or("");
        let bc_type = type_col.get(row_index).unwrap_or("");
        let pred_abundance = pred_abundance_col.get(row_index).unwrap_or(0.0);
        let gold_abundance = gold_abundance_col.get(row_index).unwrap_or(0.0);

        let key = (
            rank.to_string(),
            id.to_string(),
            allow_alternatives,
            adjusted.to_string(),
        );
        let entry = groups.entry(key).or_insert_with(|| AbundanceMetrics {
            prediction: Vec::new(),
            gold_std: Vec::new(),
            prediction_tp: Vec::new(),
            gold_std_tp: Vec::new(),
        });

        entry.prediction.push(pred_abundance);
        entry.gold_std.push(gold_abundance);
        if bc_type == "TP" {
            entry.prediction_tp.push(pred_abundance);
            entry.gold_std_tp.push(gold_abundance);
        }
    }

    let mut ranks = Vec::with_capacity(groups.len());
    let mut ids = Vec::with_capacity(groups.len());
    let mut allow_alternatives = Vec::with_capacity(groups.len());
    let mut adjusted = Vec::with_capacity(groups.len());
    let mut bray_curtis = Vec::with_capacity(groups.len());
    let mut l2_sim = Vec::with_capacity(groups.len());
    let mut pearson = Vec::with_capacity(groups.len());
    let mut spearman = Vec::with_capacity(groups.len());
    let mut l2_sim_tp = Vec::with_capacity(groups.len());
    let mut pearson_tp = Vec::with_capacity(groups.len());
    let mut spearman_tp = Vec::with_capacity(groups.len());

    for ((rank, id, allow, adjusted_val), values) in groups.into_iter() {
        let bc_sim = bray_curtis_similarity(&values.prediction, &values.gold_std);
        let l2_score = l2_similarity(&values.prediction, &values.gold_std);
        let pearson_score = pearson_correlation(&values.prediction, &values.gold_std);
        let spearman_score = spearman_correlation(&values.prediction, &values.gold_std);

        let l2_score_tp = if values.prediction_tp.is_empty() {
            f64::NAN
        } else {
            l2_similarity(&values.prediction_tp, &values.gold_std_tp)
        };
        let pearson_score_tp = if values.prediction_tp.len() < 2 {
            f64::NAN
        } else {
            pearson_correlation(&values.prediction_tp, &values.gold_std_tp)
        };
        let spearman_score_tp = if values.prediction_tp.len() < 2 {
            f64::NAN
        } else {
            spearman_correlation(&values.prediction_tp, &values.gold_std_tp)
        };
        debug!(
            "{} {}: l2 {} pc {} sc {} bc {} l2tp {} pctp {} sctp {}",
            rank,
            id,
            l2_score,
            pearson_score,
            spearman_score,
            bc_sim,
            l2_score_tp,
            pearson_score_tp,
            spearman_score_tp
        );

        ranks.push(rank);
        ids.push(id);
        allow_alternatives.push(allow);
        adjusted.push(adjusted_val);
        bray_curtis.push(bc_sim);
        l2_sim.push(l2_score);
        pearson.push(pearson_score);
        spearman.push(spearman_score);
        l2_sim_tp.push(l2_score_tp);
        pearson_tp.push(pearson_score_tp);
        spearman_tp.push(spearman_score_tp);
    }

    DataFrame::new(vec![
        Series::new("Rank".into(), ranks),
        Series::new("ID".into(), ids),
        Series::new("AllowAlternatives".into(), allow_alternatives),
        Series::new("Adjusted".into(), adjusted),
        Series::new("BrayCurtisSimilarity".into(), bray_curtis),
        Series::new("L2Similarity".into(), l2_sim),
        Series::new("PearsonCorrelation".into(), pearson),
        Series::new("SpearmanCorrelation".into(), spearman),
        Series::new("L2SimilarityTP".into(), l2_sim_tp),
        Series::new("PearsonCorrelationTP".into(), pearson_tp),
        Series::new("SpearmanCorrelationTP".into(), spearman_tp),
    ])
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
pub fn get_taxon_df(
    handler: &ProfileHandler,
    allow_alternatives: bool,
    ignore_abundance_error: bool,
) -> DataFrame {
    debug!("Number of profiles: {}", handler.prediction_map.len());
    debug!("Number of gold profiles: {}", handler.gold_std_map.len());

    let mut rescaled = HashSet::new();
    let mut dfs: Vec<DataFrame> = Vec::default();
    let _dfs_with_taxa: Vec<DataFrame> = Vec::default();

    for row in handler.meta.entries.iter() {
        let prediction = handler.prediction_map.get(row.profile.to_str().unwrap());
        let goldstd = handler.gold_std_map.get(row.goldstd.to_str().unwrap());

        let taxa = row.taxa_list.as_ref().map(|x| {
            handler
                .available_taxa
                .get(x.to_str().unwrap())
                .expect("Should be loaded...")
        });

        let not_loaded = match (prediction.is_some(), goldstd.is_some()) {
            (true, true) => None,
            (false, true) => Some(format!("prediction '{}'", row.profile.display())),
            (true, false) => Some(format!("gold standard '{}'", row.goldstd.display())),
            (false, false) => Some(format!(
                "prediction '{}' and gold standard '{}'",
                row.profile.display(),
                row.goldstd.display()
            )),
        };
        if let Some(what) = not_loaded {
            warn!(
                "Row '{}' not evaluated: {} could not be loaded (see the error above)",
                row.id, what
            );
        }

        if let (Some(prediction), Some(goldstd)) = (prediction, goldstd) {
            let mut prediction = prediction.clone();
            let mut goldstd = goldstd.clone();

            for (profile, path, label) in [
                (&mut prediction, &row.profile, "prediction"),
                (&mut goldstd, &row.goldstd, "gold standard"),
            ] {
                if normalize_profile_or_panic(profile, path, label, ignore_abundance_error) {
                    rescaled.insert(path.clone());
                }
            }

            debug!("Dataset: {} ... {:?}", row.id, row.goldstd);
            let df = prediction.binary_classification(&goldstd, allow_alternatives);
            if let Err(e) = &df {
                warn!("Row '{}' not evaluated: {}", row.id, e);
            }

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
    if !rescaled.is_empty() {
        info!(
            "Rescaled {} profile(s) from percentages to fractions",
            rescaled.len()
        );
    }
    if dfs.is_empty() {
        error!("No row could be evaluated; see the warnings above");
        exit(1);
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
/// * `args` - Parsed CLI arguments shared by subcommands
pub fn meta_based_workflow(args: &CommonArgs) {
    let path = Path::new(args.meta.as_ref().unwrap());
    let outprefix = args
        .outprefix
        .as_deref()
        .expect("--outprefix is required unless --validate-meta is used");

    let handler = match ProfileHandler::from_meta_with(path, &inference_options(args)) {
        Ok(handler) => handler,
        Err(ProfileHandlerError::TaxonomyError(details)) => {
            error!("Taxonomy problems in the meta file:");
            details.lines().for_each(|line| error!("  {}", line));
            if details.contains("cannot tell the taxonomy") {
                error!(
                    "Fix: add a Taxonomy column (GTDB or NCBI) for these rows, or use \
                     --taxonomy-index <DIR> to check names against GTDB/NCBI releases."
                );
            }
            exit(1);
        }
        Err(e) => {
            error!("Failed to load meta file '{}': {}", path.display(), e);
            exit(1);
        }
    };
    let mut warning_count = report_taxonomy(&handler);
    warning_count += check_adjusted_inputs(args, &handler.meta);

    let tree_handler = match TreeHandler::from_meta(&handler.meta) {
        Ok(th) => th,
        Err(e) => {
            error!("Failed to load GoldStdTree: {}", e);
            exit(1);
        }
    };

    ///////////////////////////////////////////////////////////////////////
    // Binary classification DF

    let complete_df = get_taxon_df(
        &handler,
        args.allow_alternatives,
        args.ignore_abundance_error,
    );

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
                    .left_join_to(
                        &dfa,
                        &[MetaColumn::Dataset, MetaColumn::Taxonomy, MetaColumn::Tool],
                        true,
                    )
                    .expect("Cannot join dfs?");

                debug!("Height DFA {}", dfa.height());

                let _ =
                    add_string_columns(&mut dfa, &[("Adjusted".to_string(), "True".to_string())]);

                debug!("Height DFA {}", dfa.height());
                let _ = add_string_columns(
                    &mut new_complete_df,
                    &[("ClosestNeighbor", ""), ("ClosestNeighborType", "")]
                        .map(|(a, b)| (a.to_string(), b.to_string())),
                );
                new_complete_df
                    .with_column(Series::new(
                        "ClosestNeighborDistance".into(),
                        std::iter::repeat(0f64)
                            .take(new_complete_df.height())
                            .collect_vec(),
                    ))
                    .unwrap();
                new_complete_df
                    .with_column(Series::new(
                        "ClosestNeighborAbundance".into(),
                        std::iter::repeat(0f64)
                            .take(new_complete_df.height())
                            .collect_vec(),
                    ))
                    .unwrap();

                let reorder = new_complete_df
                    .get_column_names()
                    .into_iter()
                    .map(|x| x.to_owned());
                let dfa = dfa.select(reorder).unwrap();

                new_complete_df = new_complete_df.vstack(&dfa).unwrap();
            }
            Err(e) => {
                warn!("Adjusted benchmarks failed: {:?}", e);
            }
        }
    }

    new_complete_df = new_complete_df
        .sort(
            [
                "ID",
                "Rank",
                "Name",
                "Type",
                "AllowAlternatives",
                "Adjusted",
            ],
            SortMultipleOptions::default(),
        )
        .expect("Cannot sort detailed output");

    // Specify the path where you want to save the CSV file

    let file_name_detailed = format!("{}_detailed.tsv", outprefix);
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

    debug!("new_complete_df:\n{:#?}", new_complete_df);

    let newdf =
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
    warning_count += check_results(&newdf);

    let file_name_bc = format!("{}.tsv", outprefix);
    let file_path: &Path = Path::new(&file_name_bc);
    let mut file = File::create(file_path).expect("Could not create file");

    // Write the DataFrame to the CSV file using CsvWriter
    let _ = CsvWriter::new(&mut file)
        .include_header(true)
        .with_separator(b'\t')
        .with_quote_style(QuoteStyle::Never)
        .finish(&mut newdf);

    debug!("NEWDF {}", newdf);
    let file_name_inputs = format!("{}_inputs.tsv", outprefix);
    if let Err(e) = write_inputs_report(&file_name_inputs, &handler, args) {
        warn!("Could not write {}: {}", file_name_inputs, e);
        warning_count += 1;
    }

    let rows = handler.meta.entries.len();
    let inferred = handler
        .taxonomy_rows
        .iter()
        .filter(|r| r.source == "inferred")
        .count();
    let adjusted = handler
        .meta
        .entries
        .iter()
        .filter(|e| e.goldstd_tree.is_some())
        .count();
    let evaluated = newdf
        .column("ID")
        .ok()
        .and_then(|c| c.str().ok())
        .map(|c| {
            c.into_iter()
                .flatten()
                .map(str::to_string)
                .collect::<HashSet<_>>()
        })
        .unwrap_or_default();
    let not_evaluated = handler
        .meta
        .entries
        .iter()
        .filter(|e| !evaluated.contains(&e.id))
        .map(|e| e.id.as_str())
        .collect::<Vec<_>>();
    info!(
        "Done: {} rows, taxonomy given for {}, inferred for {}",
        rows,
        rows - inferred,
        inferred
    );
    if !not_evaluated.is_empty() {
        warn!(
            "{} row(s) not evaluated (reasons in the warnings above): {}. For CAMI lineage \
             length errors, add ProfileColumns/GoldStdColumns or use --cami-ignore-lineage-error",
            not_evaluated.len(),
            preview(&not_evaluated)
        );
        warning_count += 1;
    }
    if args.adjusted {
        info!(
            "Adjusted: {} rows, skipped {} (no GoldStdTree)",
            adjusted,
            rows - adjusted
        );
    }
    info!("Warnings: {}", warning_count);
    info!("Summary:  {}", file_name_bc);
    info!("Per taxon: {}", file_name_detailed);
    info!(
        "Inputs (what was given, inferred and skipped): {}",
        file_name_inputs
    );
}

fn inference_options(args: &CommonArgs) -> InferenceOptions {
    InferenceOptions {
        index_dir: args.resolved_taxonomy_index(),
        allow_taxonomy_mismatch: args.allow_taxonomy_mismatch,
    }
}

/// Logs inferred taxonomies and grouped warnings; returns the number of warning lines.
fn report_taxonomy(handler: &ProfileHandler) -> usize {
    let mut inferred: Vec<(&str, &str)> = handler
        .taxonomy_rows
        .iter()
        .filter(|r| r.source == "inferred")
        .map(|r| (r.taxonomy.as_str(), r.id.as_str()))
        .collect();
    inferred.sort();
    for (taxonomy, group) in &inferred.iter().chunk_by(|(taxonomy, _)| *taxonomy) {
        let ids = group.map(|(_, id)| *id).collect::<Vec<_>>();
        info!(
            "Taxonomy inferred as {} for {} row(s): {}",
            taxonomy,
            ids.len(),
            preview(&ids)
        );
    }
    log_taxonomy_warnings(&handler.taxonomy_rows, &handler.meta)
}

/// Logs one line per warning kind (rows, tools, one example); details go to `_inputs.tsv`.
///
/// # Arguments
///
/// * `rows` - Taxonomy decisions per meta row
/// * `meta` - Parsed meta table, used to name the tools involved
///
/// # Returns
///
/// Number of warning lines logged.
fn log_taxonomy_warnings(rows: &[RowTaxonomy], meta: &Meta) -> usize {
    let tool_of = meta
        .entries
        .iter()
        .map(|e| (e.id.as_str(), e.tool.as_str()))
        .collect::<HashMap<_, _>>();
    let mut kinds: Vec<&str> = Vec::new();
    let mut by_kind: HashMap<&str, Vec<(&str, &str)>> = HashMap::new();
    for row in rows {
        for warning in &row.warnings {
            if !by_kind.contains_key(warning.kind) {
                kinds.push(warning.kind);
            }
            by_kind
                .entry(warning.kind)
                .or_default()
                .push((row.id.as_str(), warning.message.as_str()));
        }
    }
    for kind in &kinds {
        let items = &by_kind[kind];
        let tools = items
            .iter()
            .map(|(id, _)| tool_of.get(id).copied().unwrap_or("?"))
            .unique()
            .collect::<Vec<_>>();
        warn!(
            "{}: {} row(s), tools {}. Example: {}",
            kind,
            items.iter().map(|(id, _)| id).unique().count(),
            preview(&tools),
            items[0].1
        );
    }
    if !kinds.is_empty() {
        warn!("Per-row details: Warnings column of <outprefix>_inputs.tsv");
    }
    kinds.len()
}

/// Checks which rows `--adjusted` can run on; exits when none can.
///
/// # Arguments
///
/// * `args` - CLI arguments
/// * `meta` - Parsed meta table
///
/// # Returns
///
/// Number of warnings logged.
fn check_adjusted_inputs(args: &CommonArgs, meta: &Meta) -> usize {
    if !args.adjusted {
        return 0;
    }
    let without_tree = meta
        .entries
        .iter()
        .filter(|e| e.goldstd_tree.is_none())
        .map(|e| e.id.as_str())
        .collect::<Vec<_>>();
    if without_tree.len() == meta.entries.len() {
        error!("--adjusted needs a GoldStdTree for at least one row; none has one");
        exit(2);
    }
    let mut warnings = 0;
    if !without_tree.is_empty() {
        warn!(
            "--adjusted skips {} row(s) without GoldStdTree: {}",
            without_tree.len(),
            preview(&without_tree)
        );
        warnings += 1;
    }
    let lenient = meta
        .entries
        .iter()
        .filter(|e| e.goldstd_tree.is_some() && e.taxa_list.is_none())
        .map(|e| e.id.as_str())
        .collect::<Vec<_>>();
    if !lenient.is_empty() {
        warn!(
            "{} adjusted row(s) have no AvailableSpecies, so every FN counts as undetectable \
             and adjusted results are more lenient: {}",
            lenient.len(),
            preview(&lenient)
        );
        warnings += 1;
    }
    warnings
}

/// Formats up to five items, noting how many more there are.
fn preview(items: &[&str]) -> String {
    let shown = items.iter().take(5).copied().collect::<Vec<_>>().join(", ");
    match items.len() {
        n if n > 5 => format!("{} and {} more", shown, n - 5),
        _ => shown,
    }
}

/// Warns about results that usually mean a naming or taxonomy problem, not a bad tool.
///
/// # Arguments
///
/// * `summary` - Summary DataFrame
///
/// # Returns
///
/// Number of warnings logged.
fn check_results(summary: &DataFrame) -> usize {
    let zero_tp = summary
        .clone()
        .lazy()
        .filter(
            col("Rank")
                .eq(lit(TaxonomicRank::Species.to_string()))
                .and(col("AllowAlternatives").eq(lit(false)))
                .and(col("Adjusted").eq(lit("False")))
                .and(col("TP").eq(lit(0))),
        )
        .select([col("ID")])
        .collect();
    let Ok(zero_tp) = zero_tp else {
        return 0;
    };
    let ids = zero_tp
        .column("ID")
        .ok()
        .and_then(|c| c.str().ok())
        .map(|c| {
            c.into_iter()
                .flatten()
                .map(str::to_string)
                .collect::<Vec<_>>()
        })
        .unwrap_or_default();
    if ids.is_empty() {
        return 0;
    }
    warn!(
        "{} row(s) have 0 true positives at species level. This usually means the prediction \
         and gold standard name species differently (taxonomy or release mismatch): {}",
        ids.len(),
        preview(&ids.iter().map(String::as_str).collect::<Vec<_>>())
    );
    1
}

/// Writes one row per meta entry describing what was given, inferred, loaded and skipped.
///
/// # Arguments
///
/// * `path` - Output file path
/// * `handler` - Loaded profiles and taxonomy decisions
/// * `args` - CLI arguments
///
/// # Returns
///
/// `Ok(())` when the file is written.
///
/// # Errors
///
/// Returns a Polars error when the file cannot be created or written.
fn write_inputs_report(
    path: &str,
    handler: &ProfileHandler,
    args: &CommonArgs,
) -> PolarsResult<()> {
    let mut columns: Vec<(&str, Vec<String>)> = [
        "ID",
        "Sample",
        "Tool",
        "Profile",
        "GoldStd",
        "Taxonomy",
        "TaxonomySource",
        "ProfileTaxonomy",
        "ProfileConfidence",
        "ProfileRelease",
        "ProfileEvidence",
        "GoldStdTaxonomy",
        "GoldStdConfidence",
        "GoldStdRelease",
        "GoldStdEvidence",
        "ProfileLoaded",
        "ProfileSpecies",
        "GoldStdLoaded",
        "GoldStdSpecies",
        "GoldStdTree",
        "AvailableSpecies",
        "Adjusted",
        "Warnings",
    ]
    .into_iter()
    .map(|name| (name, Vec::new()))
    .collect();

    let taxonomy_by_id = handler
        .taxonomy_rows
        .iter()
        .map(|r| (r.id.as_str(), r))
        .collect::<HashMap<_, _>>();
    let describe = |f: Option<&FileInference>| -> [String; 4] {
        match f {
            None => ["unreadable".into(), "none".into(), "".into(), "".into()],
            Some(f) => [
                f.family.map_or("unknown".to_string(), |x| x.to_string()),
                f.confidence.to_string(),
                match f.release_candidates.len() {
                    0 => String::new(),
                    _ => describe_versions(&f.release_candidates),
                },
                f.evidence.clone(),
            ],
        }
    };
    let loaded = |map: &HashMap<String, ProfileWrapper>, p: &Path| -> [String; 2] {
        match map.get(p.to_str().unwrap_or_default()) {
            None => ["no".into(), "".into()],
            Some(profile) => [
                "yes".into(),
                profile
                    .taxa(&TaxonomicRank::Species)
                    .map_or(0, |t| t.len())
                    .to_string(),
            ],
        }
    };

    for entry in &handler.meta.entries {
        let row = taxonomy_by_id.get(entry.id.as_str());
        let [pt, pc, pr, pe] = describe(row.and_then(|r| r.profile.as_ref()));
        let [gt, gc, gr, ge] = describe(row.and_then(|r| r.goldstd.as_ref()));
        let [pl, ps] = loaded(&handler.prediction_map, &entry.profile);
        let [gl, gs] = loaded(&handler.gold_std_map, &entry.goldstd);
        let path_or_na = |p: &Option<std::path::PathBuf>| {
            p.as_ref()
                .map_or("NA".to_string(), |p| p.display().to_string())
        };
        let adjusted = match (args.adjusted, &entry.goldstd_tree) {
            (false, _) => "not requested",
            (true, None) => "skipped: no GoldStdTree",
            (true, Some(_)) if entry.taxa_list.is_none() => "run (lenient: no AvailableSpecies)",
            (true, Some(_)) => "run",
        };
        let values = [
            entry.id.clone(),
            entry.sample.clone(),
            entry.tool.clone(),
            entry.profile.display().to_string(),
            entry.goldstd.display().to_string(),
            entry.taxonomy.clone(),
            row.map_or(String::new(), |r| r.source.clone()),
            pt,
            pc,
            pr,
            pe,
            gt,
            gc,
            gr,
            ge,
            pl,
            ps,
            gl,
            gs,
            path_or_na(&entry.goldstd_tree),
            path_or_na(&entry.taxa_list),
            adjusted.to_string(),
            row.map_or(String::new(), |r| {
                r.warnings
                    .iter()
                    .map(|w| w.message.as_str())
                    .collect::<Vec<_>>()
                    .join(" | ")
            }),
        ];
        for ((_, column), value) in columns.iter_mut().zip(values) {
            column.push(value.replace(['\t', '\n'], " "));
        }
    }

    let mut df = DataFrame::new(
        columns
            .into_iter()
            .map(|(name, values)| Series::new(name.into(), values))
            .collect(),
    )?;
    let mut file = File::create(path)?;
    CsvWriter::new(&mut file)
        .include_header(true)
        .with_separator(b'\t')
        .with_quote_style(QuoteStyle::Never)
        .finish(&mut df)
}

#[cfg(test)]
mod tests {
    use polars::prelude::{DataFrame, NamedFrom, Series};

    use super::{add_binary_classification, check_output_prefix};

    fn approx_eq(left: f64, right: f64) -> bool {
        (left - right).abs() <= 1e-12
    }

    #[test]
    fn test_add_binary_classification_summary_metrics() {
        let df = DataFrame::new(vec![
            Series::new("Name".into(), &["s__a", "s__b", "s__c", "s__d"]),
            Series::new("Type".into(), &["TP", "TP", "FP", "FN"]),
            Series::new("PredictionAbundance".into(), &[0.5, 0.3, 0.2, 0.0]),
            Series::new("GoldStdAbundance".into(), &[0.6, 0.3, 0.0, 0.1]),
            Series::new("Rank".into(), &["species"; 4]),
            Series::new("ID".into(), &["S1"; 4]),
            Series::new("AllowAlternatives".into(), &[false; 4]),
            Series::new("Adjusted".into(), &["False"; 4]),
        ])
        .unwrap();

        let summary = add_binary_classification(df).unwrap();
        assert_eq!(summary.height(), 1);

        let get_u32 = |name: &str| summary.column(name).unwrap().u32().unwrap().get(0).unwrap();
        let get_f64 = |name: &str| summary.column(name).unwrap().f64().unwrap().get(0).unwrap();

        assert_eq!((get_u32("TP"), get_u32("FP"), get_u32("FN")), (2, 1, 1));
        assert!(approx_eq(get_f64("F1"), 4.0 / 6.0));
        assert!(approx_eq(get_f64("Sensitivity"), 2.0 / 3.0));
        assert!(approx_eq(get_f64("Precision"), 2.0 / 3.0));
        // |0.5-0.6| + 0 + 0.2 + 0.1 = 0.4 over a total of 2.0.
        assert!(approx_eq(get_f64("BrayCurtisSimilarity"), 0.8));
        // Two TP points always lie on a line.
        assert!(approx_eq(get_f64("PearsonCorrelationTP"), 1.0));
        assert!(approx_eq(get_f64("SpearmanCorrelationTP"), 1.0));
    }

    #[test]
    fn test_check_output_prefix() {
        assert!(check_output_prefix("run1").is_ok());
        assert!(check_output_prefix("src/run1").is_ok());

        let err = check_output_prefix("no_such_folder_xyz/run1").unwrap_err();
        assert!(err.contains("mkdir -p no_such_folder_xyz"), "{}", err);
    }
}
