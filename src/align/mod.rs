//! Read-alignment benchmark.
//!
//! Where [`crate::benchpro`] asks *which taxa are in a sample*, this module asks a different
//! question about a different input: given a read aligner's SAM/PAF output and per-read ground
//! truth, *did each read go to the right place, and is the alignment it reports any good?*
//!
//! It shares none of the taxonomic machinery (no `Taxonomy`, no `TaxonomicRank`) — only benchpro's
//! shape: a meta samplesheet in, aggregated TSVs out.
//!
//! The module scores files that already exist; it never runs an aligner and never measures time.
//! Orchestration (building indexes, timing runs, cache control) belongs to the harness that calls
//! it.
//!
//! The unit of work is a `(dataset, sample)` **group**, not a row: the recall denominator and the
//! base-level head-to-head both need every contender on a sample in memory at once.

pub mod cigar;
pub mod error;
pub mod mapq;
pub mod meta;
pub mod metrics;
pub mod report;
pub mod sam;
pub mod truth;

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use log::{debug, info, warn};

use crate::options::{AlignArgs, ScoringMode};

use error::{AlignError, AlignResult};
use meta::{AlignMeta, AlignRow};
use metrics::{ReadVerdict, ScoringContext};
use report::SampleResult;
use sam::ParsedAlignment;
use truth::Truth;

/// A per-read verdict tagged with the run it came from.
type TaggedVerdict = (String, String, String, ScoringMode, ReadVerdict);

/// Runs the alignment benchmark described by `args.meta`.
///
/// # Arguments
///
/// * `args` - Parsed `align` subcommand arguments
///
/// # Returns
///
/// `Ok(())` once every output table has been written.
///
/// # Errors
///
/// Returns an [`AlignError`] when the samplesheet is invalid, an input cannot be read or an output
/// cannot be written.
pub fn run_align(args: &AlignArgs) -> AlignResult<()> {
    let meta = AlignMeta::from_path(&args.meta, args)?;
    info!(
        "Loaded {} row(s) from '{}': {} dataset(s), {} tool(s)",
        meta.rows.len(),
        args.meta.display(),
        meta.datasets().len(),
        meta.tools().len()
    );

    if args.validate_meta {
        info!("Meta file is valid");
        return Ok(());
    }

    let mut results: Vec<SampleResult> = Vec::new();
    let mut per_read: Vec<TaggedVerdict> = Vec::new();

    for ((dataset, sample), rows) in meta.groups() {
        info!(
            "Scoring {}/{}: {} contender(s)",
            dataset,
            sample,
            rows.len()
        );
        let (group_results, group_reads) = score_group(dataset, sample, &rows, args)?;
        results.extend(group_results);
        per_read.extend(group_reads);
    }

    let (results, notes) = report::restrict_to_common_samples(results);
    let summaries = report::summarize(&results, &notes);
    log_summary(&summaries);

    let prefix = &args.outprefix;
    report::write(
        &mut report::summary_frame(&summaries)?,
        &output_path(prefix, "align_summary.tsv"),
    )?;
    report::write(
        &mut report::samples_frame(&results)?,
        &output_path(prefix, "align_samples.tsv"),
    )?;
    report::write(
        &mut report::mapq_frame(&summaries)?,
        &output_path(prefix, "align_mapq.tsv"),
    )?;
    if args.per_read {
        report::write(
            &mut report::reads_frame(&per_read)?,
            &output_path(prefix, "align_reads.tsv"),
        )?;
    }

    Ok(())
}

/// Scores every contender on one `(dataset, sample)`.
///
/// The group is the unit because the recall denominator is a property of the whole field, not of
/// one tool: it is the most any single contender got right (see [`mapq::mappable_base`]).
fn score_group(
    dataset: &str,
    sample: &str,
    rows: &[&AlignRow],
    args: &AlignArgs,
) -> AlignResult<(Vec<SampleResult>, Vec<TaggedVerdict>)> {
    // One truth and one contig map per path, however many contenders share them.
    let mut truths: HashMap<PathBuf, Truth> = HashMap::new();
    let mut maps: HashMap<PathBuf, HashMap<Box<str>, Box<str>>> = HashMap::new();
    let mut parsed: Vec<ParsedAlignment> = Vec::with_capacity(rows.len());

    for row in rows {
        if !truths.contains_key(&row.truth) {
            let truth = truth::load_truth(&row.truth)?;
            debug!("{}: {} truth reads", row.truth.display(), truth.len());
            truths.insert(row.truth.clone(), truth);
        }
        if let Some(path) = &row.contig2genome {
            if !maps.contains_key(path) {
                maps.insert(path.clone(), truth::load_contig2genome(path)?);
            }
        }

        let keep_seq = !args.no_replay && row.reference.is_some() && row.format.has_alignments();
        let alignment = sam::parse_alignment(
            &row.alignment,
            row.format,
            keep_seq,
            args.verify_sample,
            args.threads,
            args.seed,
        )?;
        info!(
            "  {}: {} primary alignment(s) from {}",
            row.tool,
            alignment.records.len(),
            row.alignment.display()
        );
        parsed.push(alignment);
    }

    // Score each contender, then fill in the recall denominator once the whole field is known.
    let mut scored = Vec::with_capacity(rows.len());
    let mut reads = Vec::new();

    for (row, alignment) in rows.iter().zip(parsed.iter()) {
        let truth = &truths[&row.truth];
        let context = ScoringContext {
            scoring: row.scoring,
            sep: &row.sep,
            contig2genome: row.contig2genome.as_ref().and_then(|p| maps.get(p)),
            tolerance: args.tolerance,
        };

        // Zero overlap is almost always a read-id convention mismatch rather than a genuinely
        // useless aligner, and it would otherwise be reported as a flawless run of all zeroes.
        if !alignment.records.is_empty()
            && !truth.is_empty()
            && !alignment.records.keys().any(|key| truth.contains_key(key))
        {
            return Err(AlignError::NoSharedReads {
                alignment: row.alignment.clone(),
                truth: row.truth.clone(),
                aligned: alignment.records.len(),
                truth_reads: truth.len(),
            });
        }

        let score = metrics::score(&alignment.records, truth, &context);
        let mapq_counts = mapq::counts(&alignment.records, truth, &context);

        if args.per_read {
            for verdict in metrics::read_verdicts(&alignment.records, truth, &context) {
                reads.push((
                    dataset.to_string(),
                    sample.to_string(),
                    row.tool.clone(),
                    row.scoring,
                    verdict,
                ));
            }
        }

        scored.push(SampleResult {
            dataset: dataset.to_string(),
            sample: sample.to_string(),
            tool: row.tool.clone(),
            scoring: row.scoring,
            score,
            counters: alignment.counters.clone(),
            mapq_counts,
            mappable: 0,
        });
    }

    let mappable = mapq::mappable_base(scored.iter().map(|s| s.mapq_counts.as_slice()));
    if mappable == 0 {
        warn!(
            "{}/{}: no contender placed a single read correctly, so recall has no denominator",
            dataset, sample
        );
    }
    for result in &mut scored {
        result.mappable = mappable;
    }

    Ok((scored, reads))
}

/// Appends a suffix to the output prefix.
fn output_path(prefix: &str, suffix: &str) -> PathBuf {
    Path::new(&format!("{prefix}.{suffix}")).to_path_buf()
}

/// Logs the metric table the Python tool prints: metric label down the left, one column per tool.
fn log_summary(summaries: &[report::ToolSummary]) {
    type Render = Box<dyn Fn(&report::ToolSummary) -> String>;

    for dataset in unique_datasets(summaries) {
        let group: Vec<&report::ToolSummary> =
            summaries.iter().filter(|s| s.dataset == dataset).collect();

        let mut header = format!("{:<40}", format!("[{dataset}]"));
        for summary in &group {
            header.push_str(&format!("{:>18}", summary.tool));
        }
        info!("{}", header);

        let rows: Vec<(&str, Render)> = vec![
            (
                "samples",
                Box::new(|s: &report::ToolSummary| s.samples.to_string()),
            ),
            (
                "reads in truth",
                Box::new(|s: &report::ToolSummary| s.score.total.to_string()),
            ),
            (
                "aligned",
                Box::new(|s: &report::ToolSummary| s.score.aligned.to_string()),
            ),
            (
                "aligned (% of truth)",
                Box::new(|s: &report::ToolSummary| format!("{:.2}%", s.score.align_pct())),
            ),
            (
                "correct (% of aligned)",
                Box::new(|s: &report::ToolSummary| format!("{:.2}%", s.score.correct_pct())),
            ),
            (
                "recall (% of mappable)",
                Box::new(|s: &report::ToolSummary| format!("{:.2}%", s.recall_pct())),
            ),
            (
                "best MAPQ cutoff (F1)",
                Box::new(|s: &report::ToolSummary| {
                    s.best_f1()
                        .map(|p| format!("{} ({:.1})", p.mapq, p.f1()))
                        .unwrap_or_else(|| "n/a".to_string())
                }),
            ),
        ];

        for (label, render) in rows {
            let mut line = format!("{label:<40}");
            for summary in &group {
                line.push_str(&format!("{:>18}", render(summary)));
            }
            info!("{}", line);
        }

        if let Some(note) = group.iter().find_map(|s| s.note.as_ref()) {
            info!("note: {}", note);
        }
    }
}

/// Dataset labels in first-seen order.
fn unique_datasets(summaries: &[report::ToolSummary]) -> Vec<String> {
    let mut seen = std::collections::HashSet::new();
    summaries
        .iter()
        .filter(|s| seen.insert(s.dataset.clone()))
        .map(|s| s.dataset.clone())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn output_paths_hang_off_the_prefix() {
        assert_eq!(
            output_path("results/run1", "align_summary.tsv"),
            PathBuf::from("results/run1.align_summary.tsv")
        );
    }
}
