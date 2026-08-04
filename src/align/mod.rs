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

pub mod base;
pub mod cigar;
pub mod clip;
pub mod error;
pub mod mapq;
pub mod meta;
pub mod metrics;
pub mod reference;
pub mod report;
pub mod sam;
pub mod truth;

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use log::{debug, info, warn};

use crate::options::AlignArgs;

use error::{AlignError, AlignResult};
use meta::{AlignMeta, AlignRow};
use metrics::ScoringContext;
use report::{SampleResult, TaggedVerdict};
use sam::ParsedAlignment;
use truth::Truth;

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

    // Which samples survive is decided from the samplesheet, before any scoring: a sample that
    // will be dropped from the aggregates is not worth scoring, and settling it up front lets the
    // per-read table be written group by group instead of accumulated.
    let (retained, notes) = report::common_samples(&meta.rows);

    let mut results: Vec<SampleResult> = Vec::new();
    let mut reads_writer = args
        .per_read
        .then(|| report::ReadsWriter::create(&output_path(&args.outprefix, "align_reads.tsv")))
        .transpose()?;

    for ((dataset, sample), rows) in meta.groups() {
        if !retained.contains(&(dataset.to_string(), sample.to_string())) {
            debug!(
                "Skipping {}/{}: not every contender produced it",
                dataset, sample
            );
            continue;
        }
        info!(
            "Scoring {}/{}: {} contender(s)",
            dataset,
            sample,
            rows.len()
        );
        let (group_results, group_reads) = score_group(dataset, sample, &rows, args)?;
        results.extend(group_results);
        if let Some(writer) = &mut reads_writer {
            writer.append(&group_reads)?;
        }
    }

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
    if let Some(writer) = reads_writer {
        writer.finish();
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
        // The contig map is loaded first: a gold-standard SAM truth needs it to label each read's
        // genome, exactly as build_truth.py does when it turns one into a truth TSV.
        if let Some(path) = &row.contig2genome {
            if !maps.contains_key(path) {
                maps.insert(path.clone(), truth::load_contig2genome(path)?);
            }
        }
        if !truths.contains_key(&row.truth) {
            let truth = truth::load_truth_any(
                &row.truth,
                row.contig2genome.as_ref().and_then(|p| maps.get(p)),
                args.threads,
            )?;
            debug!("{}: {} truth reads", row.truth.display(), truth.len());
            truths.insert(row.truth.clone(), truth);
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

    replay_against_references(rows, &mut parsed, args)?;

    let contig_lengths = if args.clip_geometry {
        contig_lengths_for(rows, &parsed)?
    } else {
        Vec::new()
    };

    // Score each contender, then fill in the recall denominator once the whole field is known.
    let mut scored = Vec::with_capacity(rows.len());
    let mut reads = Vec::new();

    for (index, (row, alignment)) in rows.iter().zip(parsed.iter()).enumerate() {
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

        // Base-level metrics need a CIGAR and a SEQ, which PAF does not carry. Computing them
        // anyway would fill the row with zeroes that read as measurements: 0% coverage, 0
        // unclipped, 0 replayed.
        let base = row
            .format
            .has_alignments()
            .then(|| {
                base::summarize(
                    &alignment.records,
                    &alignment.counters,
                    Some(truth),
                    &context,
                )
            })
            .flatten();

        // The head-to-head needs the peer's records, which is why the group -- not the row -- is
        // the unit of work.
        let h2h = row.peer.as_ref().and_then(|peer_name| {
            rows.iter()
                .position(|candidate| &candidate.tool == peer_name)
                .filter(|position| *position != index)
                .map(|position| {
                    base::head_to_head(
                        &alignment.records,
                        &parsed[position].records,
                        peer_name,
                        args.tolerance,
                    )
                })
        });

        let clip = args
            .clip_geometry
            .then(|| {
                contig_lengths
                    .get(index)
                    .and_then(|lengths| clip::summarize(&alignment.records, lengths))
            })
            .flatten();

        scored.push(SampleResult {
            dataset: dataset.to_string(),
            sample: sample.to_string(),
            tool: row.tool.clone(),
            scoring: row.scoring,
            score,
            counters: alignment.counters.clone(),
            mapq_counts,
            mappable: 0,
            base,
            h2h,
            clip,
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

/// Replays each contender's retained alignments against its reference.
///
/// Contenders that name the same reference share one open handle and one filtered index: on the
/// marker reference the index alone is the expensive part, and every tool on a sample is aligned
/// against the same FASTA.
fn replay_against_references(
    rows: &[&AlignRow],
    parsed: &mut [ParsedAlignment],
    args: &AlignArgs,
) -> AlignResult<()> {
    if args.no_replay {
        return Ok(());
    }

    // Group the contenders by the reference they were aligned against.
    let mut by_reference: HashMap<&Path, Vec<usize>> = HashMap::new();
    for (index, row) in rows.iter().enumerate() {
        if let Some(reference) = &row.reference {
            if row.format.has_alignments() {
                by_reference
                    .entry(reference.as_path())
                    .or_default()
                    .push(index);
            }
        }
    }

    for (path, indices) in by_reference {
        // Only the contigs the alignments actually name are loaded: a marker reference has
        // millions of sequences and a scoring run touches a small fraction of them.
        let wanted: HashSet<Box<str>> = indices
            .iter()
            .flat_map(|i| parsed[*i].records.values())
            .filter(|record| record.seq.is_some())
            .map(|record| record.target.clone())
            .collect();
        if wanted.is_empty() {
            continue;
        }

        let mut reference = reference::Reference::open(path, Some(&wanted))?;
        for index in indices {
            let alignment = &mut parsed[index];
            let verified = reference::verify(
                &mut alignment.records,
                &mut reference,
                &mut alignment.counters,
            )?;
            info!(
                "  {}: {} alignment(s) replayed against {}",
                rows[index].tool,
                verified,
                path.display()
            );
        }
    }

    Ok(())
}

/// Contig lengths per contender, for the clip geometry.
///
/// Prefers the reference's `.fai`, which is authoritative; falls back to the SAM's own `@SQ`
/// header, which every well-formed SAM carries. A contig in neither is reported as unknown rather
/// than guessed at.
///
/// Contenders sharing a reference resolve it once, against the union of the contigs they actually
/// name — the marker reference has 14.5 M sequences, and loading all of them (once per contender)
/// would cost gigabytes for a run that touches a small fraction.
fn contig_lengths_for(
    rows: &[&AlignRow],
    parsed: &[ParsedAlignment],
) -> AlignResult<Vec<HashMap<Box<str>, u64>>> {
    let mut by_reference: HashMap<&Path, Vec<usize>> = HashMap::new();
    for (index, row) in rows.iter().enumerate() {
        if let Some(reference) = &row.reference {
            by_reference
                .entry(reference.as_path())
                .or_default()
                .push(index);
        }
    }

    let mut lengths: Vec<HashMap<Box<str>, u64>> = vec![HashMap::new(); rows.len()];

    for (path, indices) in by_reference {
        let wanted: HashSet<Box<str>> = indices
            .iter()
            .flat_map(|i| parsed[*i].records.values())
            .map(|record| record.target.clone())
            .collect();
        if wanted.is_empty() {
            continue;
        }
        let reference = reference::Reference::open(path, Some(&wanted))?;
        if reference.is_empty() {
            continue;
        }
        let shared = reference.lengths();
        for index in indices {
            lengths[index] = shared.clone();
        }
    }

    // Anything still empty had no reference, or one that named none of its contigs; the SAM's own
    // header is the fallback.
    for (index, row) in rows.iter().enumerate() {
        if lengths[index].is_empty() && row.format.has_alignments() {
            lengths[index] = sam::header_lengths(&row.alignment)?;
        }
    }

    Ok(lengths)
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
                "recall (% of truth)",
                Box::new(|s: &report::ToolSummary| format!("{:.2}%", s.recall_pct())),
            ),
            (
                "best MAPQ cutoff (of mappable)",
                Box::new(|s: &report::ToolSummary| {
                    s.best_f1_mappable()
                        .map(|p| format!("{} (F1 {:.1})", p.mapq, p.f1_mappable()))
                        .unwrap_or_else(|| "n/a".to_string())
                }),
            ),
            (
                "best MAPQ cutoff (of truth)",
                Box::new(|s: &report::ToolSummary| {
                    s.best_f1_total()
                        .map(|p| format!("{} (F1 {:.1})", p.mapq, p.f1_total()))
                        .unwrap_or_else(|| "n/a".to_string())
                }),
            ),
        ];

        // A metric no contender defines is left out entirely rather than printed as a row of
        // "n/a" -- a PAF-only field has no base-level section at all.
        let base_rows: Vec<(&str, Render)> = vec![
            (
                "alignments emitted",
                Box::new(|s: &report::ToolSummary| base_count(s, |b| b.records)),
            ),
            (
                "of those, replayed vs reference",
                Box::new(|s: &report::ToolSummary| base_count(s, |b| b.verified)),
            ),
            (
                "identity vs reference (replayed)",
                Box::new(|s: &report::ToolSummary| base_pct(s, |b| b.identity_v)),
            ),
            (
                "replayed identity >= 95%",
                Box::new(|s: &report::ToolSummary| base_pct(s, |b| b.identity_high_v)),
            ),
            (
                "NM tag matches the reference",
                Box::new(|s: &report::ToolSummary| base_pct(s, |b| b.nm_agree)),
            ),
            (
                "mean identity (as reported)",
                Box::new(|s: &report::ToolSummary| base_pct(s, |b| b.identity)),
            ),
            (
                "mean identity, right target",
                Box::new(|s: &report::ToolSummary| base_pct(s, |b| b.identity_correct)),
            ),
            (
                "mean identity, wrong target",
                Box::new(|s: &report::ToolSummary| base_pct(s, |b| b.identity_wrong)),
            ),
            (
                "query coverage",
                Box::new(|s: &report::ToolSummary| base_pct(s, |b| Some(b.coverage))),
            ),
            (
                "unclipped alignments",
                Box::new(|s: &report::ToolSummary| base_pct(s, |b| Some(b.full_length))),
            ),
            (
                "alignments with an indel",
                Box::new(|s: &report::ToolSummary| base_pct(s, |b| Some(b.indel))),
            ),
            (
                "flagged proper pair",
                Box::new(|s: &report::ToolSummary| base_pct(s, |b| Some(b.proper_pair))),
            ),
            (
                "malformed (CIGAR != SEQ len)",
                Box::new(|s: &report::ToolSummary| base_pct(s, |b| Some(b.malformed))),
            ),
            (
                "target not in the reference",
                Box::new(|s: &report::ToolSummary| base_pct(s, |b| b.unknown_ref)),
            ),
            (
                "missing NM tag",
                Box::new(|s: &report::ToolSummary| base_pct(s, |b| Some(b.no_nm))),
            ),
        ];

        let h2h_rows: Vec<(&str, Render)> = vec![
            (
                "peer",
                Box::new(|s: &report::ToolSummary| {
                    s.h2h.as_ref().map(|h| h.peer.clone()).unwrap_or_default()
                }),
            ),
            (
                "reads aligned by both",
                Box::new(|s: &report::ToolSummary| h2h_count(s, |h| h.common)),
            ),
            (
                "reads only this tool aligned",
                Box::new(|s: &report::ToolSummary| h2h_count(s, |h| h.only_self)),
            ),
            (
                "reads only the peer aligned",
                Box::new(|s: &report::ToolSummary| h2h_count(s, |h| h.only_peer)),
            ),
            (
                "same locus as peer",
                Box::new(|s: &report::ToolSummary| h2h_pct(s, |h| h.same_locus)),
            ),
            (
                "lower edit distance than peer",
                Box::new(|s: &report::ToolSummary| h2h_pct(s, |h| h.better)),
            ),
            (
                "equal edit distance",
                Box::new(|s: &report::ToolSummary| h2h_pct(s, |h| h.equal)),
            ),
            (
                "higher edit distance",
                Box::new(|s: &report::ToolSummary| h2h_pct(s, |h| h.worse)),
            ),
            (
                "mean NM delta vs peer (<0 = better)",
                Box::new(|s: &report::ToolSummary| {
                    s.h2h
                        .as_ref()
                        .and_then(|h| h.nm_delta)
                        .map(|v| format!("{v:.3}"))
                        .unwrap_or_default()
                }),
            ),
        ];

        let print = |rows: Vec<(&str, Render)>| {
            for (label, render) in rows {
                let cells: Vec<String> = group.iter().map(|s| render(s)).collect();
                if cells.iter().all(String::is_empty) {
                    continue;
                }
                let mut line = format!("{label:<40}");
                for cell in cells {
                    line.push_str(&format!("{cell:>18}"));
                }
                info!("{}", line);
            }
        };

        print(rows);
        if group.iter().any(|s| s.base.is_some()) {
            print(base_rows);
        }
        if group.iter().any(|s| s.h2h.is_some()) {
            print(h2h_rows);
        }

        if let Some(note) = group.iter().find_map(|s| s.note.as_ref()) {
            info!("note: {}", note);
        }
    }
}

/// A base-level count for the log table, empty when the tool has no base-level metrics.
fn base_count(summary: &report::ToolSummary, f: fn(&base::BaseMetrics) -> u64) -> String {
    summary
        .base
        .as_ref()
        .map(|b| f(b).to_string())
        .unwrap_or_default()
}

/// A base-level percentage for the log table, empty when it is not defined for this tool.
fn base_pct(summary: &report::ToolSummary, f: fn(&base::BaseMetrics) -> Option<f64>) -> String {
    summary
        .base
        .as_ref()
        .and_then(f)
        .map(|v| format!("{v:.2}%"))
        .unwrap_or_default()
}

/// A head-to-head count for the log table.
fn h2h_count(summary: &report::ToolSummary, f: fn(&base::HeadToHead) -> u64) -> String {
    summary
        .h2h
        .as_ref()
        .map(|h| f(h).to_string())
        .unwrap_or_default()
}

/// A head-to-head percentage for the log table.
fn h2h_pct(summary: &report::ToolSummary, f: fn(&base::HeadToHead) -> Option<f64>) -> String {
    summary
        .h2h
        .as_ref()
        .and_then(f)
        .map(|v| format!("{v:.2}%"))
        .unwrap_or_default()
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
