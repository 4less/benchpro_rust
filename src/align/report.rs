//! Result aggregation and the output tables.
//!
//! Aggregation rules here are decisions, not defaults, and two of them are load-bearing:
//!
//! * **Pooled rates, never means of rates.** Numerators and denominators are summed across samples
//!   and divided once. A mean of per-sample percentages weights a 1000-read sample like a
//!   million-read one.
//! * **Restrict to common samples.** Aggregates are sums, and a sum over ten samples is not
//!   comparable with a sum over five. Left alone that publishes nonsense that looks fine, so the
//!   sample set is intersected across contenders and the note goes into the output.

use std::collections::{HashMap, HashSet};
use std::path::Path;

use log::{info, warn};
use polars::prelude::{DataFrame, NamedFrom, Series};

use crate::options::ScoringMode;
use crate::utils::write_df;

use super::base::{BaseMetrics, HeadToHead};
use super::clip::ClipGeometry;
use super::error::{AlignError, AlignResult};
use super::mapq::{self, MapqCount};
use super::metrics::{pct, MappingScore, ReadVerdict};
use super::sam::ParseCounters;

/// Everything scored for one contender on one sample.
#[derive(Debug, Clone)]
pub struct SampleResult {
    /// Dataset label.
    pub dataset: String,
    /// Sample id.
    pub sample: String,
    /// Contender name.
    pub tool: String,
    /// Which correctness definition was used.
    pub scoring: ScoringMode,
    /// Mapping-level score.
    pub score: MappingScore,
    /// What the parse saw and skipped.
    pub counters: ParseCounters,
    /// Cumulative MAPQ counts, ascending in cutoff.
    pub mapq_counts: Vec<MapqCount>,
    /// Denominator of the MAPQ curve's recall axis, shared by every contender on this sample.
    /// Not used by [`Self::recall_pct`].
    pub mappable: u64,
    /// Base-level metrics, absent for a contender that emitted no alignment.
    pub base: Option<BaseMetrics>,
    /// Base-level comparison against the row's peer, when one was named.
    pub h2h: Option<HeadToHead>,
    /// Clip geometry, when `--clip-geometry` asked for it.
    pub clip: Option<ClipGeometry>,
}

impl SampleResult {
    /// Recall: correct placements as a share of *every* read in the truth.
    ///
    /// The companion of `correct_pct`, and both are reported because neither alone is honest.
    /// `correct_pct` divides by the reads the tool placed, which is precision — it rewards a tool
    /// for mapping less, since one that places its most confident 1% and gets them all right
    /// scores 100%. Recall uses a denominator the tool cannot influence, so it is comparable
    /// however permissive the tool is.
    ///
    /// Note this is *not* normalised by [`Self::mappable`]: that base exists for the MAPQ curve's
    /// axis, where an all-reads denominator would flatten every contender (see
    /// [`mapq::mappable_base`]).
    ///
    /// # Returns
    ///
    /// `100 * correct / total`, or 0 for an empty truth.
    pub fn recall_pct(&self) -> f64 {
        pct(self.score.correct, self.score.total)
    }
}

/// A contender's results pooled over every sample of one dataset.
#[derive(Debug, Clone)]
pub struct ToolSummary {
    /// Dataset label.
    pub dataset: String,
    /// Contender name.
    pub tool: String,
    /// Which correctness definition was used.
    pub scoring: ScoringMode,
    /// Samples that contributed.
    pub samples: usize,
    /// Pooled mapping score.
    pub score: MappingScore,
    /// Pooled parse counters.
    pub counters: ParseCounters,
    /// Pooled MAPQ counts.
    pub mapq_counts: Vec<MapqCount>,
    /// Pooled denominator of the MAPQ curve's recall axis.
    pub mappable: u64,
    /// Standard deviation of `correct_pct` across samples, when there is more than one.
    pub correct_pct_sd: Option<f64>,
    /// Base-level metrics averaged over the samples that have them.
    pub base: Option<BaseMetrics>,
    /// Head-to-head averaged over the samples that have one.
    pub h2h: Option<HeadToHead>,
    /// Clip geometry pooled over the samples that have it.
    pub clip: Option<ClipGeometry>,
    /// Note explaining any restriction applied to the sample set.
    pub note: Option<String>,
}

impl ToolSummary {
    /// Recall: correct placements as a share of every read in the truth.
    ///
    /// See [`SampleResult::recall_pct`] for why this denominator, and not the mappable base.
    ///
    /// # Returns
    ///
    /// `100 * correct / total`.
    pub fn recall_pct(&self) -> f64 {
        pct(self.score.correct, self.score.total)
    }

    /// The MAPQ cutoff maximising F1, with its recall and precision.
    ///
    /// # Returns
    ///
    /// The best point of the pooled curve, or `None` when there is no curve.
    pub fn best_f1(&self) -> Option<mapq::MapqPoint> {
        mapq::best_f1(&mapq::curve(&self.mapq_counts, self.mappable))
    }
}

/// Drops samples that not every contender produced, per dataset.
///
/// Aggregates are sums, and a sum over ten samples is not comparable with a sum over five. This can
/// only ever make a comparison narrower, never wrong.
///
/// # Arguments
///
/// * `results` - Every scored `(dataset, sample, tool)` result
///
/// # Returns
///
/// The retained results, and one note per dataset that lost samples.
pub fn restrict_to_common_samples(
    results: Vec<SampleResult>,
) -> (Vec<SampleResult>, HashMap<String, String>) {
    let mut per_dataset: HashMap<&str, HashMap<&str, HashSet<&str>>> = HashMap::new();
    for result in &results {
        per_dataset
            .entry(&result.dataset)
            .or_default()
            .entry(&result.tool)
            .or_default()
            .insert(&result.sample);
    }

    let mut keep: HashMap<String, HashSet<String>> = HashMap::new();
    let mut notes: HashMap<String, String> = HashMap::new();

    for (dataset, per_tool) in &per_dataset {
        let mut sets = per_tool.values();
        let Some(first) = sets.next() else { continue };
        let common: HashSet<&str> = sets.fold(first.clone(), |acc, set| &acc & set);
        let union: HashSet<&str> = per_tool.values().flatten().copied().collect();

        if per_tool.len() > 1 && common.len() < union.len() {
            let mut dropped: Vec<&str> = union.difference(&common).copied().collect();
            dropped.sort_unstable();
            let note = format!(
                "restricted to the {} sample(s) every contender produced; dropped {}",
                common.len(),
                dropped.join(", ")
            );
            warn!(
                "{}: sample(s) {} are missing for at least one contender — aggregating over the {} \
                 sample(s) every contender ran, so the totals stay comparable",
                dataset,
                dropped.join(", "),
                common.len()
            );
            notes.insert((*dataset).to_string(), note);
        }
        keep.insert(
            (*dataset).to_string(),
            common.into_iter().map(str::to_string).collect(),
        );
    }

    let retained = results
        .into_iter()
        .filter(|r| {
            keep.get(&r.dataset)
                .is_none_or(|samples| samples.contains(&r.sample))
        })
        .collect();

    (retained, notes)
}

/// Pools per-sample results into one row per `(dataset, tool)`.
///
/// # Arguments
///
/// * `results` - Per-sample results, already restricted to the common sample set
/// * `notes` - Per-dataset restriction notes from [`restrict_to_common_samples`]
///
/// # Returns
///
/// One summary per contender per dataset, in first-seen order.
pub fn summarize(results: &[SampleResult], notes: &HashMap<String, String>) -> Vec<ToolSummary> {
    let mut order: Vec<(String, String)> = Vec::new();
    let mut grouped: HashMap<(String, String), Vec<&SampleResult>> = HashMap::new();

    for result in results {
        let key = (result.dataset.clone(), result.tool.clone());
        if !grouped.contains_key(&key) {
            order.push(key.clone());
        }
        grouped.entry(key).or_default().push(result);
    }

    order
        .into_iter()
        .map(|key| {
            let group = &grouped[&key];
            let mut summary = ToolSummary {
                dataset: key.0.clone(),
                tool: key.1.clone(),
                scoring: group[0].scoring,
                samples: group.len(),
                score: MappingScore::default(),
                counters: ParseCounters::default(),
                mapq_counts: Vec::new(),
                mappable: 0,
                correct_pct_sd: None,
                base: pool_base(group),
                h2h: pool_h2h(group),
                clip: pool_clip(group),
                note: notes.get(&key.0).cloned(),
            };

            for result in group {
                summary.score.add(&result.score);
                summary.counters.add(&result.counters);
                mapq::add_counts(&mut summary.mapq_counts, &result.mapq_counts);
                summary.mappable += result.mappable;
            }

            if group.len() > 1 {
                let values: Vec<f64> = group.iter().map(|r| r.score.correct_pct()).collect();
                summary.correct_pct_sd = Some(sample_sd(&values));
            }

            summary
        })
        .collect()
}

/// Pools base-level metrics across the samples that have them.
///
/// Counts are summed and shares are weighted by each sample's alignment count — a sample with a
/// thousand alignments must not weigh the same as one with a million. A share that no sample
/// defines stays `None`: a tool that emits no `NM` has no reported identity, and rendering that as
/// zero would read as "every alignment is a complete mismatch".
fn pool_base(group: &[&SampleResult]) -> Option<BaseMetrics> {
    let present: Vec<(&BaseMetrics, u64)> = group
        .iter()
        .filter_map(|r| r.base.as_ref().map(|b| (b, b.records)))
        .collect();
    if present.is_empty() {
        return None;
    }

    let total: u64 = present.iter().map(|(_, n)| *n).sum();
    let verified: u64 = present.iter().map(|(b, _)| b.verified).sum();
    // A share over the replayed subset must be weighted by that subset, not by all alignments.
    let by_records = |f: fn(&BaseMetrics) -> f64| -> f64 {
        weighted(present.iter().map(|(b, n)| (f(b), *n))).unwrap_or(0.0)
    };
    let by_records_opt = |f: fn(&BaseMetrics) -> Option<f64>| -> Option<f64> {
        weighted(present.iter().filter_map(|(b, n)| f(b).map(|v| (v, *n))))
    };
    let by_verified = |f: fn(&BaseMetrics) -> Option<f64>| -> Option<f64> {
        weighted(
            present
                .iter()
                .filter_map(|(b, _)| f(b).map(|v| (v, b.verified))),
        )
    };

    Some(BaseMetrics {
        records: total,
        verified,
        identity: by_records_opt(|b| b.identity),
        identity_high: by_records_opt(|b| b.identity_high),
        identity_v: by_verified(|b| b.identity_v),
        identity_high_v: by_verified(|b| b.identity_high_v),
        nm_agree: by_verified(|b| b.nm_agree),
        identity_correct: by_records_opt(|b| b.identity_correct),
        identity_wrong: by_records_opt(|b| b.identity_wrong),
        coverage: by_records(|b| b.coverage),
        full_length: by_records(|b| b.full_length),
        indel: by_records(|b| b.indel),
        malformed: by_records(|b| b.malformed),
        unknown_ref: by_records_opt(|b| b.unknown_ref),
        no_nm: by_records(|b| b.no_nm),
        proper_pair: by_records(|b| b.proper_pair),
    })
}

/// Pools head-to-head results across the samples that have one.
fn pool_h2h(group: &[&SampleResult]) -> Option<HeadToHead> {
    let present: Vec<&HeadToHead> = group.iter().filter_map(|r| r.h2h.as_ref()).collect();
    let first = present.first()?;

    let common: u64 = present.iter().map(|h| h.common).sum();
    let weighted_by_common = |f: fn(&HeadToHead) -> Option<f64>| -> Option<f64> {
        weighted(present.iter().filter_map(|h| f(h).map(|v| (v, h.common))))
    };

    Some(HeadToHead {
        peer: first.peer.clone(),
        common,
        only_self: present.iter().map(|h| h.only_self).sum(),
        only_peer: present.iter().map(|h| h.only_peer).sum(),
        same_locus: weighted_by_common(|h| h.same_locus),
        better: weighted_by_common(|h| h.better),
        equal: weighted_by_common(|h| h.equal),
        worse: weighted_by_common(|h| h.worse),
        nm_delta: weighted_by_common(|h| h.nm_delta),
    })
}

/// Pools clip geometry across the samples that have it.
fn pool_clip(group: &[&SampleResult]) -> Option<ClipGeometry> {
    let mut pooled: Option<ClipGeometry> = None;
    for geometry in group.iter().filter_map(|r| r.clip.as_ref()) {
        match &mut pooled {
            Some(total) => total.add(geometry),
            None => pooled = Some(geometry.clone()),
        }
    }
    pooled
}

/// Weighted mean, or `None` when nothing carries weight.
fn weighted(values: impl Iterator<Item = (f64, u64)>) -> Option<f64> {
    let (mut sum, mut weight) = (0.0, 0u64);
    for (value, w) in values {
        sum += value * w as f64;
        weight += w;
    }
    (weight > 0).then(|| sum / weight as f64)
}

/// Sample standard deviation (Bessel-corrected).
fn sample_sd(values: &[f64]) -> f64 {
    if values.len() < 2 {
        return 0.0;
    }
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let variance =
        values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (values.len() - 1) as f64;
    variance.sqrt()
}

/// Builds the `<prefix>.align_summary.tsv` frame.
///
/// # Arguments
///
/// * `summaries` - Pooled per-contender results
///
/// # Returns
///
/// One row per `(dataset, tool)`.
///
/// # Errors
///
/// Returns [`AlignError::Output`] when the frame cannot be assembled.
pub fn summary_frame(summaries: &[ToolSummary]) -> AlignResult<DataFrame> {
    let stratified: Vec<Option<u64>> = summaries.iter().map(|s| s.score.reference).collect();
    let has_strata = stratified.iter().any(Option::is_some);
    let best: Vec<Option<mapq::MapqPoint>> = summaries.iter().map(|s| s.best_f1()).collect();

    let mut columns = vec![
        Series::new("dataset".into(), strings(summaries, |s| s.dataset.clone())),
        Series::new("tool".into(), strings(summaries, |s| s.tool.clone())),
        Series::new(
            "scoring".into(),
            strings(summaries, |s| scoring_label(s.scoring).to_string()),
        ),
        Series::new(
            "samples".into(),
            summaries
                .iter()
                .map(|s| s.samples as u64)
                .collect::<Vec<_>>(),
        ),
        Series::new("total".into(), numbers(summaries, |s| s.score.total)),
        Series::new("aligned".into(), numbers(summaries, |s| s.score.aligned)),
        Series::new("correct".into(), numbers(summaries, |s| s.score.correct)),
        Series::new("mappable".into(), numbers(summaries, |s| s.mappable)),
        Series::new(
            "align_pct".into(),
            floats(summaries, |s| s.score.align_pct()),
        ),
        Series::new(
            "correct_pct".into(),
            floats(summaries, |s| s.score.correct_pct()),
        ),
        Series::new(
            "correct_pct_sd".into(),
            summaries
                .iter()
                .map(|s| s.correct_pct_sd)
                .collect::<Vec<_>>(),
        ),
        Series::new("recall_pct".into(), floats(summaries, |s| s.recall_pct())),
    ];

    if has_strata {
        // The strata are shares of ALL reads in the truth, as report.py renders them
        // ("reference-correct (of total)", "position-correct (of total)"). Dividing by `aligned`
        // instead would make a tool look better for placing fewer reads, which is the trap
        // `correct_pct` already carries and why it is reported beside `recall_pct`.
        columns.push(Series::new(
            "reference_pct".into(),
            summaries
                .iter()
                .map(|s| s.score.reference.map(|v| pct(v, s.score.total)))
                .collect::<Vec<_>>(),
        ));
        columns.push(Series::new(
            "position_pct".into(),
            summaries
                .iter()
                .map(|s| s.score.position.map(|v| pct(v, s.score.total)))
                .collect::<Vec<_>>(),
        ));
        // ...and the precision counterpart the Python reports alongside them.
        columns.push(Series::new(
            "position_precision_pct".into(),
            summaries
                .iter()
                .map(|s| s.score.position.map(|v| pct(v, s.score.aligned)))
                .collect::<Vec<_>>(),
        ));
    }

    // Base-level metrics. A column stays entirely null when no contender defines it, which is not
    // the same as zero -- a tool that emits no NM has no reported identity.
    let base = |f: fn(&BaseMetrics) -> Option<f64>| -> Vec<Option<f64>> {
        summaries
            .iter()
            .map(|s| s.base.as_ref().and_then(f))
            .collect()
    };
    let base_share = |f: fn(&BaseMetrics) -> f64| -> Vec<Option<f64>> {
        summaries.iter().map(|s| s.base.as_ref().map(f)).collect()
    };
    let h2h = |f: fn(&HeadToHead) -> Option<f64>| -> Vec<Option<f64>> {
        summaries
            .iter()
            .map(|s| s.h2h.as_ref().and_then(f))
            .collect()
    };
    let h2h_count = |f: fn(&HeadToHead) -> u64| -> Vec<Option<u64>> {
        summaries.iter().map(|s| s.h2h.as_ref().map(f)).collect()
    };

    columns.extend([
        // Every aln_* column below comes from `base`, so they all share one denominator: the
        // alignments the tool emitted. Deriving some of them from `score.aligned` instead -- the
        // truth reads it placed -- silently mixes two denominators in one row, and produces
        // percentages over 100 as soon as a tool places a read the truth does not know about.
        Series::new(
            "aln_records".into(),
            summaries
                .iter()
                .map(|s| s.base.as_ref().map(|b| b.records))
                .collect::<Vec<Option<u64>>>(),
        ),
        Series::new(
            "aln_unmapped".into(),
            numbers(summaries, |s| s.counters.unmapped),
        ),
        Series::new(
            "aln_secondary".into(),
            numbers(summaries, |s| s.counters.secondary),
        ),
        Series::new(
            "aln_verified".into(),
            summaries
                .iter()
                .map(|s| s.base.as_ref().map(|b| b.verified))
                .collect::<Vec<Option<u64>>>(),
        ),
        Series::new("aln_identity".into(), base(|b| b.identity)),
        Series::new("aln_identity_high".into(), base(|b| b.identity_high)),
        Series::new("aln_identity_v".into(), base(|b| b.identity_v)),
        Series::new("aln_identity_high_v".into(), base(|b| b.identity_high_v)),
        Series::new("aln_nm_agree".into(), base(|b| b.nm_agree)),
        Series::new("aln_identity_correct".into(), base(|b| b.identity_correct)),
        Series::new("aln_identity_wrong".into(), base(|b| b.identity_wrong)),
        Series::new("aln_coverage".into(), base_share(|b| b.coverage)),
        Series::new("aln_full_length".into(), base_share(|b| b.full_length)),
        Series::new("aln_indel".into(), base_share(|b| b.indel)),
        Series::new("aln_malformed".into(), base_share(|b| b.malformed)),
        Series::new("aln_unknown_ref".into(), base(|b| b.unknown_ref)),
        Series::new("aln_no_nm".into(), base_share(|b| b.no_nm)),
        Series::new("aln_proper_pair".into(), base_share(|b| b.proper_pair)),
        Series::new(
            "h2h_peer".into(),
            summaries
                .iter()
                .map(|s| s.h2h.as_ref().map(|h| h.peer.clone()))
                .collect::<Vec<Option<String>>>(),
        ),
        Series::new("h2h_common".into(), h2h_count(|h| h.common)),
        Series::new("h2h_only_self".into(), h2h_count(|h| h.only_self)),
        Series::new("h2h_only_peer".into(), h2h_count(|h| h.only_peer)),
        Series::new("h2h_same_locus".into(), h2h(|h| h.same_locus)),
        Series::new("h2h_better".into(), h2h(|h| h.better)),
        Series::new("h2h_equal".into(), h2h(|h| h.equal)),
        Series::new("h2h_worse".into(), h2h(|h| h.worse)),
        Series::new("h2h_nm_delta".into(), h2h(|h| h.nm_delta)),
        Series::new(
            "clip_dovetail_pct".into(),
            summaries
                .iter()
                .map(|s| s.clip.as_ref().map(|c| c.dovetail_pct()))
                .collect::<Vec<Option<f64>>>(),
        ),
        Series::new(
            "clip_contained_pct".into(),
            summaries
                .iter()
                .map(|s| s.clip.as_ref().map(|c| c.contained_pct()))
                .collect::<Vec<Option<f64>>>(),
        ),
        Series::new(
            "clip_dovetail_mean_bases".into(),
            summaries
                .iter()
                .map(|s| s.clip.as_ref().and_then(|c| c.dovetail_mean_bases))
                .collect::<Vec<Option<f64>>>(),
        ),
        Series::new(
            "clip_contained_mean_bases".into(),
            summaries
                .iter()
                .map(|s| s.clip.as_ref().and_then(|c| c.contained_mean_bases))
                .collect::<Vec<Option<f64>>>(),
        ),
        Series::new(
            "clip_unknown_contig".into(),
            summaries
                .iter()
                .map(|s| s.clip.as_ref().map(|c| c.unknown))
                .collect::<Vec<Option<u64>>>(),
        ),
        Series::new(
            "mapq_best_cutoff".into(),
            best.iter()
                .map(|b| b.map(|p| p.mapq as u32))
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "mapq_best_f1".into(),
            best.iter().map(|b| b.map(|p| p.f1())).collect::<Vec<_>>(),
        ),
        Series::new(
            "note".into(),
            summaries
                .iter()
                .map(|s| s.note.clone())
                .collect::<Vec<Option<String>>>(),
        ),
    ]);

    DataFrame::new(columns).map_err(|e| AlignError::Output {
        path: "align_summary.tsv".into(),
        message: e.to_string(),
    })
}

/// Builds the `<prefix>.align_samples.tsv` frame: the rows behind every aggregate.
///
/// # Arguments
///
/// * `results` - Per-sample results
///
/// # Returns
///
/// One row per `(dataset, sample, tool)`.
///
/// # Errors
///
/// Returns [`AlignError::Output`] when the frame cannot be assembled.
pub fn samples_frame(results: &[SampleResult]) -> AlignResult<DataFrame> {
    let has_strata = results.iter().any(|r| r.score.reference.is_some());

    let mut columns = vec![
        Series::new("dataset".into(), strings(results, |r| r.dataset.clone())),
        Series::new("sample".into(), strings(results, |r| r.sample.clone())),
        Series::new("tool".into(), strings(results, |r| r.tool.clone())),
        Series::new("total".into(), numbers(results, |r| r.score.total)),
        Series::new("aligned".into(), numbers(results, |r| r.score.aligned)),
        Series::new("correct".into(), numbers(results, |r| r.score.correct)),
        Series::new("mappable".into(), numbers(results, |r| r.mappable)),
        Series::new("align_pct".into(), floats(results, |r| r.score.align_pct())),
        Series::new(
            "correct_pct".into(),
            floats(results, |r| r.score.correct_pct()),
        ),
        Series::new("recall_pct".into(), floats(results, |r| r.recall_pct())),
    ];

    if has_strata {
        columns.push(Series::new(
            "reference".into(),
            results
                .iter()
                .map(|r| r.score.reference)
                .collect::<Vec<_>>(),
        ));
        columns.push(Series::new(
            "position".into(),
            results.iter().map(|r| r.score.position).collect::<Vec<_>>(),
        ));
    }

    let base = |f: fn(&BaseMetrics) -> Option<f64>| -> Vec<Option<f64>> {
        results
            .iter()
            .map(|r| r.base.as_ref().and_then(f))
            .collect()
    };

    columns.extend([
        Series::new("unmapped".into(), numbers(results, |r| r.counters.unmapped)),
        Series::new(
            "secondary".into(),
            numbers(results, |r| r.counters.secondary),
        ),
        Series::new("no_nm".into(), numbers(results, |r| r.counters.no_nm)),
        Series::new(
            "proper_pair".into(),
            numbers(results, |r| r.counters.proper_pair),
        ),
        Series::new(
            "aln_verified".into(),
            results
                .iter()
                .map(|r| r.base.as_ref().map(|b| b.verified))
                .collect::<Vec<Option<u64>>>(),
        ),
        Series::new("aln_identity".into(), base(|b| b.identity)),
        Series::new("aln_identity_v".into(), base(|b| b.identity_v)),
        Series::new("aln_nm_agree".into(), base(|b| b.nm_agree)),
        Series::new(
            "h2h_nm_delta".into(),
            results
                .iter()
                .map(|r| r.h2h.as_ref().and_then(|h| h.nm_delta))
                .collect::<Vec<Option<f64>>>(),
        ),
    ]);

    DataFrame::new(columns).map_err(|e| AlignError::Output {
        path: "align_samples.tsv".into(),
        message: e.to_string(),
    })
}

/// Builds the `<prefix>.align_mapq.tsv` frame: the pooled precision/recall curve per contender.
///
/// # Arguments
///
/// * `summaries` - Pooled per-contender results
///
/// # Returns
///
/// One row per `(dataset, tool, mapq)`.
///
/// # Errors
///
/// Returns [`AlignError::Output`] when the frame cannot be assembled.
pub fn mapq_frame(summaries: &[ToolSummary]) -> AlignResult<DataFrame> {
    let mut dataset = Vec::new();
    let mut tool = Vec::new();
    let mut cutoff = Vec::new();
    let mut correct = Vec::new();
    let mut kept = Vec::new();
    let mut recall = Vec::new();
    let mut precision = Vec::new();
    let mut f1 = Vec::new();

    for summary in summaries {
        let points = mapq::curve(&summary.mapq_counts, summary.mappable);
        for (count, point) in summary.mapq_counts.iter().zip(points.iter()) {
            dataset.push(summary.dataset.clone());
            tool.push(summary.tool.clone());
            cutoff.push(point.mapq as u32);
            correct.push(count.correct);
            kept.push(count.kept);
            recall.push(point.recall_pct);
            precision.push(point.precision_pct);
            f1.push(point.f1());
        }
    }

    DataFrame::new(vec![
        Series::new("dataset".into(), dataset),
        Series::new("tool".into(), tool),
        Series::new("mapq".into(), cutoff),
        Series::new("correct".into(), correct),
        Series::new("kept".into(), kept),
        Series::new("recall_pct".into(), recall),
        Series::new("precision_pct".into(), precision),
        Series::new("f1".into(), f1),
    ])
    .map_err(|e| AlignError::Output {
        path: "align_mapq.tsv".into(),
        message: e.to_string(),
    })
}

/// Builds the `<prefix>.align_reads.tsv` frame: one row per truth read per contender.
///
/// # Arguments
///
/// * `rows` - Per-read verdicts, each tagged with the run that produced it
///
/// # Returns
///
/// The per-read table.
///
/// # Errors
///
/// Returns [`AlignError::Output`] when the frame cannot be assembled.
pub fn reads_frame(
    rows: &[(String, String, String, ScoringMode, ReadVerdict)],
) -> AlignResult<DataFrame> {
    DataFrame::new(vec![
        Series::new("dataset".into(), strings(rows, |r| r.0.clone())),
        Series::new("sample".into(), strings(rows, |r| r.1.clone())),
        Series::new("tool".into(), strings(rows, |r| r.2.clone())),
        Series::new("read_id".into(), strings(rows, |r| r.4.key.id.to_string())),
        Series::new(
            "mate".into(),
            rows.iter().map(|r| r.4.key.mate as u32).collect::<Vec<_>>(),
        ),
        Series::new(
            "truth_contig".into(),
            strings(rows, |r| r.4.truth.contig.to_string()),
        ),
        Series::new("truth_pos".into(), numbers(rows, |r| r.4.truth.pos0)),
        Series::new(
            "truth_genome".into(),
            strings(rows, |r| r.4.truth.genome.to_string()),
        ),
        Series::new(
            "aln_target".into(),
            rows.iter()
                .map(|r| r.4.target.as_ref().map(|t| t.to_string()))
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "aln_pos".into(),
            rows.iter().map(|r| r.4.pos0).collect::<Vec<_>>(),
        ),
        Series::new(
            "mapq".into(),
            rows.iter()
                .map(|r| r.4.mapq.map(|q| q as u32))
                .collect::<Vec<_>>(),
        ),
        Series::new("nm".into(), rows.iter().map(|r| r.4.nm).collect::<Vec<_>>()),
        Series::new(
            "vnm".into(),
            rows.iter().map(|r| r.4.vnm).collect::<Vec<_>>(),
        ),
        Series::new(
            "verdict".into(),
            strings(rows, |r| r.4.verdict.label(r.3).to_string()),
        ),
    ])
    .map_err(|e| AlignError::Output {
        path: "align_reads.tsv".into(),
        message: e.to_string(),
    })
}

/// Writes a frame as a TSV.
///
/// # Arguments
///
/// * `frame` - Table to write
/// * `path` - Destination
///
/// # Returns
///
/// `Ok(())` once written.
///
/// # Errors
///
/// Returns [`AlignError::Output`] when the file cannot be created or serialised.
pub fn write(frame: &mut DataFrame, path: &Path) -> AlignResult<()> {
    write_df(frame, path).map_err(|e| AlignError::Output {
        path: path.to_path_buf(),
        message: e.to_string(),
    })?;
    info!("Wrote {} ({} rows)", path.display(), frame.height());
    Ok(())
}

/// The scoring mode as it appears in the output.
fn scoring_label(scoring: ScoringMode) -> &'static str {
    match scoring {
        ScoringMode::Full => "full",
        ScoringMode::Species => "species",
    }
}

/// Collects a string column.
fn strings<T>(items: &[T], f: impl Fn(&T) -> String) -> Vec<String> {
    items.iter().map(f).collect()
}

/// Collects a `u64` column.
fn numbers<T>(items: &[T], f: impl Fn(&T) -> u64) -> Vec<u64> {
    items.iter().map(f).collect()
}

/// Collects an `f64` column.
fn floats<T>(items: &[T], f: impl Fn(&T) -> f64) -> Vec<f64> {
    items.iter().map(f).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn result(
        dataset: &str,
        sample: &str,
        tool: &str,
        total: u64,
        aligned: u64,
        correct: u64,
    ) -> SampleResult {
        SampleResult {
            dataset: dataset.to_string(),
            sample: sample.to_string(),
            tool: tool.to_string(),
            scoring: ScoringMode::Species,
            score: MappingScore {
                total,
                aligned,
                correct,
                reference: None,
                position: None,
            },
            counters: ParseCounters::default(),
            mapq_counts: vec![MapqCount {
                mapq: 0,
                correct,
                kept: aligned,
            }],
            mappable: correct,
            base: None,
            h2h: None,
            clip: None,
        }
    }

    #[test]
    fn common_sample_restriction_drops_the_odd_sample_out() {
        let results = vec![
            result("ds", "s1", "a", 10, 10, 10),
            result("ds", "s2", "a", 10, 10, 10),
            result("ds", "s1", "b", 10, 10, 5),
        ];
        let (kept, notes) = restrict_to_common_samples(results);

        assert_eq!(kept.len(), 2, "s2 has no 'b' result, so it is dropped");
        assert!(kept.iter().all(|r| r.sample == "s1"));
        assert!(notes["ds"].contains("dropped s2"), "{:?}", notes);
    }

    #[test]
    fn a_single_contender_keeps_every_sample() {
        let results = vec![
            result("ds", "s1", "a", 10, 10, 10),
            result("ds", "s2", "a", 10, 10, 10),
        ];
        let (kept, notes) = restrict_to_common_samples(results);
        assert_eq!(kept.len(), 2);
        assert!(notes.is_empty());
    }

    #[test]
    fn summaries_pool_counts_rather_than_averaging_rates() {
        // 50% of 100 and 100% of 100 pool to 75%, which is also the mean here -- so make the
        // sample sizes differ, where a mean of rates would give the wrong answer.
        let results = vec![
            result("ds", "s1", "a", 10, 10, 10),      // 100% of 10
            result("ds", "s2", "a", 1000, 1000, 500), // 50% of 1000
        ];
        let (kept, notes) = restrict_to_common_samples(results);
        let summaries = summarize(&kept, &notes);

        assert_eq!(summaries.len(), 1);
        assert_eq!(summaries[0].samples, 2);
        assert_eq!(summaries[0].score.aligned, 1010);
        assert_eq!(summaries[0].score.correct, 510);
        let pooled = summaries[0].score.correct_pct();
        assert!((pooled - 100.0 * 510.0 / 1010.0).abs() < 1e-9);
        assert!(
            (pooled - 75.0).abs() > 1.0,
            "a mean of rates would have said 75%"
        );
    }

    #[test]
    fn summaries_report_spread_only_when_there_is_more_than_one_sample() {
        let one = summarize(&[result("ds", "s1", "a", 10, 10, 10)], &HashMap::new());
        assert_eq!(one[0].correct_pct_sd, None);

        let two = summarize(
            &[
                result("ds", "s1", "a", 10, 10, 10),
                result("ds", "s2", "a", 10, 10, 5),
            ],
            &HashMap::new(),
        );
        assert!(two[0].correct_pct_sd.unwrap() > 0.0);
    }

    #[test]
    fn frames_have_the_expected_shape() {
        let results = vec![
            result("ds", "s1", "a", 10, 10, 10),
            result("ds", "s1", "b", 10, 8, 4),
        ];
        let (kept, notes) = restrict_to_common_samples(results);
        let summaries = summarize(&kept, &notes);

        let summary = summary_frame(&summaries).unwrap();
        assert_eq!(summary.height(), 2);
        assert!(summary
            .get_column_names()
            .iter()
            .any(|n| n.as_str() == "correct_pct"));
        // species scoring has no reference/position strata, so those columns stay out entirely
        assert!(!summary
            .get_column_names()
            .iter()
            .any(|n| n.as_str() == "position_pct"));

        assert_eq!(samples_frame(&kept).unwrap().height(), 2);
        assert_eq!(mapq_frame(&summaries).unwrap().height(), 2);
    }

    #[test]
    fn full_scoring_adds_the_stratum_columns() {
        let mut results = vec![result("ds", "s1", "a", 10, 10, 10)];
        results[0].scoring = ScoringMode::Full;
        results[0].score.reference = Some(8);
        results[0].score.position = Some(6);

        let summaries = summarize(&results, &HashMap::new());
        let frame = summary_frame(&summaries).unwrap();
        let names: Vec<String> = frame
            .get_column_names()
            .iter()
            .map(|n| n.to_string())
            .collect();

        assert!(names.contains(&"reference_pct".to_string()));
        assert!(names.contains(&"position_pct".to_string()));
    }

    #[test]
    fn sd_is_bessel_corrected() {
        assert_eq!(sample_sd(&[1.0]), 0.0);
        assert!((sample_sd(&[2.0, 4.0]) - 2.0f64.sqrt()).abs() < 1e-12);
    }
}
