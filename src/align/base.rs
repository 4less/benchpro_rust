//! Base-level scoring: is the alignment itself any good?
//!
//! The mapping score only asks *where* a read went. Once a tool emits real alignments the
//! interesting questions are about the alignment:
//!
//! * is the record well formed? — the CIGAR must consume exactly the read
//! * does it cover the read? — query coverage, the rest is clipped away
//! * how close to the reference? — identity, preferably replayed rather than self-reported
//! * is it the *best* alignment? — a head-to-head against a peer aligner at the same locus
//!
//! On a marker database there is no base-level truth: truth knows the read's source species, not
//! which marker locus it should land on. So quality is judged two ways — identity split by
//! right-versus-wrong target (a healthy aligner's wrong-target hits should be its worse ones), and
//! the head-to-head on the loci both tools found.

use std::collections::HashMap;

use super::metrics::{pct, ScoringContext};
use super::sam::{AlnRecord, ParseCounters};
use super::truth::{ReadKey, Truth};

/// Identity at or above which an alignment counts as near-perfect.
pub const IDENTITY_HIGH: f64 = 0.95;

/// Base-level metrics of one contender on one sample.
///
/// Every field is a mean or a share over the emitted primary alignments, except the counts. `None`
/// means the quantity is not defined for this tool — a tool that emits no `NM` has no reported
/// identity, and that must not be rendered as zero.
#[derive(Debug, Clone, Default)]
pub struct BaseMetrics {
    /// Primary alignments emitted.
    pub records: u64,
    /// Alignments replayed against the reference.
    pub verified: u64,
    /// Mean identity as the tool reports it.
    pub identity: Option<f64>,
    /// Share of reported identities at or above [`IDENTITY_HIGH`].
    pub identity_high: Option<f64>,
    /// Mean identity from the replayed edit distance.
    pub identity_v: Option<f64>,
    /// Share of replayed identities at or above [`IDENTITY_HIGH`].
    pub identity_high_v: Option<f64>,
    /// Share of replayed records whose `NM` tag matches the replayed distance.
    pub nm_agree: Option<f64>,
    /// Mean identity of alignments on the right target.
    pub identity_correct: Option<f64>,
    /// Mean identity of alignments on the wrong target.
    pub identity_wrong: Option<f64>,
    /// Mean query coverage.
    pub coverage: f64,
    /// Share of unclipped alignments.
    pub full_length: f64,
    /// Share of alignments containing an indel.
    pub indel: f64,
    /// Share of malformed records (CIGAR does not match SEQ length).
    pub malformed: f64,
    /// Share whose target is absent from the reference.
    pub unknown_ref: Option<f64>,
    /// Share missing an `NM` tag.
    pub no_nm: f64,
    /// Share flagged as a proper pair.
    pub proper_pair: f64,
}

/// Computes the base-level metrics of one contender's alignments.
///
/// # Arguments
///
/// * `records` - Primary alignments, keyed by read mate
/// * `counters` - Parse counters, for the shares that come from skipped records
/// * `truth` - The sample's truth, used only to split identity by right-vs-wrong target
/// * `context` - Correctness definition, for that same split
///
/// # Returns
///
/// The metrics, or `None` when the tool emitted no alignment at all.
pub fn summarize(
    records: &HashMap<ReadKey, AlnRecord>,
    counters: &ParseCounters,
    truth: Option<&Truth>,
    context: &ScoringContext,
) -> Option<BaseMetrics> {
    let total = records.len() as u64;
    if total == 0 {
        return None;
    }

    let mut metrics = BaseMetrics {
        records: total,
        coverage: 100.0 * mean(records.values().map(AlnRecord::coverage)).unwrap_or(0.0),
        full_length: pct(
            records.values().filter(|r| r.counts.clip() == 0).count() as u64,
            total,
        ),
        indel: pct(
            records.values().filter(|r| r.has_indel()).count() as u64,
            total,
        ),
        malformed: pct(
            records.values().filter(|r| r.malformed).count() as u64,
            total,
        ),
        no_nm: pct(counters.no_nm, total),
        proper_pair: pct(counters.proper_pair, total),
        ..Default::default()
    };

    // Identity as the tool reports it, over the records that carry an NM tag.
    let reported: Vec<f64> = records.values().filter_map(AlnRecord::identity).collect();
    metrics.identity = mean(reported.iter().copied()).map(|v| 100.0 * v);
    metrics.identity_high = share_at_least(&reported, IDENTITY_HIGH);

    // The same alignments scored against the reference bases rather than the tool's bookkeeping.
    let replayed: Vec<&AlnRecord> = records.values().filter(|r| r.vnm.is_some()).collect();
    if !replayed.is_empty() {
        let identities: Vec<f64> = replayed
            .iter()
            .filter_map(|r| r.identity_verified())
            .collect();
        metrics.verified = replayed.len() as u64;
        metrics.identity_v = mean(identities.iter().copied()).map(|v| 100.0 * v);
        metrics.identity_high_v = share_at_least(&identities, IDENTITY_HIGH);

        // The tool's NM checked against the replay: a disagreement means the record's edit
        // distance is not the alignment it describes.
        let comparable: Vec<&&AlnRecord> = replayed.iter().filter(|r| r.nm.is_some()).collect();
        if !comparable.is_empty() {
            let agreeing = comparable
                .iter()
                .filter(|r| r.nm == r.vnm.map(|v| v as i64))
                .count();
            metrics.nm_agree = Some(pct(agreeing as u64, comparable.len() as u64));
        }
    }

    if counters.unknown_ref > 0 {
        metrics.unknown_ref = Some(pct(counters.unknown_ref, total));
    }

    // Identity split by whether the alignment is on the right target. Wrong-target alignments
    // should be the worse ones; if they are not, the score used to pick the target is not tracking
    // alignment quality. The replayed identity is preferred, which every tool has once verified.
    if let Some(truth) = truth {
        let (mut right, mut wrong) = (Vec::new(), Vec::new());
        for (key, record) in records {
            let Some(entry) = truth.get(key) else {
                continue;
            };
            let Some(identity) = record.identity_verified().or_else(|| record.identity()) else {
                continue;
            };
            if context.verdict(record, entry).is_correct() {
                right.push(identity);
            } else {
                wrong.push(identity);
            }
        }
        metrics.identity_correct = mean(right.into_iter()).map(|v| 100.0 * v);
        metrics.identity_wrong = mean(wrong.into_iter()).map(|v| 100.0 * v);
    }

    Some(metrics)
}

/// Base-level comparison against a peer aligner on the loci both tools found.
#[derive(Debug, Clone, Default)]
pub struct HeadToHead {
    /// The peer's name.
    pub peer: String,
    /// Read mates both tools placed.
    pub common: u64,
    /// Mates only this tool placed.
    pub only_self: u64,
    /// Mates only the peer placed.
    pub only_peer: u64,
    /// Share of shared mates placed at the same locus.
    pub same_locus: Option<f64>,
    /// Share where this tool reached the lower edit distance.
    pub better: Option<f64>,
    /// Share where the two tied.
    pub equal: Option<f64>,
    /// Share where this tool was worse.
    pub worse: Option<f64>,
    /// Mean `self_nm - peer_nm`; negative means this tool aligns better.
    pub nm_delta: Option<f64>,
}

/// Compares one tool's alignments against a peer's at the loci both found.
///
/// Only mates where both landed on the same reference within `tolerance` bp are compared —
/// elsewhere the two are solving different problems and their edit distances are not comparable.
///
/// Replayed edit distances are used where *both* records were verified and self-reported ones
/// otherwise, never a mix: a mixed comparison would measure the bookkeeping rather than the
/// alignment.
///
/// # Arguments
///
/// * `records` - This tool's primary alignments
/// * `peer` - The peer's primary alignments
/// * `peer_name` - The peer's name, for the output
/// * `tolerance` - bp of slack for "same locus"
///
/// # Returns
///
/// The head-to-head summary.
pub fn head_to_head(
    records: &HashMap<ReadKey, AlnRecord>,
    peer: &HashMap<ReadKey, AlnRecord>,
    peer_name: &str,
    tolerance: u64,
) -> HeadToHead {
    let mut result = HeadToHead {
        peer: peer_name.to_string(),
        ..Default::default()
    };

    let (mut same, mut better, mut equal, mut worse) = (0u64, 0u64, 0u64, 0u64);
    let mut delta_sum = 0i64;

    for (key, record) in records {
        let Some(other) = peer.get(key) else { continue };
        result.common += 1;
        if record.target != other.target || record.pos0.abs_diff(other.pos0) > tolerance {
            continue;
        }
        same += 1;

        let (mine, theirs) = match (record.vnm, other.vnm) {
            (Some(a), Some(b)) => (Some(a as i64), Some(b as i64)),
            _ => (record.nm, other.nm),
        };
        let (Some(mine), Some(theirs)) = (mine, theirs) else {
            continue;
        };
        delta_sum += mine - theirs;
        match mine.cmp(&theirs) {
            std::cmp::Ordering::Less => better += 1,
            std::cmp::Ordering::Equal => equal += 1,
            std::cmp::Ordering::Greater => worse += 1,
        }
    }

    result.only_self = records.len() as u64 - result.common;
    result.only_peer = peer.len() as u64 - result.common;
    if result.common > 0 {
        result.same_locus = Some(pct(same, result.common));
    }

    let compared = better + equal + worse;
    if compared > 0 {
        result.better = Some(pct(better, compared));
        result.equal = Some(pct(equal, compared));
        result.worse = Some(pct(worse, compared));
        result.nm_delta = Some(delta_sum as f64 / compared as f64);
    }

    result
}

/// Fixed-point scale for [`mean`]. 2^52 keeps ~15 significant digits for the [0, 1] ratios these
/// means are taken over, which is f64's own precision there.
const MEAN_SCALE: f64 = (1u64 << 52) as f64;

/// Mean of an iterator, or `None` when it is empty.
///
/// Accumulates in fixed point rather than f64. These means are taken over `records.values()`, and a
/// `HashMap` iterates in an order seeded per process, so a floating-point sum would give a slightly
/// different answer on every run — f64 addition is not associative. Integer addition is, so the
/// result depends only on the multiset of values, which is what makes a run reproducible.
fn mean(values: impl Iterator<Item = f64>) -> Option<f64> {
    let mut sum: i128 = 0;
    let mut count = 0u64;
    for value in values {
        sum += (value * MEAN_SCALE) as i128;
        count += 1;
    }
    (count > 0).then(|| sum as f64 / MEAN_SCALE / count as f64)
}

/// Share of values at or above a threshold, as a percentage; `None` for an empty slice.
fn share_at_least(values: &[f64], threshold: f64) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    Some(pct(
        values.iter().filter(|v| **v >= threshold).count() as u64,
        values.len() as u64,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::cigar;
    use crate::align::truth::TruthEntry;
    use crate::options::ScoringMode;

    fn record(
        target: &str,
        pos0: u64,
        cigar_str: &str,
        nm: Option<i64>,
        vnm: Option<u64>,
    ) -> AlnRecord {
        AlnRecord {
            target: target.into(),
            pos0,
            mapq: 60,
            nm,
            counts: cigar::count(cigar_str.as_bytes()),
            clip_ends: cigar::clip_ends(cigar_str.as_bytes()),
            cigar_fingerprint: Some(cigar::fingerprint(cigar_str.as_bytes())),
            malformed: false,
            proper_pair: false,
            cigar: None,
            seq: None,
            vnm,
            offset: 0,
        }
    }

    fn key(id: &str) -> ReadKey {
        ReadKey {
            id: id.into(),
            mate: 1,
        }
    }

    fn context() -> ScoringContext<'static> {
        ScoringContext {
            scoring: ScoringMode::Species,
            sep: "_",
            contig2genome: None,
            tolerance: 100,
        }
    }

    #[test]
    fn mean_does_not_depend_on_summation_order() {
        // The values a HashMap hands out come in a per-process order, so the same multiset summed
        // two ways must still give bit-identical answers.
        let values: Vec<f64> = (0..1000).map(|i| 0.5 + i as f64 / 3000.0).collect();
        let forward = mean(values.iter().copied()).unwrap();
        let backward = mean(values.iter().rev().copied()).unwrap();
        let shuffled: Vec<f64> = values
            .iter()
            .skip(377)
            .chain(values.iter().take(377))
            .copied()
            .collect();

        assert_eq!(forward.to_bits(), backward.to_bits());
        assert_eq!(
            forward.to_bits(),
            mean(shuffled.into_iter()).unwrap().to_bits()
        );
        // ...and it is still the right answer.
        let naive: f64 = values.iter().sum::<f64>() / values.len() as f64;
        assert!((forward - naive).abs() < 1e-12);
    }

    #[test]
    fn mean_of_nothing_is_none() {
        assert_eq!(mean(std::iter::empty()), None);
    }

    #[test]
    fn no_alignments_means_no_metrics_rather_than_zeroes() {
        let empty = HashMap::new();
        assert!(summarize(&empty, &ParseCounters::default(), None, &context()).is_none());
    }

    #[test]
    fn reported_identity_covers_only_the_records_that_carry_nm() {
        let records: HashMap<ReadKey, AlnRecord> = [
            (key("r1"), record("ctg_1", 0, "100M", Some(0), None)),
            (key("r2"), record("ctg_1", 0, "100M", Some(10), None)),
            (key("r3"), record("ctg_1", 0, "100M", None, None)),
        ]
        .into_iter()
        .collect();
        let counters = ParseCounters {
            no_nm: 1,
            ..Default::default()
        };

        let metrics = summarize(&records, &counters, None, &context()).unwrap();

        // Mean of 100% and 90%; the NM-less record contributes nothing rather than a zero.
        assert!((metrics.identity.unwrap() - 95.0).abs() < 1e-9);
        assert!((metrics.identity_high.unwrap() - 50.0).abs() < 1e-9);
        assert!((metrics.no_nm - 100.0 / 3.0).abs() < 1e-9);
        assert!(metrics.identity_v.is_none(), "nothing was replayed");
    }

    #[test]
    fn replayed_identity_and_the_nm_agreement_check() {
        let records: HashMap<ReadKey, AlnRecord> = [
            // tool says 0 mismatches, the reference says 0: agrees
            (key("honest"), record("ctg_1", 0, "100M", Some(0), Some(0))),
            // tool says 0, the reference says 10: its CIGAR describes a worse alignment
            (key("liar"), record("ctg_1", 0, "100M", Some(0), Some(10))),
        ]
        .into_iter()
        .collect();

        let metrics = summarize(&records, &ParseCounters::default(), None, &context()).unwrap();

        assert_eq!(metrics.verified, 2);
        assert!(
            (metrics.identity.unwrap() - 100.0).abs() < 1e-9,
            "as reported"
        );
        assert!(
            (metrics.identity_v.unwrap() - 95.0).abs() < 1e-9,
            "as replayed"
        );
        assert!((metrics.nm_agree.unwrap() - 50.0).abs() < 1e-9);
    }

    #[test]
    fn identity_splits_by_right_and_wrong_target() {
        let truth: Truth = ["r1", "r2"]
            .iter()
            .map(|id| {
                (
                    key(id),
                    TruthEntry {
                        contig: "src".into(),
                        pos0: 0,
                        genome: "1234".into(),
                        cigar_fingerprint: None,
                        shape: None,
                    },
                )
            })
            .collect();
        let records: HashMap<ReadKey, AlnRecord> = [
            // right species, near perfect
            (key("r1"), record("1234_1", 0, "100M", Some(1), None)),
            // wrong species, visibly worse -- which is what a healthy aligner looks like
            (key("r2"), record("9999_1", 0, "100M", Some(20), None)),
        ]
        .into_iter()
        .collect();

        let metrics = summarize(
            &records,
            &ParseCounters::default(),
            Some(&truth),
            &context(),
        )
        .unwrap();

        assert!((metrics.identity_correct.unwrap() - 99.0).abs() < 1e-9);
        assert!((metrics.identity_wrong.unwrap() - 80.0).abs() < 1e-9);
    }

    #[test]
    fn coverage_and_clipping_shares() {
        let records: HashMap<ReadKey, AlnRecord> = [
            (key("clean"), record("ctg_1", 0, "100M", Some(0), None)),
            (key("clipped"), record("ctg_1", 0, "50S50M", Some(0), None)),
            (key("indel"), record("ctg_1", 0, "50M2I48M", Some(0), None)),
        ]
        .into_iter()
        .collect();

        let metrics = summarize(&records, &ParseCounters::default(), None, &context()).unwrap();

        assert!((metrics.full_length - 200.0 / 3.0).abs() < 1e-9);
        assert!((metrics.indel - 100.0 / 3.0).abs() < 1e-9);
        // (100% + 50% + 100%) / 3
        assert!((metrics.coverage - 250.0 / 3.0).abs() < 1e-9);
    }

    #[test]
    fn unknown_ref_is_reported_only_when_it_happened() {
        let records: HashMap<ReadKey, AlnRecord> =
            [(key("r1"), record("ctg_1", 0, "100M", Some(0), None))]
                .into_iter()
                .collect();

        let quiet = summarize(&records, &ParseCounters::default(), None, &context()).unwrap();
        assert!(quiet.unknown_ref.is_none());

        let noisy = summarize(
            &records,
            &ParseCounters {
                unknown_ref: 1,
                ..Default::default()
            },
            None,
            &context(),
        )
        .unwrap();
        assert_eq!(noisy.unknown_ref, Some(100.0));
    }

    #[test]
    fn head_to_head_compares_only_shared_loci() {
        let mine: HashMap<ReadKey, AlnRecord> = [
            (
                key("both_same"),
                record("ctg1", 1000, "100M", Some(1), None),
            ),
            (
                key("both_apart"),
                record("ctg1", 1000, "100M", Some(1), None),
            ),
            (
                key("mine_only"),
                record("ctg1", 1000, "100M", Some(1), None),
            ),
        ]
        .into_iter()
        .collect();
        let theirs: HashMap<ReadKey, AlnRecord> = [
            (
                key("both_same"),
                record("ctg1", 1010, "100M", Some(5), None),
            ),
            // 50 kb away: the two are solving different problems here
            (
                key("both_apart"),
                record("ctg1", 51_000, "100M", Some(0), None),
            ),
            (
                key("theirs_only"),
                record("ctg1", 1000, "100M", Some(1), None),
            ),
        ]
        .into_iter()
        .collect();

        let h2h = head_to_head(&mine, &theirs, "peer", 100);

        assert_eq!(h2h.common, 2);
        assert_eq!(h2h.only_self, 1);
        assert_eq!(h2h.only_peer, 1);
        assert_eq!(h2h.same_locus, Some(50.0));
        // Only the shared locus is compared, and there this tool wins 1 vs 5.
        assert_eq!(h2h.better, Some(100.0));
        assert_eq!(h2h.nm_delta, Some(-4.0));
    }

    #[test]
    fn head_to_head_never_mixes_replayed_and_reported_distances() {
        // This tool was replayed, the peer was not. Comparing vnm against nm would measure the
        // bookkeeping, so both fall back to their reported NM.
        let mine: HashMap<ReadKey, AlnRecord> =
            [(key("r1"), record("ctg1", 0, "100M", Some(9), Some(1)))]
                .into_iter()
                .collect();
        let theirs: HashMap<ReadKey, AlnRecord> =
            [(key("r1"), record("ctg1", 0, "100M", Some(9), None))]
                .into_iter()
                .collect();

        let h2h = head_to_head(&mine, &theirs, "peer", 100);
        assert_eq!(h2h.equal, Some(100.0), "9 vs 9, not 1 vs 9");
        assert_eq!(h2h.nm_delta, Some(0.0));
    }

    #[test]
    fn head_to_head_skips_pairs_with_no_edit_distance_at_all() {
        let mine: HashMap<ReadKey, AlnRecord> =
            [(key("r1"), record("ctg1", 0, "100M", None, None))]
                .into_iter()
                .collect();
        let theirs = mine.clone();

        let h2h = head_to_head(&mine, &theirs, "peer", 100);
        assert_eq!(h2h.common, 1);
        assert_eq!(h2h.same_locus, Some(100.0));
        assert_eq!(h2h.better, None, "nothing comparable");
        assert_eq!(h2h.nm_delta, None);
    }
}
