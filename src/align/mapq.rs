//! Precision and recall as a MAPQ cutoff slides from 0 to the highest MAPQ emitted.
//!
//! At threshold `q` only alignments with `MAPQ >= q` are kept, so recall is what a stricter filter
//! costs you and precision is what it buys. A tool whose MAPQ is informative trades a little recall
//! for a lot of precision; one whose MAPQ is noise moves along a flat line.
//!
//! The counts come from the same pass that produces the mapping score, so the curve and the
//! headline numbers cannot disagree — and they cover PAF contenders too, not only the tools that
//! happen to emit SAM.

use std::collections::HashMap;

use super::metrics::{pct, ScoringContext};
use super::sam::AlnRecord;
use super::truth::{ReadKey, Truth};

/// Cumulative counts at one MAPQ cutoff: how many alignments survive it, and how many are right.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MapqCount {
    /// The cutoff.
    pub mapq: u8,
    /// Correct alignments with `MAPQ >= mapq`.
    pub correct: u64,
    /// All alignments with `MAPQ >= mapq`.
    pub kept: u64,
}

/// One plotted point of the curve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MapqPoint {
    /// The cutoff.
    pub mapq: u8,
    /// `100 * correct / mappable_base`.
    pub recall_pct: f64,
    /// `100 * correct / kept`.
    pub precision_pct: f64,
    /// Alignments surviving the cutoff.
    pub kept: u64,
}

impl MapqPoint {
    /// Harmonic mean of recall and precision.
    ///
    /// # Returns
    ///
    /// The F1 score, or 0 when both are 0.
    pub fn f1(&self) -> f64 {
        let sum = self.recall_pct + self.precision_pct;
        if sum == 0.0 {
            0.0
        } else {
            2.0 * self.recall_pct * self.precision_pct / sum
        }
    }
}

/// Counts correct and kept alignments cumulatively from the highest MAPQ downward.
///
/// Raw counts rather than percentages, because recall's denominator is the *mappable* read count,
/// which is not known until every contender on the sample has been parsed (see [`mappable_base`]).
///
/// # Arguments
///
/// * `records` - Primary alignments, keyed by read mate
/// * `truth` - The sample's per-read truth
/// * `context` - Correctness definition and its parameters
///
/// # Returns
///
/// One entry per cutoff from 0 up to the highest MAPQ emitted, ascending. Empty when there is
/// nothing to score.
pub fn counts(
    records: &HashMap<ReadKey, AlnRecord>,
    truth: &Truth,
    context: &ScoringContext,
) -> Vec<MapqCount> {
    if truth.is_empty() || records.is_empty() {
        return Vec::new();
    }

    let mut kept_at = [0u64; 256];
    let mut correct_at = [0u64; 256];
    let mut any = false;

    for (key, record) in records {
        let Some(entry) = truth.get(key) else {
            continue;
        };
        any = true;
        kept_at[record.mapq as usize] += 1;
        // The MAPQ curve is about placement quality, so under `full` scoring "correct" has to mean
        // the read also landed near its true locus -- otherwise a MAPQ filter would look good for
        // merely picking the right genome.
        let verdict = context.verdict(record, entry);
        let correct = match context.scoring {
            crate::options::ScoringMode::Species => verdict.is_correct(),
            crate::options::ScoringMode::Full => verdict >= super::metrics::Verdict::Position,
        };
        correct_at[record.mapq as usize] += correct as u64;
    }
    if !any {
        return Vec::new();
    }

    let top = kept_at
        .iter()
        .rposition(|&count| count > 0)
        .expect("at least one alignment was counted");

    let mut out = Vec::with_capacity(top + 1);
    let (mut kept, mut correct) = (0u64, 0u64);
    for q in (0..=top).rev() {
        kept += kept_at[q];
        correct += correct_at[q];
        if kept > 0 {
            out.push(MapqCount {
                mapq: q as u8,
                correct,
                kept,
            });
        }
    }
    out.reverse();
    out
}

/// How many reads are actually mappable, taken as the most any *one* contender got right.
///
/// Recall is deliberately not measured against every read in the truth. On a marker-gene reference
/// only a few percent of whole-genome reads can land in a marker at all, so all-reads recall is
/// capped near 5%, dominates F1, and pins every tool's F1-optimal cutoff at 0 — the marker stops
/// discriminating exactly where it should be most useful.
///
/// The base is the maximum *correct* count at MAPQ 0 across contenders, not the union of everything
/// they placed: a tool that aligns almost every read at chance accuracy would otherwise drag the
/// union up to nearly the whole read set and reinstate exactly the ceiling being removed.
/// Best-achieved-correct is the honest empirical answer to "how many of these reads can be placed
/// at all", and it makes the leading tool peak at 100% recall with everyone else measured against
/// it.
///
/// # Arguments
///
/// * `per_tool` - Each contender's counts on the same sample
///
/// # Returns
///
/// The recall denominator, or 0 when no contender placed anything.
pub fn mappable_base<'a>(per_tool: impl Iterator<Item = &'a [MapqCount]>) -> u64 {
    per_tool
        .filter_map(|counts| counts.first().map(|c| c.correct))
        .max()
        .unwrap_or(0)
}

/// Turns cumulative counts into plottable points.
///
/// # Arguments
///
/// * `counts` - Cumulative counts from [`counts`]
/// * `mappable` - Recall denominator from [`mappable_base`]
///
/// # Returns
///
/// One point per cutoff, ascending in MAPQ. Empty when there is no denominator.
pub fn curve(counts: &[MapqCount], mappable: u64) -> Vec<MapqPoint> {
    if mappable == 0 {
        return Vec::new();
    }
    counts
        .iter()
        .map(|c| MapqPoint {
            mapq: c.mapq,
            recall_pct: pct(c.correct, mappable),
            precision_pct: pct(c.correct, c.kept),
            kept: c.kept,
        })
        .collect()
}

/// The cutoff that maximises F1.
///
/// Ties go to the *lower* cutoff: when filtering harder buys nothing, do not throw reads away.
///
/// # Arguments
///
/// * `points` - A curve from [`curve`]
///
/// # Returns
///
/// The best point, or `None` for an empty curve.
pub fn best_f1(points: &[MapqPoint]) -> Option<MapqPoint> {
    let mut best: Option<MapqPoint> = None;
    for point in points {
        if best.is_none_or(|b| point.f1() > b.f1() + 1e-12) {
            best = Some(*point);
        }
    }
    best
}

/// Sums two contenders' or samples' cumulative counts cutoff by cutoff.
///
/// Used to pool a tool's curve across the samples of a dataset: the counts are additive, the rates
/// derived from them are not.
///
/// # Arguments
///
/// * `into` - Accumulator, extended when `from` reaches a higher MAPQ
/// * `from` - Counts to fold in
pub fn add_counts(into: &mut Vec<MapqCount>, from: &[MapqCount]) {
    if from.is_empty() {
        return;
    }
    if into.len() < from.len() {
        // A cumulative count at a cutoff the accumulator never reached is 0 correct of 0 kept.
        for q in into.len()..from.len() {
            into.push(MapqCount {
                mapq: q as u8,
                correct: 0,
                kept: 0,
            });
        }
    }
    for (target, source) in into.iter_mut().zip(from.iter()) {
        target.correct += source.correct;
        target.kept += source.kept;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::cigar::CigarCounts;
    use crate::align::truth::TruthEntry;
    use crate::options::ScoringMode;

    fn record(target: &str, mapq: u8) -> AlnRecord {
        AlnRecord {
            target: target.into(),
            pos0: 0,
            mapq,
            nm: Some(0),
            counts: CigarCounts {
                matches: 100,
                ..Default::default()
            },
            clip_ends: (0, 0),
            malformed: false,
            proper_pair: false,
            cigar: None,
            seq: None,
            vnm: None,
            offset: 0,
        }
    }

    fn key(id: &str) -> ReadKey {
        ReadKey {
            id: id.into(),
            mate: 1,
        }
    }

    fn species_context() -> ScoringContext<'static> {
        ScoringContext {
            scoring: ScoringMode::Species,
            sep: "_",
            contig2genome: None,
            tolerance: 100,
        }
    }

    /// Four reads: two right at high MAPQ, two wrong at low MAPQ.
    fn fixture() -> (HashMap<ReadKey, AlnRecord>, Truth) {
        let truth: Truth = ["r1", "r2", "r3", "r4"]
            .iter()
            .map(|id| {
                (
                    key(id),
                    TruthEntry {
                        contig: "src".into(),
                        pos0: 0,
                        genome: "1234".into(),
                    },
                )
            })
            .collect();
        let records = [
            (key("r1"), record("1234_1", 60)),
            (key("r2"), record("1234_2", 60)),
            (key("r3"), record("9999_1", 0)),
            (key("r4"), record("9999_2", 0)),
        ]
        .into_iter()
        .collect();
        (records, truth)
    }

    #[test]
    fn counts_are_cumulative_from_the_top_down() {
        let (records, truth) = fixture();
        let counts = counts(&records, &truth, &species_context());

        assert_eq!(counts.first().unwrap().mapq, 0);
        assert_eq!(counts.last().unwrap().mapq, 60);
        // At cutoff 0 everything is kept; at 60 only the two correct ones survive.
        assert_eq!(
            counts[0],
            MapqCount {
                mapq: 0,
                correct: 2,
                kept: 4
            }
        );
        assert_eq!(
            *counts.last().unwrap(),
            MapqCount {
                mapq: 60,
                correct: 2,
                kept: 2
            }
        );
    }

    #[test]
    fn an_informative_mapq_trades_recall_for_precision() {
        let (records, truth) = fixture();
        let counts = counts(&records, &truth, &species_context());
        let points = curve(&counts, mappable_base([counts.as_slice()].into_iter()));

        assert_eq!(points[0].precision_pct, 50.0);
        assert_eq!(points[0].recall_pct, 100.0);
        assert_eq!(points.last().unwrap().precision_pct, 100.0);
        assert_eq!(points.last().unwrap().recall_pct, 100.0);
    }

    #[test]
    fn mappable_base_is_the_best_contenders_correct_count_not_the_union() {
        let weak = vec![MapqCount {
            mapq: 0,
            correct: 10,
            kept: 1000,
        }];
        let strong = vec![MapqCount {
            mapq: 0,
            correct: 80,
            kept: 100,
        }];
        // The weak tool placed 1000 reads; that must not become the denominator.
        assert_eq!(
            mappable_base([weak.as_slice(), strong.as_slice()].into_iter()),
            80
        );
    }

    #[test]
    fn mappable_base_of_nothing_is_zero_and_yields_no_curve() {
        assert_eq!(mappable_base(std::iter::empty()), 0);
        assert!(curve(
            &[MapqCount {
                mapq: 0,
                correct: 0,
                kept: 5
            }],
            0
        )
        .is_empty());
    }

    #[test]
    fn recall_is_measured_against_the_leader_not_the_truth_size() {
        // 100 reads in the truth, but the best any tool manages is 40 -- so 40 is 100% recall.
        let counts = vec![MapqCount {
            mapq: 0,
            correct: 40,
            kept: 50,
        }];
        let points = curve(&counts, 40);
        assert_eq!(points[0].recall_pct, 100.0);
        assert_eq!(points[0].precision_pct, 80.0);
    }

    #[test]
    fn best_f1_breaks_ties_toward_the_lower_cutoff() {
        let points = vec![
            MapqPoint {
                mapq: 0,
                recall_pct: 90.0,
                precision_pct: 90.0,
                kept: 10,
            },
            MapqPoint {
                mapq: 1,
                recall_pct: 90.0,
                precision_pct: 90.0,
                kept: 10,
            },
        ];
        assert_eq!(best_f1(&points).unwrap().mapq, 0);
    }

    #[test]
    fn best_f1_picks_the_maximum() {
        let points = vec![
            MapqPoint {
                mapq: 0,
                recall_pct: 100.0,
                precision_pct: 50.0,
                kept: 20,
            },
            MapqPoint {
                mapq: 30,
                recall_pct: 90.0,
                precision_pct: 95.0,
                kept: 11,
            },
            MapqPoint {
                mapq: 60,
                recall_pct: 10.0,
                precision_pct: 100.0,
                kept: 1,
            },
        ];
        assert_eq!(best_f1(&points).unwrap().mapq, 30);
    }

    #[test]
    fn full_scoring_requires_the_right_locus_not_just_the_right_genome() {
        let map: HashMap<Box<str>, Box<str>> = [(Box::from("ctg1"), Box::from("genomeA"))]
            .into_iter()
            .collect();
        let context = ScoringContext {
            scoring: ScoringMode::Full,
            sep: "_",
            contig2genome: Some(&map),
            tolerance: 100,
        };
        let truth: Truth = [(
            key("r1"),
            TruthEntry {
                contig: "ctg1".into(),
                pos0: 10_000,
                genome: "genomeA".into(),
            },
        )]
        .into_iter()
        .collect();
        // Right genome and contig, but 10 kb off: not correct for the curve.
        let records = [(key("r1"), record("ctg1", 60))].into_iter().collect();

        let counts = counts(&records, &truth, &context);
        assert_eq!(counts[0].correct, 0);
        assert_eq!(counts[0].kept, 1);
    }

    #[test]
    fn adding_counts_pools_samples_and_extends_the_range() {
        let mut into = vec![
            MapqCount {
                mapq: 0,
                correct: 5,
                kept: 10,
            },
            MapqCount {
                mapq: 1,
                correct: 4,
                kept: 8,
            },
        ];
        add_counts(
            &mut into,
            &[
                MapqCount {
                    mapq: 0,
                    correct: 1,
                    kept: 2,
                },
                MapqCount {
                    mapq: 1,
                    correct: 1,
                    kept: 1,
                },
                MapqCount {
                    mapq: 2,
                    correct: 1,
                    kept: 1,
                },
            ],
        );

        assert_eq!(into.len(), 3);
        assert_eq!(
            into[0],
            MapqCount {
                mapq: 0,
                correct: 6,
                kept: 12
            }
        );
        assert_eq!(
            into[2],
            MapqCount {
                mapq: 2,
                correct: 1,
                kept: 1
            }
        );
    }

    #[test]
    fn empty_inputs_produce_an_empty_curve() {
        let truth = Truth::default();
        let records = HashMap::new();
        assert!(counts(&records, &truth, &species_context()).is_empty());
    }
}
