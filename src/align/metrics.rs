//! Mapping-level scoring: where did each read go, and was that right?
//!
//! Two modes, and the difference is the reference, not a preference:
//!
//! * [`ScoringMode::Full`] — a whole-genome reference, so truth knows the read's exact locus and
//!   correctness is stratified `genome -> reference -> position`, each stratum implying the one
//!   before it.
//! * [`ScoringMode::Species`] — a marker-gene reference, where a read may legitimately land on any
//!   of its species' markers, so the only meaningful question is whether the marker belongs to the
//!   right species. Reference and position strata are not defined and are left out.

use std::collections::HashMap;

use crate::options::ScoringMode;

use super::sam::AlnRecord;
use super::truth::{ReadKey, Truth, TruthEntry};

/// How well a read was placed. Each level implies every level below it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Verdict {
    /// The tool placed no primary alignment for this read.
    Unaligned = 0,
    /// Placed, but on the wrong genome/species.
    Wrong = 1,
    /// Placed on the right genome (or, for `species` scoring, the right species).
    Genome = 2,
    /// Placed on the right contig of the right genome.
    Reference = 3,
    /// Placed within `tolerance` bp of the true locus.
    Position = 4,
}

impl Verdict {
    /// The label written to the per-read table.
    ///
    /// # Arguments
    ///
    /// * `scoring` - Which mode produced the verdict; `species` has no reference/position strata
    ///
    /// # Returns
    ///
    /// A short lowercase label.
    pub fn label(self, scoring: ScoringMode) -> &'static str {
        match (scoring, self) {
            (_, Verdict::Unaligned) => "unaligned",
            (ScoringMode::Species, Verdict::Wrong) => "wrong",
            (ScoringMode::Species, _) => "correct",
            (ScoringMode::Full, Verdict::Wrong) => "aligned",
            (ScoringMode::Full, Verdict::Genome) => "genome",
            (ScoringMode::Full, Verdict::Reference) => "reference",
            (ScoringMode::Full, Verdict::Position) => "position",
        }
    }

    /// Whether this verdict counts as correct in the headline score.
    ///
    /// # Returns
    ///
    /// `true` from [`Verdict::Genome`] upward — landing on the right genome is what "correct"
    /// means in both modes; the finer strata are diagnostics.
    pub fn is_correct(self) -> bool {
        self >= Verdict::Genome
    }
}

/// The mapping score of one contender on one sample.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MappingScore {
    /// Reads in the truth. The denominator of `align_pct`.
    pub total: u64,
    /// Reads the tool placed anywhere.
    pub aligned: u64,
    /// Reads placed on the right genome or species.
    pub correct: u64,
    /// Reads placed on the right contig; `full` scoring only.
    pub reference: Option<u64>,
    /// Reads placed within tolerance of the true locus; `full` scoring only.
    pub position: Option<u64>,
}

impl MappingScore {
    /// Share of the truth the tool placed anywhere.
    ///
    /// # Returns
    ///
    /// `100 * aligned / total`, or 0 for an empty truth.
    pub fn align_pct(&self) -> f64 {
        pct(self.aligned, self.total)
    }

    /// Share of *placed* reads that went to the right target.
    ///
    /// The denominator is `aligned`, not `total` — matching the Python benchmark, where this reads
    /// as "when the tool commits to a placement, how often is it right".
    ///
    /// # Returns
    ///
    /// `100 * correct / aligned`, or 0 when nothing was placed.
    pub fn correct_pct(&self) -> f64 {
        pct(self.correct, self.aligned)
    }

    /// Adds another sample's counts into this one, for pooled aggregation.
    ///
    /// # Arguments
    ///
    /// * `other` - The score to fold in
    pub fn add(&mut self, other: &MappingScore) {
        self.total += other.total;
        self.aligned += other.aligned;
        self.correct += other.correct;
        self.reference = add_option(self.reference, other.reference);
        self.position = add_option(self.position, other.position);
    }
}

/// Sums two optional counters, keeping `None` only when neither side has a value.
fn add_option(a: Option<u64>, b: Option<u64>) -> Option<u64> {
    match (a, b) {
        (None, None) => None,
        (a, b) => Some(a.unwrap_or(0) + b.unwrap_or(0)),
    }
}

/// A percentage, with a zero denominator yielding 0 rather than NaN.
///
/// # Arguments
///
/// * `numerator` - Count of interest
/// * `denominator` - Count it is measured against
///
/// # Returns
///
/// `100 * numerator / denominator`, or 0.
pub fn pct(numerator: u64, denominator: u64) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        100.0 * numerator as f64 / denominator as f64
    }
}

/// Everything needed to judge whether an alignment landed on the right target.
#[derive(Debug, Clone, Copy)]
pub struct ScoringContext<'a> {
    /// Which correctness definition applies.
    pub scoring: ScoringMode,
    /// Marker-contig separator, for `species` scoring.
    pub sep: &'a str,
    /// `contig -> genome`, for `full` scoring.
    pub contig2genome: Option<&'a HashMap<Box<str>, Box<str>>>,
    /// bp of slack allowed on the position stratum.
    pub tolerance: u64,
}

impl ScoringContext<'_> {
    /// Judges one alignment against one read's truth.
    ///
    /// # Arguments
    ///
    /// * `record` - The primary alignment
    /// * `truth` - Where the read really came from
    ///
    /// # Returns
    ///
    /// The strongest verdict the alignment earns.
    pub fn verdict(&self, record: &AlnRecord, truth: &TruthEntry) -> Verdict {
        match self.scoring {
            ScoringMode::Species => {
                // The marker's species is encoded in the contig-name prefix, so no per-contig map
                // is needed: split the mapped contig on `sep`.
                let prefix = record
                    .target
                    .split(self.sep)
                    .next()
                    .unwrap_or(&record.target);
                if prefix == &*truth.genome {
                    Verdict::Genome
                } else {
                    Verdict::Wrong
                }
            }
            ScoringMode::Full => {
                let genome = self
                    .contig2genome
                    .and_then(|map| map.get(&record.target))
                    .map(|g| &**g)
                    .unwrap_or("NA");
                if genome != &*truth.genome {
                    return Verdict::Wrong;
                }
                if *record.target != *truth.contig {
                    return Verdict::Genome;
                }
                if record.pos0.abs_diff(truth.pos0) <= self.tolerance {
                    Verdict::Position
                } else {
                    Verdict::Reference
                }
            }
        }
    }
}

/// Scores one contender's alignments against one sample's truth.
///
/// Iteration is over the *truth*, not over the alignments: a read the tool placed but the truth
/// does not know about is not scoreable, and a read the truth knows but the tool did not place is
/// a miss that must still count against it.
///
/// # Arguments
///
/// * `records` - Primary alignments, keyed by read mate
/// * `truth` - The sample's per-read truth
/// * `context` - Correctness definition and its parameters
///
/// # Returns
///
/// The mapping score. `reference`/`position` are `None` under `species` scoring, where those
/// strata are not defined.
pub fn score(
    records: &HashMap<ReadKey, AlnRecord>,
    truth: &Truth,
    context: &ScoringContext,
) -> MappingScore {
    let mut score = MappingScore {
        total: truth.len() as u64,
        ..Default::default()
    };
    let stratified = context.scoring == ScoringMode::Full;
    if stratified {
        score.reference = Some(0);
        score.position = Some(0);
    }

    for (key, entry) in truth {
        let Some(record) = records.get(key) else {
            continue;
        };
        score.aligned += 1;
        let verdict = context.verdict(record, entry);
        if verdict.is_correct() {
            score.correct += 1;
        }
        if stratified {
            if verdict >= Verdict::Reference {
                score.reference = score.reference.map(|v| v + 1);
            }
            if verdict >= Verdict::Position {
                score.position = score.position.map(|v| v + 1);
            }
        }
    }

    score
}

/// One row of the per-read verdict table.
#[derive(Debug, Clone)]
pub struct ReadVerdict {
    /// The read mate.
    pub key: ReadKey,
    /// Its true origin.
    pub truth: TruthEntry,
    /// Where the tool put it, if anywhere.
    pub target: Option<Box<str>>,
    /// The alignment's 0-based position.
    pub pos0: Option<u64>,
    /// The alignment's MAPQ.
    pub mapq: Option<u8>,
    /// The tool's reported edit distance.
    pub nm: Option<i64>,
    /// The replayed edit distance, when the record was verified.
    pub vnm: Option<u64>,
    /// How well the read was placed.
    pub verdict: Verdict,
}

/// Builds the per-read verdict rows behind a mapping score.
///
/// This is what makes a surprising aggregate debuggable, and it is also the largest output the
/// module produces — hence opt-in behind `--per-read`.
///
/// # Arguments
///
/// * `records` - Primary alignments, keyed by read mate
/// * `truth` - The sample's per-read truth
/// * `context` - Correctness definition and its parameters
///
/// # Returns
///
/// One row per truth read, sorted by read key so output is reproducible.
pub fn read_verdicts(
    records: &HashMap<ReadKey, AlnRecord>,
    truth: &Truth,
    context: &ScoringContext,
) -> Vec<ReadVerdict> {
    let mut rows: Vec<ReadVerdict> = truth
        .iter()
        .map(|(key, entry)| match records.get(key) {
            None => ReadVerdict {
                key: key.clone(),
                truth: entry.clone(),
                target: None,
                pos0: None,
                mapq: None,
                nm: None,
                vnm: None,
                verdict: Verdict::Unaligned,
            },
            Some(record) => ReadVerdict {
                key: key.clone(),
                truth: entry.clone(),
                target: Some(record.target.clone()),
                pos0: Some(record.pos0),
                mapq: Some(record.mapq),
                nm: record.nm,
                vnm: record.vnm,
                verdict: context.verdict(record, entry),
            },
        })
        .collect();
    rows.sort_by(|a, b| a.key.cmp(&b.key));
    rows
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::cigar::CigarCounts;

    fn record(target: &str, pos0: u64) -> AlnRecord {
        AlnRecord {
            target: target.into(),
            pos0,
            mapq: 60,
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

    fn truth_entry(contig: &str, pos0: u64, genome: &str) -> TruthEntry {
        TruthEntry {
            contig: contig.into(),
            pos0,
            genome: genome.into(),
        }
    }

    fn c2g(pairs: &[(&str, &str)]) -> HashMap<Box<str>, Box<str>> {
        pairs
            .iter()
            .map(|(c, g)| ((*c).into(), (*g).into()))
            .collect()
    }

    #[test]
    fn species_scoring_compares_the_contig_prefix() {
        let context = ScoringContext {
            scoring: ScoringMode::Species,
            sep: "_",
            contig2genome: None,
            tolerance: 100,
        };
        // The marker's species is the prefix before the first '_'; the gene id after it is noise.
        assert_eq!(
            context.verdict(&record("1234_567", 0), &truth_entry("src", 0, "1234")),
            Verdict::Genome
        );
        assert_eq!(
            context.verdict(&record("9999_567", 0), &truth_entry("src", 0, "1234")),
            Verdict::Wrong
        );
    }

    #[test]
    fn species_scoring_leaves_the_finer_strata_undefined() {
        let truth: Truth = [(key("r1"), truth_entry("src", 0, "1234"))]
            .into_iter()
            .collect();
        let records = [(key("r1"), record("1234_5", 0))].into_iter().collect();
        let score = score(
            &records,
            &truth,
            &ScoringContext {
                scoring: ScoringMode::Species,
                sep: "_",
                contig2genome: None,
                tolerance: 100,
            },
        );

        assert_eq!(score.correct, 1);
        assert_eq!(score.reference, None);
        assert_eq!(score.position, None);
    }

    #[test]
    fn full_scoring_strata_nest() {
        let map = c2g(&[
            ("ctg1", "genomeA"),
            ("ctg2", "genomeA"),
            ("ctg9", "genomeB"),
        ]);
        let context = ScoringContext {
            scoring: ScoringMode::Full,
            sep: "_",
            contig2genome: Some(&map),
            tolerance: 100,
        };
        let truth = truth_entry("ctg1", 1000, "genomeA");

        assert_eq!(
            context.verdict(&record("ctg1", 1050), &truth),
            Verdict::Position
        );
        assert_eq!(
            context.verdict(&record("ctg1", 5000), &truth),
            Verdict::Reference,
            "right contig, wrong locus"
        );
        assert_eq!(
            context.verdict(&record("ctg2", 1000), &truth),
            Verdict::Genome,
            "right genome, wrong contig"
        );
        assert_eq!(
            context.verdict(&record("ctg9", 1000), &truth),
            Verdict::Wrong
        );
    }

    #[test]
    fn an_unmapped_contig_is_wrong_not_a_crash() {
        let map = c2g(&[("ctg1", "genomeA")]);
        let context = ScoringContext {
            scoring: ScoringMode::Full,
            sep: "_",
            contig2genome: Some(&map),
            tolerance: 100,
        };
        assert_eq!(
            context.verdict(
                &record("unlisted", 1000),
                &truth_entry("ctg1", 1000, "genomeA")
            ),
            Verdict::Wrong
        );
    }

    #[test]
    fn position_tolerance_is_inclusive_and_symmetric() {
        let map = c2g(&[("ctg1", "genomeA")]);
        let context = ScoringContext {
            scoring: ScoringMode::Full,
            sep: "_",
            contig2genome: Some(&map),
            tolerance: 100,
        };
        let truth = truth_entry("ctg1", 1000, "genomeA");

        assert_eq!(
            context.verdict(&record("ctg1", 1100), &truth),
            Verdict::Position
        );
        assert_eq!(
            context.verdict(&record("ctg1", 900), &truth),
            Verdict::Position
        );
        assert_eq!(
            context.verdict(&record("ctg1", 1101), &truth),
            Verdict::Reference
        );
    }

    #[test]
    fn reads_the_tool_never_placed_still_count_against_it() {
        let truth: Truth = [
            (key("r1"), truth_entry("ctg1", 0, "genomeA")),
            (key("r2"), truth_entry("ctg1", 0, "genomeA")),
            (key("r3"), truth_entry("ctg1", 0, "genomeA")),
        ]
        .into_iter()
        .collect();
        let map = c2g(&[("ctg1", "genomeA")]);
        let records = [(key("r1"), record("ctg1", 0))].into_iter().collect();

        let score = score(
            &records,
            &truth,
            &ScoringContext {
                scoring: ScoringMode::Full,
                sep: "_",
                contig2genome: Some(&map),
                tolerance: 100,
            },
        );

        assert_eq!(score.total, 3);
        assert_eq!(score.aligned, 1);
        assert_eq!(score.correct, 1);
        assert!((score.align_pct() - 100.0 / 3.0).abs() < 1e-9);
        assert_eq!(score.correct_pct(), 100.0, "of ALIGNED, not of total");
    }

    #[test]
    fn alignments_absent_from_the_truth_are_not_scoreable() {
        let truth: Truth = [(key("r1"), truth_entry("ctg1", 0, "genomeA"))]
            .into_iter()
            .collect();
        let map = c2g(&[("ctg1", "genomeA")]);
        let records = [
            (key("r1"), record("ctg1", 0)),
            (key("ghost"), record("ctg1", 0)),
        ]
        .into_iter()
        .collect();

        let score = score(
            &records,
            &truth,
            &ScoringContext {
                scoring: ScoringMode::Full,
                sep: "_",
                contig2genome: Some(&map),
                tolerance: 100,
            },
        );
        assert_eq!(score.total, 1);
        assert_eq!(score.aligned, 1);
    }

    #[test]
    fn pooled_addition_sums_counts_not_rates() {
        let mut a = MappingScore {
            total: 100,
            aligned: 50,
            correct: 25,
            reference: Some(20),
            position: Some(10),
        };
        a.add(&MappingScore {
            total: 100,
            aligned: 100,
            correct: 75,
            reference: Some(70),
            position: Some(60),
        });

        assert_eq!(a.total, 200);
        assert_eq!(a.aligned, 150);
        assert_eq!(a.correct, 100);
        assert_eq!(a.reference, Some(90));
        // Pooled, not the mean of 50% and 75%.
        assert!((a.correct_pct() - 100.0 * 100.0 / 150.0).abs() < 1e-9);
    }

    #[test]
    fn verdict_labels_differ_by_mode() {
        assert_eq!(Verdict::Genome.label(ScoringMode::Species), "correct");
        assert_eq!(Verdict::Genome.label(ScoringMode::Full), "genome");
        assert_eq!(Verdict::Wrong.label(ScoringMode::Species), "wrong");
        assert_eq!(Verdict::Unaligned.label(ScoringMode::Full), "unaligned");
    }

    #[test]
    fn per_read_rows_cover_every_truth_read_and_are_sorted() {
        let truth: Truth = [
            (key("r2"), truth_entry("ctg1", 0, "genomeA")),
            (key("r1"), truth_entry("ctg1", 0, "genomeA")),
        ]
        .into_iter()
        .collect();
        let map = c2g(&[("ctg1", "genomeA")]);
        let records = [(key("r1"), record("ctg1", 0))].into_iter().collect();

        let rows = read_verdicts(
            &records,
            &truth,
            &ScoringContext {
                scoring: ScoringMode::Full,
                sep: "_",
                contig2genome: Some(&map),
                tolerance: 100,
            },
        );

        assert_eq!(rows.len(), 2);
        assert_eq!(&*rows[0].key.id, "r1");
        assert_eq!(rows[0].verdict, Verdict::Position);
        assert_eq!(rows[1].verdict, Verdict::Unaligned);
        assert!(rows[1].target.is_none());
    }
}
