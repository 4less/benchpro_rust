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

use std::collections::{HashMap, HashSet};

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
    /// The reported alignment is identical to the gold standard's: same contig, same start, same
    /// CIGAR. Only reachable when the truth came from a gold-standard SAM.
    Exact = 5,
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
            (ScoringMode::Full, Verdict::Exact) => "exact",
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
    /// Reads whose alignment is identical to the gold standard's; only when the truth is a
    /// gold-standard SAM.
    pub exact: Option<u64>,
    /// Recovery of the reads whose gold alignment carries an indel; gold-standard SAM only.
    pub indel: Option<SubsetScore>,
    /// Recovery of the reads whose gold alignment is clipped; gold-standard SAM only.
    pub clipped: Option<SubsetScore>,
}

/// Counts over one subset of the truth — reads sharing a source genome, or a property of their
/// gold alignment.
///
/// A total hides both of the things a subset shows. Overall accuracy is dominated by the easy
/// majority, so whether a tool handles gapped reads is invisible in it; and a tool that is 95%
/// correct may be 100% correct on most genomes and 0% on one, which is a different problem from
/// being uniformly slightly wrong.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SubsetScore {
    /// Truth reads in this subset.
    pub total: u64,
    /// ...that the tool placed anywhere.
    pub aligned: u64,
    /// ...that it placed on the right genome.
    pub correct: u64,
    /// ...that it placed within `tolerance` of the true locus.
    pub position: u64,
    /// ...whose alignment it reproduced exactly.
    pub exact: u64,
}

impl SubsetScore {
    /// Share of the subset the tool placed anywhere.
    ///
    /// # Returns
    ///
    /// `100 * aligned / total`, or 0 for an empty subset.
    pub fn aligned_pct(&self) -> f64 {
        pct(self.aligned, self.total)
    }

    /// Share of the subset placed on the right genome.
    ///
    /// # Returns
    ///
    /// `100 * correct / total`, or 0 for an empty subset.
    pub fn correct_pct(&self) -> f64 {
        pct(self.correct, self.total)
    }

    /// Share of the subset placed within tolerance of the true locus.
    ///
    /// # Returns
    ///
    /// `100 * position / total`, or 0 for an empty subset.
    pub fn position_pct(&self) -> f64 {
        pct(self.position, self.total)
    }

    /// Share of the subset whose alignment was reproduced exactly.
    ///
    /// # Returns
    ///
    /// `100 * exact / total`, or 0 for an empty subset.
    pub fn exact_pct(&self) -> f64 {
        pct(self.exact, self.total)
    }

    /// Folds another sample's counts in, for pooled aggregation.
    ///
    /// # Arguments
    ///
    /// * `other` - The counts to add
    pub fn add(&mut self, other: &SubsetScore) {
        self.total += other.total;
        self.aligned += other.aligned;
        self.correct += other.correct;
        self.position += other.position;
        self.exact += other.exact;
    }
}

/// Sums two optional subset scores, keeping `None` only when neither side has one.
fn add_shape(a: Option<SubsetScore>, b: Option<SubsetScore>) -> Option<SubsetScore> {
    match (a, b) {
        (None, None) => None,
        (a, b) => {
            let mut total = a.unwrap_or_default();
            total.add(&b.unwrap_or_default());
            Some(total)
        }
    }
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
        self.exact = add_option(self.exact, other.exact);
        self.indel = add_shape(self.indel, other.indel);
        self.clipped = add_shape(self.clipped, other.clipped);
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
    /// The genome label this mode derives from an alignment's target.
    ///
    /// # Arguments
    ///
    /// * `target` - The reference sequence an alignment names
    ///
    /// # Returns
    ///
    /// The label that will be compared against the truth's genome.
    pub fn label_of<'a>(&self, target: &'a str) -> &'a str
    where
        Self: 'a,
    {
        match self.scoring {
            ScoringMode::Species => target.split(self.sep).next().unwrap_or(target),
            ScoringMode::Full => self
                .contig2genome
                .and_then(|map| map.get(target))
                .map(|g| &**g)
                .unwrap_or("NA"),
        }
    }

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
        // Both modes ask the same question of a different label; `label_of` is the single place
        // that derives it, so the vocabulary check cannot disagree with the scorer.
        if self.label_of(&record.target) != &*truth.genome {
            return Verdict::Wrong;
        }

        match self.scoring {
            // A marker read may legitimately land on any of its species' markers, so the finer
            // strata are not defined.
            ScoringMode::Species => Verdict::Genome,
            ScoringMode::Full => {
                if *record.target != *truth.contig {
                    return Verdict::Genome;
                }
                if record.pos0.abs_diff(truth.pos0) > self.tolerance {
                    return Verdict::Reference;
                }
                // The third benchmark, and the only one a truth TSV cannot answer: is this the
                // alignment the simulator produced, base for base? Same start and same CIGAR is
                // exactly that -- a shifted start or a different gap placement is a different
                // alignment of the same read to the same locus.
                match (truth.cigar_fingerprint, record.cigar_fingerprint) {
                    (Some(gold), Some(reported))
                        if record.pos0 == truth.pos0 && reported == gold =>
                    {
                        Verdict::Exact
                    }
                    _ => Verdict::Position,
                }
            }
        }
    }
}

/// Per-genome counts, keyed by the truth's genome label.
pub type PerGenome = HashMap<std::sync::Arc<str>, SubsetScore>;

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
    score_detailed(records, truth, context).0
}

/// Scores one contender, and additionally breaks the counts down by source genome.
///
/// # Arguments
///
/// * `records` - Primary alignments, keyed by read mate
/// * `truth` - The sample's per-read truth
/// * `context` - Correctness definition and its parameters
///
/// # Returns
///
/// The mapping score, and one [`SubsetScore`] per genome the truth names.
pub fn score_detailed(
    records: &HashMap<ReadKey, AlnRecord>,
    truth: &Truth,
    context: &ScoringContext,
) -> (MappingScore, PerGenome) {
    let mut score = MappingScore {
        total: truth.len() as u64,
        ..Default::default()
    };
    let stratified = context.scoring == ScoringMode::Full;
    if stratified {
        score.reference = Some(0);
        score.position = Some(0);
    }
    // The exact stratum and the hard-read subsets exist only when the truth can express them.
    let gold_alignments = truth.values().any(|t| t.cigar_fingerprint.is_some());
    if stratified && gold_alignments {
        score.exact = Some(0);
    }
    if gold_alignments {
        score.indel = Some(SubsetScore::default());
        score.clipped = Some(SubsetScore::default());
    }

    let mut per_genome: PerGenome = HashMap::new();

    for (key, entry) in truth {
        // Every subset counts each truth read, placed or not: a read the tool never placed is
        // precisely the kind of failure these are meant to surface.
        let genome = per_genome.entry(entry.genome.clone()).or_default();
        genome.total += 1;
        let shape = entry.shape.unwrap_or_default();
        if shape.indel {
            if let Some(indel) = &mut score.indel {
                indel.total += 1;
            }
        }
        if shape.clipped {
            if let Some(clipped) = &mut score.clipped {
                clipped.total += 1;
            }
        }

        let Some(record) = records.get(key) else {
            continue;
        };
        score.aligned += 1;
        let genome = per_genome
            .get_mut(&entry.genome)
            .expect("inserted above for this entry");
        genome.aligned += 1;
        let verdict = context.verdict(record, entry);
        if verdict.is_correct() {
            score.correct += 1;
            genome.correct += 1;
        }
        if verdict >= Verdict::Position {
            genome.position += 1;
        }
        if verdict >= Verdict::Exact {
            genome.exact += 1;
        }
        if stratified {
            if verdict >= Verdict::Reference {
                score.reference = score.reference.map(|v| v + 1);
            }
            if verdict >= Verdict::Position {
                score.position = score.position.map(|v| v + 1);
            }
            if verdict >= Verdict::Exact {
                score.exact = score.exact.map(|v| v + 1);
            }
        }

        for (in_subset, subset) in [
            (shape.indel, &mut score.indel),
            (shape.clipped, &mut score.clipped),
        ] {
            if !in_subset {
                continue;
            }
            let Some(subset) = subset else { continue };
            subset.aligned += 1;
            if verdict.is_correct() {
                subset.correct += 1;
            }
            if verdict >= Verdict::Position {
                subset.position += 1;
            }
            if verdict >= Verdict::Exact {
                subset.exact += 1;
            }
        }
    }

    (score, per_genome)
}

/// How many records and truth entries to sample when checking label vocabularies.
const VOCABULARY_SAMPLE: usize = 2000;

/// Checks that the truth's genome labels and the labels the scorer derives from predictions are
/// drawn from the same vocabulary.
///
/// A truth written for one scoring mode and scored under another produces a plausible-looking zero
/// rather than an error: `species` compares the contig's prefix (`1001`) while a `full`-style truth
/// names a genome (`genomeA`), and the two never match. The read ids line up, every read is placed,
/// and every read is wrong. That is the failure this catches.
///
/// Only a warning, and only when the two sets are *entirely* disjoint: a tool that places every
/// read on the wrong genome would also produce disjoint sets, and it is not this function's job to
/// call that a misconfiguration.
///
/// # Arguments
///
/// * `records` - Primary alignments
/// * `truth` - The sample's truth
/// * `context` - The correctness definition whose vocabulary is in question
///
/// # Returns
///
/// A message naming an example from each side, or `None` when the vocabularies overlap or either
/// side is empty.
pub fn vocabulary_mismatch(
    records: &HashMap<ReadKey, AlnRecord>,
    truth: &Truth,
    context: &ScoringContext,
) -> Option<String> {
    let predicted: HashSet<&str> = records
        .values()
        .take(VOCABULARY_SAMPLE)
        .map(|record| context.label_of(&record.target))
        .collect();
    let expected: HashSet<&str> = truth
        .values()
        .take(VOCABULARY_SAMPLE)
        .map(|entry| &*entry.genome)
        .collect();

    if predicted.is_empty() || expected.is_empty() || !predicted.is_disjoint(&expected) {
        return None;
    }

    let example = |set: &HashSet<&str>| {
        let mut names: Vec<&str> = set.iter().copied().collect();
        names.sort_unstable();
        names.truncate(3);
        names.join(", ")
    };
    Some(format!(
        "no genome label the alignments imply appears in the truth at all (alignments give {}; \
         truth has {}). Every read will be scored wrong. Check that Scoring matches how the truth \
         was written: 'species' compares the target's prefix before --sep, 'full' maps it through \
         Contig2Genome",
        example(&predicted),
        example(&expected)
    ))
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
            cigar_fingerprint: None,
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
            cigar_fingerprint: None,
            shape: None,
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
            exact: None,
            indel: None,
            clipped: None,
        };
        a.add(&MappingScore {
            total: 100,
            aligned: 100,
            correct: 75,
            reference: Some(70),
            position: Some(60),
            exact: None,
            indel: None,
            clipped: None,
        });

        assert_eq!(a.total, 200);
        assert_eq!(a.aligned, 150);
        assert_eq!(a.correct, 100);
        assert_eq!(a.reference, Some(90));
        // Pooled, not the mean of 50% and 75%.
        assert!((a.correct_pct() - 100.0 * 100.0 / 150.0).abs() < 1e-9);
    }

    #[test]
    fn a_truth_written_for_the_other_mode_is_reported_not_scored_as_zero() {
        // The trap: read ids line up, every read is placed, and every read is wrong -- because the
        // truth names genomes while species scoring compares contig prefixes.
        let truth: Truth = [(key("r1"), truth_entry("src", 0, "genomeA"))]
            .into_iter()
            .collect();
        let records = [(key("r1"), record("1001_geneA", 0))].into_iter().collect();
        let context = ScoringContext {
            scoring: ScoringMode::Species,
            sep: "_",
            contig2genome: None,
            tolerance: 100,
        };

        let problem = vocabulary_mismatch(&records, &truth, &context).expect("mismatch missed");
        assert!(problem.contains("1001"), "{problem}");
        assert!(problem.contains("genomeA"), "{problem}");
        // ...and the score really is the zero the warning is about.
        assert_eq!(score(&records, &truth, &context).correct, 0);
    }

    #[test]
    fn matching_vocabularies_are_not_reported() {
        let truth: Truth = [(key("r1"), truth_entry("src", 0, "1001"))]
            .into_iter()
            .collect();
        let context = ScoringContext {
            scoring: ScoringMode::Species,
            sep: "_",
            contig2genome: None,
            tolerance: 100,
        };
        // Right species, and also a tool that placed this read on the WRONG species: neither is a
        // vocabulary problem, because both labels are drawn from the same namespace.
        for target in ["1001_geneA", "9999_geneA"] {
            let records: HashMap<ReadKey, AlnRecord> = [
                (key("r1"), record(target, 0)),
                (key("r2"), record("1001_geneB", 0)),
            ]
            .into_iter()
            .collect();
            assert!(
                vocabulary_mismatch(&records, &truth, &context).is_none(),
                "false alarm on target {target}"
            );
        }
    }

    #[test]
    fn an_empty_side_is_not_a_mismatch() {
        let context = ScoringContext {
            scoring: ScoringMode::Species,
            sep: "_",
            contig2genome: None,
            tolerance: 100,
        };
        let truth: Truth = [(key("r1"), truth_entry("src", 0, "1001"))]
            .into_iter()
            .collect();
        assert!(vocabulary_mismatch(&HashMap::new(), &truth, &context).is_none());
        assert!(vocabulary_mismatch(
            &[(key("r1"), record("1001_g", 0))].into_iter().collect(),
            &Truth::default(),
            &context
        )
        .is_none());
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
