//! Is the clipped part of an alignment hanging off the reference, or sitting inside it?
//!
//! A clipped alignment has two very different explanations, and only one makes the clipping
//! legitimate:
//!
//! * **Dovetail** — the read runs off the start or end of the contig, so the clipped bases have
//!   nowhere to align. The clipping is correct and the read should be accepted.
//! * **Contained** — the alignment sits strictly inside the contig with room on both sides, so the
//!   clipped bases *could* have aligned and did not: divergence, a chimera, or a junction.
//!
//! The distinction matters when a coverage gate is discarding reads: a gate that rejects dovetails
//! is throwing away reads for being near a contig edge.

use std::collections::HashMap;

use super::metrics::pct;
use super::sam::AlnRecord;
use super::truth::ReadKey;

/// Bases of tolerance when calling an alignment "flush" with a contig edge.
pub const EDGE_SLACK: u64 = 5;

/// How an alignment's clipping is explained by its geometry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClipKind {
    /// Nothing was clipped.
    None,
    /// The clipped bases would fall off the contig: legitimate.
    Dovetail,
    /// The clipped bases would have fit inside the contig: suspicious.
    Contained,
    /// The contig length is unknown, so the geometry cannot be judged.
    Unknown,
}

/// Distribution of clip explanations across a contender's alignments.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ClipGeometry {
    /// Alignments whose geometry could be judged (clipped or not).
    pub judged: u64,
    /// Alignments with no clipping at all.
    pub unclipped: u64,
    /// Clipped alignments hanging off a contig edge.
    pub dovetail: u64,
    /// Clipped alignments sitting strictly inside their contig.
    pub contained: u64,
    /// Clipped alignments on a contig of unknown length.
    pub unknown: u64,
    /// Mean clipped bases of the dovetailing alignments.
    pub dovetail_mean_bases: Option<f64>,
    /// Mean clipped bases of the contained alignments.
    pub contained_mean_bases: Option<f64>,
}

impl ClipGeometry {
    /// Share of judged alignments that dovetail off a contig edge.
    ///
    /// # Returns
    ///
    /// `100 * dovetail / judged`.
    pub fn dovetail_pct(&self) -> f64 {
        pct(self.dovetail, self.judged)
    }

    /// Share of judged alignments clipped strictly inside their contig.
    ///
    /// # Returns
    ///
    /// `100 * contained / judged`.
    pub fn contained_pct(&self) -> f64 {
        pct(self.contained, self.judged)
    }

    /// Folds another sample's geometry into this one, for pooled aggregation.
    ///
    /// # Arguments
    ///
    /// * `other` - The geometry to add
    pub fn add(&mut self, other: &ClipGeometry) {
        // Means are re-derived from the totals they summarise, so pooling stays a sum of counts.
        let total_bases = |mean: Option<f64>, n: u64| mean.unwrap_or(0.0) * n as f64;
        let dovetail_bases = total_bases(self.dovetail_mean_bases, self.dovetail)
            + total_bases(other.dovetail_mean_bases, other.dovetail);
        let contained_bases = total_bases(self.contained_mean_bases, self.contained)
            + total_bases(other.contained_mean_bases, other.contained);

        self.judged += other.judged;
        self.unclipped += other.unclipped;
        self.dovetail += other.dovetail;
        self.contained += other.contained;
        self.unknown += other.unknown;

        self.dovetail_mean_bases =
            (self.dovetail > 0).then(|| dovetail_bases / self.dovetail as f64);
        self.contained_mean_bases =
            (self.contained > 0).then(|| contained_bases / self.contained as f64);
    }
}

/// Classifies one alignment's clipping.
///
/// # Arguments
///
/// * `record` - The alignment
/// * `contig_length` - Length of the contig it sits on, or `None` when unknown
///
/// # Returns
///
/// Which explanation the geometry supports.
pub fn classify(record: &AlnRecord, contig_length: Option<u64>) -> ClipKind {
    let (lead, tail) = record.clip_ends;
    if lead == 0 && tail == 0 {
        return ClipKind::None;
    }
    let Some(length) = contig_length else {
        return ClipKind::Unknown;
    };

    let start = record.pos0;
    let end = start + record.counts.reference_span();

    // Would the clipped bases fall outside the contig if laid down end to end...
    let off_left = lead > 0 && start < lead;
    let off_right = tail > 0 && end + tail > length;
    // ...and is the alignment actually flush against that edge?
    let flush_left = start <= EDGE_SLACK;
    let flush_right = end + EDGE_SLACK >= length;

    if (off_left && flush_left) || (off_right && flush_right) {
        ClipKind::Dovetail
    } else {
        ClipKind::Contained
    }
}

/// Classifies every alignment of one contender.
///
/// # Arguments
///
/// * `records` - Primary alignments, keyed by read mate
/// * `contig_lengths` - Contig name to length; contigs absent from it are counted as unknown
///
/// # Returns
///
/// The distribution, or `None` when there is nothing to classify.
pub fn summarize(
    records: &HashMap<ReadKey, AlnRecord>,
    contig_lengths: &HashMap<Box<str>, u64>,
) -> Option<ClipGeometry> {
    if records.is_empty() {
        return None;
    }

    let mut geometry = ClipGeometry {
        judged: records.len() as u64,
        ..Default::default()
    };
    let (mut dovetail_bases, mut contained_bases) = (0u64, 0u64);

    for record in records.values() {
        let clipped = record.clip_ends.0 + record.clip_ends.1;
        match classify(record, contig_lengths.get(&record.target).copied()) {
            ClipKind::None => geometry.unclipped += 1,
            ClipKind::Dovetail => {
                geometry.dovetail += 1;
                dovetail_bases += clipped;
            }
            ClipKind::Contained => {
                geometry.contained += 1;
                contained_bases += clipped;
            }
            ClipKind::Unknown => geometry.unknown += 1,
        }
    }

    geometry.dovetail_mean_bases =
        (geometry.dovetail > 0).then(|| dovetail_bases as f64 / geometry.dovetail as f64);
    geometry.contained_mean_bases =
        (geometry.contained > 0).then(|| contained_bases as f64 / geometry.contained as f64);

    Some(geometry)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::cigar;

    fn record(pos0: u64, cigar_str: &str) -> AlnRecord {
        AlnRecord {
            target: "ctg1".into(),
            pos0,
            mapq: 60,
            nm: Some(0),
            counts: cigar::count(cigar_str.as_bytes()),
            clip_ends: cigar::clip_ends(cigar_str.as_bytes()),
            malformed: false,
            cigar: None,
            seq: None,
            vnm: None,
            offset: 0,
        }
    }

    const LENGTH: u64 = 1000;

    #[test]
    fn an_unclipped_alignment_has_no_geometry_to_explain() {
        assert_eq!(classify(&record(500, "100M"), Some(LENGTH)), ClipKind::None);
    }

    #[test]
    fn a_read_hanging_off_the_contig_start_is_a_dovetail() {
        // Aligned from base 0 with 20 bases clipped: those 20 would sit at -20..0, off the contig.
        assert_eq!(
            classify(&record(0, "20S80M"), Some(LENGTH)),
            ClipKind::Dovetail
        );
    }

    #[test]
    fn a_read_hanging_off_the_contig_end_is_a_dovetail() {
        // Ends exactly at the contig end, with 20 clipped bases that would run past it.
        assert_eq!(
            classify(&record(920, "80M20S"), Some(LENGTH)),
            ClipKind::Dovetail
        );
    }

    #[test]
    fn a_clip_in_the_middle_of_the_contig_is_contained() {
        // Room on both sides: those 20 bases could have aligned and did not.
        assert_eq!(
            classify(&record(500, "20S80M"), Some(LENGTH)),
            ClipKind::Contained
        );
        assert_eq!(
            classify(&record(500, "80M20S"), Some(LENGTH)),
            ClipKind::Contained
        );
    }

    #[test]
    fn being_near_an_edge_is_not_enough_the_clip_must_actually_overhang() {
        // Flush with the start, but only 5 bases clipped and 10 bases of contig before it, so the
        // clip would have fitted.
        assert_eq!(
            classify(&record(10, "5S80M"), Some(LENGTH)),
            ClipKind::Contained
        );
    }

    #[test]
    fn the_edge_slack_tolerates_a_few_bases() {
        // Starts 5 bases in -- within EDGE_SLACK of the edge -- and clips 20, so 15 overhang.
        assert_eq!(
            classify(&record(EDGE_SLACK, "20S80M"), Some(LENGTH)),
            ClipKind::Dovetail
        );
        // One base further in is no longer flush.
        assert_eq!(
            classify(&record(EDGE_SLACK + 1, "20S80M"), Some(LENGTH)),
            ClipKind::Contained
        );
    }

    #[test]
    fn an_unknown_contig_length_is_reported_rather_than_guessed() {
        assert_eq!(classify(&record(0, "20S80M"), None), ClipKind::Unknown);
        // ...but an unclipped record needs no length at all.
        assert_eq!(classify(&record(0, "100M"), None), ClipKind::None);
    }

    #[test]
    fn summarize_counts_each_category_and_its_mean_clip() {
        let lengths: HashMap<Box<str>, u64> = [(Box::from("ctg1"), LENGTH)].into_iter().collect();
        let records: HashMap<ReadKey, AlnRecord> = [
            (
                ReadKey {
                    id: "a".into(),
                    mate: 1,
                },
                record(500, "100M"),
            ),
            (
                ReadKey {
                    id: "b".into(),
                    mate: 1,
                },
                record(0, "20S80M"),
            ),
            (
                ReadKey {
                    id: "c".into(),
                    mate: 1,
                },
                record(920, "80M40S"),
            ),
            (
                ReadKey {
                    id: "d".into(),
                    mate: 1,
                },
                record(500, "30S70M"),
            ),
        ]
        .into_iter()
        .collect();

        let geometry = summarize(&records, &lengths).unwrap();

        assert_eq!(geometry.judged, 4);
        assert_eq!(geometry.unclipped, 1);
        assert_eq!(geometry.dovetail, 2);
        assert_eq!(geometry.contained, 1);
        assert_eq!(geometry.dovetail_mean_bases, Some(30.0)); // (20 + 40) / 2
        assert_eq!(geometry.contained_mean_bases, Some(30.0));
        assert_eq!(geometry.dovetail_pct(), 50.0);
        assert_eq!(geometry.contained_pct(), 25.0);
    }

    #[test]
    fn pooling_re_derives_the_means_from_the_totals() {
        let mut a = ClipGeometry {
            judged: 10,
            unclipped: 8,
            dovetail: 2,
            contained: 0,
            unknown: 0,
            dovetail_mean_bases: Some(10.0),
            contained_mean_bases: None,
        };
        a.add(&ClipGeometry {
            judged: 10,
            unclipped: 4,
            dovetail: 6,
            contained: 0,
            unknown: 0,
            dovetail_mean_bases: Some(20.0),
            contained_mean_bases: None,
        });

        assert_eq!(a.dovetail, 8);
        // (2*10 + 6*20) / 8 = 17.5, not the mean of 10 and 20.
        assert_eq!(a.dovetail_mean_bases, Some(17.5));
        assert_eq!(a.contained_mean_bases, None);
    }

    #[test]
    fn nothing_to_classify_yields_nothing() {
        assert!(summarize(&HashMap::new(), &HashMap::new()).is_none());
    }
}
