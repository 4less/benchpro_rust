//! CIGAR arithmetic.
//!
//! Two things are computed here and nowhere else: how many bases a CIGAR consumes of the query and
//! of the reference, and — during the reference replay — the edit distance the CIGAR actually
//! implies when the read is laid against real reference bases.

/// Base counts of a CIGAR, split by what each operation consumes.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CigarCounts {
    /// `M`, `=` and `X`: aligned columns, consuming both query and reference.
    pub matches: u64,
    /// `I`: inserted into the query.
    pub ins: u64,
    /// `D` and `N`: deleted from the query. Skips are counted with deletions, as in the Python.
    pub del: u64,
    /// `S`: soft clipped — present in SEQ but not aligned.
    pub soft: u64,
    /// `H`: hard clipped — absent from SEQ entirely.
    pub hard: u64,
}

impl CigarCounts {
    /// Total clipped bases, soft and hard.
    ///
    /// # Returns
    ///
    /// `soft + hard`.
    pub fn clip(self) -> u64 {
        self.soft + self.hard
    }

    /// Read length as the CIGAR sees it.
    ///
    /// Hard clips are included because they were part of the read even though SEQ no longer holds
    /// them — this is the denominator query coverage is measured against, and SEQ may be `*`.
    ///
    /// # Returns
    ///
    /// `matches + ins + clip`.
    pub fn read_len(self) -> u64 {
        self.matches + self.ins + self.clip()
    }

    /// Alignment length: the usual identity denominator.
    ///
    /// # Returns
    ///
    /// `matches + ins + del` — aligned columns with gaps included.
    pub fn aln_len(self) -> u64 {
        self.matches + self.ins + self.del
    }

    /// Bases of SEQ the CIGAR accounts for.
    ///
    /// SAM requires SEQ to be exactly this long. Hard clips are excluded: they are not in SEQ.
    ///
    /// # Returns
    ///
    /// `matches + ins + soft`.
    pub fn query_consumed(self) -> u64 {
        self.matches + self.ins + self.soft
    }

    /// Bases of the reference the alignment spans.
    ///
    /// # Returns
    ///
    /// `matches + del`.
    pub fn reference_span(self) -> u64 {
        self.matches + self.del
    }
}

/// One CIGAR operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CigarOp {
    /// Repeat count.
    pub len: u64,
    /// Operation character, one of `MIDNSHP=X`.
    pub op: u8,
}

/// Iterator over the operations of a CIGAR string.
///
/// Stops at the first byte that cannot belong to a CIGAR, so a `*` (no CIGAR) yields nothing.
pub struct CigarIter<'a> {
    rest: &'a [u8],
}

impl<'a> Iterator for CigarIter<'a> {
    type Item = CigarOp;

    fn next(&mut self) -> Option<CigarOp> {
        let mut len: u64 = 0;
        let mut digits = 0;
        for (i, &byte) in self.rest.iter().enumerate() {
            if byte.is_ascii_digit() {
                len = len.saturating_mul(10).saturating_add((byte - b'0') as u64);
                digits += 1;
            } else if digits > 0
                && matches!(
                    byte,
                    b'M' | b'I' | b'D' | b'N' | b'S' | b'H' | b'P' | b'=' | b'X'
                )
            {
                self.rest = &self.rest[i + 1..];
                return Some(CigarOp { len, op: byte });
            } else {
                // Not a CIGAR: stop rather than silently skipping bytes.
                self.rest = &[];
                return None;
            }
        }
        self.rest = &[];
        None
    }
}

/// Walks the operations of a CIGAR.
///
/// # Arguments
///
/// * `cigar` - The CIGAR field, e.g. `b"10S90M"`. `b"*"` yields no operations.
///
/// # Returns
///
/// An iterator over `(len, op)` pairs.
pub fn ops(cigar: &[u8]) -> CigarIter<'_> {
    CigarIter { rest: cigar }
}

/// Counts what a CIGAR consumes.
///
/// # Arguments
///
/// * `cigar` - The CIGAR field
///
/// # Returns
///
/// Per-category base counts. `P` (padding) consumes nothing and is ignored.
pub fn count(cigar: &[u8]) -> CigarCounts {
    let mut c = CigarCounts::default();
    for CigarOp { len, op } in ops(cigar) {
        match op {
            b'M' | b'=' | b'X' => c.matches += len,
            b'I' => c.ins += len,
            b'D' | b'N' => c.del += len,
            b'S' => c.soft += len,
            b'H' => c.hard += len,
            _ => {}
        }
    }
    c
}

/// Recomputes an alignment's edit distance against the reference bases it claims.
///
/// This is what makes the identity numbers tool independent. protal emits no `NM` tag at all, and a
/// tool's own `NM` is derived from its own CIGAR, so it cannot reveal a CIGAR that describes no
/// real alignment. Walking SEQ against the actual bases settles both.
///
/// Reference that runs out — the CIGAR walks past the end of the contig — counts as mismatching,
/// which is what it is: the record claims bases the reference does not have.
///
/// # Arguments
///
/// * `cigar` - The record's CIGAR
/// * `seq` - The record's SEQ, uppercase or not
/// * `reference` - The reference bases starting at the alignment's POS, uppercase
///
/// # Returns
///
/// The edit distance: mismatches plus inserted plus deleted bases.
pub fn replay(cigar: &[u8], seq: &[u8], reference: &[u8]) -> u64 {
    let mut nm = 0u64;
    let mut qi = 0usize;
    let mut ri = 0usize;

    for CigarOp { len, op } in ops(cigar) {
        let len_usize = len as usize;
        match op {
            b'M' | b'=' | b'X' => {
                let query = &seq[qi.min(seq.len())..(qi + len_usize).min(seq.len())];
                let refr =
                    &reference[ri.min(reference.len())..(ri + len_usize).min(reference.len())];
                let compared = query.len().min(refr.len());
                nm += query
                    .iter()
                    .zip(refr.iter())
                    .filter(|(a, b)| !a.eq_ignore_ascii_case(b))
                    .count() as u64;
                // Whatever the reference could not supply is a mismatch by definition.
                nm += len - compared as u64;
                qi += len_usize;
                ri += len_usize;
            }
            b'I' => {
                nm += len;
                qi += len_usize;
            }
            b'D' | b'N' => {
                nm += len;
                ri += len_usize;
            }
            b'S' => qi += len_usize,
            _ => {}
        }
    }
    nm
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_every_operation_class() {
        let c = count(b"5S10M2I3D4N6H1P7=8X");
        assert_eq!(c.matches, 10 + 7 + 8);
        assert_eq!(c.ins, 2);
        assert_eq!(c.del, 3 + 4); // N counts with D
        assert_eq!(c.soft, 5);
        assert_eq!(c.hard, 6);
    }

    #[test]
    fn derived_lengths_follow_the_sam_definitions() {
        let c = count(b"10S80M10I");
        assert_eq!(c.clip(), 10);
        assert_eq!(c.read_len(), 80 + 10 + 10);
        assert_eq!(c.aln_len(), 80 + 10);
        assert_eq!(c.query_consumed(), 80 + 10 + 10);
        assert_eq!(c.reference_span(), 80);
    }

    #[test]
    fn hard_clips_are_in_read_len_but_not_in_query_consumed() {
        // The distinction is what the malformed check hangs on: SEQ holds soft clips, not hard.
        let c = count(b"10H90M");
        assert_eq!(c.read_len(), 100);
        assert_eq!(c.query_consumed(), 90);
    }

    #[test]
    fn a_star_cigar_yields_nothing() {
        assert_eq!(count(b"*"), CigarCounts::default());
        assert_eq!(ops(b"*").count(), 0);
    }

    #[test]
    fn malformed_cigar_stops_rather_than_skipping_bytes() {
        assert_eq!(ops(b"10M?5I").count(), 1);
        assert_eq!(ops(b"M").count(), 0);
    }

    #[test]
    fn replay_counts_mismatches_against_the_reference() {
        //            ACGT ACGT
        // reference: ACGTAAGT  -> one mismatch at position 5
        assert_eq!(replay(b"8M", b"ACGTACGT", b"ACGTAAGT"), 1);
        assert_eq!(replay(b"8M", b"ACGTACGT", b"ACGTACGT"), 0);
    }

    #[test]
    fn replay_is_case_insensitive() {
        assert_eq!(replay(b"4M", b"acgt", b"ACGT"), 0);
    }

    #[test]
    fn replay_charges_indels_and_ignores_soft_clips() {
        // 2S: skipped in the query. 2I: two inserted bases. 3D: three deleted reference bases.
        assert_eq!(replay(b"2S4M2I4M", b"XXACGTTTACGT", b"ACGTACGT"), 2);
        assert_eq!(replay(b"4M3D4M", b"ACGTACGT", b"ACGTNNNACGT"), 3);
    }

    #[test]
    fn reference_running_out_counts_as_mismatching() {
        // The record claims 8 bases; the contig only has 4. The missing 4 are not free.
        assert_eq!(replay(b"8M", b"ACGTACGT", b"ACGT"), 4);
        assert_eq!(replay(b"8M", b"ACGTACGT", b""), 8);
    }
}
