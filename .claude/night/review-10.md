# Audit: `benchpro align` (branch `alignment`) vs `TASK.md` + Python source of truth

Checkpoint 10. Auditor: `night-reviewer` (agent `a3949a5cdc62ddb34`). Verdict: **STOP**.

Scope covered: all of `src/align/*`, `src/options.rs`, `src/main.rs`, `tests/integration_align.rs`,
the fixture, `TASK.md`, `CLAUDE.md`, cross-read against `align_metrics.py`, `report.py`, `bench.py`,
`clip_geometry.py`. Build, clippy and tests run; two findings verified empirically.

The work is **not** drifting from TASK.md: module structure, CLI, phases and outputs match
§2/§5/§9/§10. The problems are correctness and two false claims in §13.

---

## CRITICAL

### 1. `aln_records`, `aln_no_nm`, `aln_proper_pair` use the wrong denominator — percentages exceed 100%

`report.rs:458-461` (`aln_records` = `s.score.aligned`), `report.rs:489-502` (`aln_no_nm`,
`aln_proper_pair` = `pct(counters.X, s.score.aligned)`).

`score.aligned` is *truth reads the tool placed*; the Python (`align_metrics.py:341,353-354`,
`n = len(alns)`) divides by *all primary alignments emitted*. `base::summarize` already computes
these correctly over `records.len()` and `pool_base` pools them — `summary_frame` throws that away
and recomputes from raw counters against the wrong denominator. The TSV and the stderr table become
two different computations of the same number.

Verified, reachable (any run where the aligner places reads absent from the truth):

```
truth = first 4 mates only; alignment = sloppy.sam (19 primaries, 10 without NM)
TSV:    aln_records 4   aln_no_nm 250.0   aln_proper_pair 475.0   aln_malformed 5.26
stderr: alignments emitted 19   missing NM tag 52.63%   proper pair 100.00%   malformed 5.26%
```

Also falsifies TASK §13.2 ("PAF rows carry no base-level columns"): a PAF row gets `aln_records`,
`aln_no_nm`, `aln_proper_pair`.

Fix: emit the pooled `base` values; null for PAF. Add a fixture where the alignment contains reads
outside the truth — the current fixture has `aligned == records` for every tool, which is why no
test catches this.

### 2. Counters are thread-count dependent — determinism claim is false

`sam.rs:457-478`: the duplicate check returns *before* `no_nm`/`proper_pair` are counted, so a
worker that has not yet seen a mate's primary counts its duplicate's counters, and `merge_from` sums
them unconditionally.

```
31 MB SAM, 200k mates, each duplicated later (2nd copy lacks NM):
  -t 1 -> aln_no_nm 0.0     md5 b008cdf2...
  -t 2 -> aln_no_nm 96.496  md5 9b88fec8...
  -t 8 -> aln_no_nm 96.496  md5 9b88fec8...
```

TASK §12 and §13.3 are both untrue. `thread_count_does_not_change_the_result` compares only
`pos0/nm/target` on a 6-line file that never spans two batches.

Fix: derive counters from the winning record set, not arrival order. Test with a file larger than
`BUFFER_SIZE` containing cross-batch duplicates, comparing full output bytes across `-t 1/2/8`.

### 3. Reservoir/records coupling breaks on duplicate keys during merge

`sam.rs:207-224`, `admit`/`strip` at `sam.rs:233-259`. `merge_from` drops `other`'s record when
`self` holds a lower-offset one, then still admits that key from `other.reservoir`. The heap can
hold two entries for one key, or an entry whose record was already stripped.

```
--verify-sample 1000:
  -t 1 -> aln_verified 1000, aln_identity_v 24.875
  -t 2 -> aln_verified  507, aln_identity_v 25.178
```

Fix: make `admit` a no-op for a key already in the heap, and only admit a key whose record still
holds a seq.

---

## WARNING

### 4. `--clip-geometry` loads the entire `.fai` — what TASK §8.3 forbids

`mod.rs:351-362`: `Reference::open(reference, None)`, once per row. The replay path filters
(`mod.rs:316-321`); this one does not, and clip geometry is the marker-DB feature. On the 14.5 M
sequence marker reference that is a multi-GB map, then `lengths()` clones it into a second one.

### 5. `judged` includes unjudgeable records, deflating clip percentages

`clip.rs:155` vs the doc at `clip.rs:47`. `clip_geometry.py:71-72` skips unknown-length contigs and
divides by the remainder. Reachable whenever `--clip-geometry` falls back to the `@SQ` header.

### 6. Replay charges query-underrun as mismatch where Python does not

`cigar.rs:209-219` charges `len - min(query, ref)`; `align_metrics.py:292` charges reference
shortfall only. Affects malformed records, which this module deliberately keeps and replays.
Undocumented divergence.

### 7. `mappable_base` follows TASK.md but the Python has since reversed this decision

`mapq.rs:127-153`. `bench.py:1697-1701` now uses `total`, rejecting `mappable_base` as "a
tool-dependent, moving denominator: adding or removing a contender silently rescales everyone's
recall." The Rust doc carries the superseded argument as current. Spec-vs-source conflict — escalate
rather than silently change.

### 8. `--per-read` rows are not restricted to the common sample set

`mod.rs:85-90` collects before `restrict_to_common_samples`. The per-read table then cannot
reconstruct the aggregates it exists to explain.

### 9. `read_exact` on a ragged-line FASTA aborts the run

`reference.rs:309-315`. Python uses `read(e - s)` and tolerates a short read.

### 10. `records` counter and `pos0` conversion diverge at the edges

`sam.rs:438`, `sam.rs:482` (POS 0 → 0 where Python gives −1). Both theoretical.

---

## NOTE

11. clippy not clean: `sam.rs:832` `let mut with_seed`.
12. `genome_pct` missing from `align_summary.tsv` (TASK §9 lists it); cosmetic, `correct_pct` is the
    same quantity.
13. `reference_pct`/`position_pct` denominators unverified against the Python renderer.
14. h2h tolerance default 100 vs Python's 50 — parity runs must pass `--tolerance 50`.
15. Truth loader silently skips `#`-prefixed lines. Theoretical.
16. Verified sound: `sample_hash` determinism, `Sampled::Ord` totality, `verify`'s sort key,
    `mean`'s fixed-point accumulation (no overflow; bias ~5e-16), first-seen ordering.

---

## Verified correct against the Python

`parse_cigar` incl. `N`→del; identity/coverage definitions; the malformed check; `ReadKey` for both
conventions; PAF 0-based start and mate-from-suffix; `score()` iterating truth; nested full strata
and inclusive symmetric tolerance; `mapq::counts`; `curve_from_counts`; `best_f1` tie-break;
`nm_agree`; the h2h never-a-mix rule; `clip_geometry`'s predicate and the `end + SLACK >= len`
rearrangement (which avoids an underflow Python's form would have in Rust); `.fai` arithmetic;
`MultiGzDecoder`; `restrict_to_common_samples`; pooled aggregation and Bessel-corrected sd. No `u64`
underflow found.

The 6 unit + 2 `integration_normalize` failures are genuinely pre-existing.

## Not examined

- `bioreader`'s `read_sam_state_par_ids` internals (findings #2/#3 confirmed empirically anyway).
- The `.fai` race under genuinely concurrent processes.
- TASK §13 parity numbers — not re-run.
- Polars memory on a tens-of-millions-row `--per-read` run.
