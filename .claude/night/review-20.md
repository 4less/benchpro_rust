# Audit: `benchpro align`, `44fdcd3..HEAD` (checkpoint 20)

Auditor: `night-reviewer` (agent `a6f57840733b57546`). Verdict: **STOP**.

No drift from TASK.md. `score()` is provably unchanged: the auditor built `44fdcd3` in a worktree
and diffed all four output tables against HEAD on the full fixture set and on a `species` run — all
byte-identical. Every finding is in the *new* `align_genomes.tsv`.

---

## CRITICAL

### C1. `align_genomes.tsv`'s `correct_pct` uses a different denominator from the same column elsewhere

`report.rs:1104-1106` renders `correct/total`; `align_summary.tsv` and `align_samples.tsv` render
`correct/aligned` ("of ALIGNED, not of total"). `align_pct` and `position_pct` *are* consistent
between the tables, so nothing signals the switch.

```
tool    align_summary correct_pct   align_genomes correct_pct (genomeA)
good    100.0  (18/18)              90.0  (18/20)
sloppy   68.42 (13/19)              65.0  (13/20)
```

The ranking flips. A conservative tool placing 10% of reads, all right, reads 100.0 in one table and
10.0 in the other. `correct/total` is exactly the summary's `recall_pct`.

### C2. Under `species` scoring, `position_pct` is a hard `0.0` for every genome

`report.rs:1107-1112` emits it unconditionally, but `species` mode can only return `Genome`/`Wrong`,
so `position` is always 0. `summary_frame` gets this right via `has_strata`. A perfect tool is
reported as landing zero reads near the true locus, in the marker-DB mode the module exists for.

### C3. `exact_pct` is fabricated on non-gold rows and vanishes on gold rows that scored zero

`report.rs:1072` decides column presence from a *value* (`any(|r| r.3.exact > 0)`), not from
measurability. (a) With mixed truths in one meta, a TSV-truth row renders `0.0` where the summary
correctly emits an empty cell — the same defect review-15 C3 raised for PAF. (b) With a gold truth
that no alignment reproduces, the column disappears entirely, so the output schema depends on how
good the tool was.

---

## WARNING

- **W1** The new test cannot see C1 or C2: `two_genomes.pred.sam` places every read
  (`aligned == total`) and runs in `full` mode, so both denominators coincide and both undefined
  strata are defined. Verbatim the lesson TASK §13 records from checkpoint 10.
  (`tests/integration_align.rs:590-625`)
- **W2** `per_genome` is held per `(dataset, sample, tool)` for the whole run though only the pooled
  form is used. `datasets × samples × tools × genomes × ~64 B`: 5,000 TIIDs over 20 samples × 5
  tools ≈ 32 MB; 100k labels ≈ 640 MB. TASK §13 rejected group-level Rayon on RSS grounds; this
  deserves the same scrutiny. (`mod.rs:86,107,315`, `report.rs:61`)
- **W3** `truth.rs:323-329` still duplicates `label_of` — the very drift `d6b5896` set out to
  remove, and the cause of review-15's C1/C2.
- **W4** `vocabulary_mismatch` samples 2000 arbitrary `HashMap` entries, so the warning is not
  reproducible: with 99.99% of records on one label, sampling misses the rest ~82% of the time and
  the warning appears on some runs and not others for identical input. (`metrics.rs:493-503`)

---

## NOTE

- `metrics.rs:445-447`: `subset.correct` is computed for indel/clipped and never rendered.
- `metrics.rs:391,409-411`: the genome entry is hashed twice per aligned read.
- `docs/align.md:120-135`: the worked example has `aligned == reads`, so it never states which
  denominator `correct_pct` uses.

---

## Verified clean

`score()` byte-identical to `44fdcd3` across both modes and all four tables. `verdict()`'s rewrite
reproduces both old arms exactly. Per-genome counters fire on exactly the global conditions
(`Σ per_genome.total == truth.len()`, none double-counted); `expect("inserted above")` holds on
every path. The worst-first sort is a total order, no NaN, identical across `-t 4`/`-t 8`. Adding
`correct` to `SubsetScore` does not perturb the indel/clipped renderings. `scratch_path` is a pure
extraction and its new test asserts the property directly, answering review-15's note that the
threaded test might pass with the old code. `mapq.rs`'s corrected inequality is right. Review-15 W2
(truth cache key) fixed. 21/21 integration tests pass; clippy clean in `src/align/`.

## Not examined

Parity against the real `SRR8797712` SAMs; `bioreader` internals; actual RSS for W2 (argued
structurally); `clip_geometry`/`base.rs` (unchanged in range); the 8 pre-existing failures.
