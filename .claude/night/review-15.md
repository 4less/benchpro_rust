# Audit: `benchpro align`, commits `79a1602..HEAD` (checkpoint 15)

Auditor: `night-reviewer` (agent `adc44ca87b9399fcb`). Verdict: **STOP**.

Scope: `d7084a1` (gold SAM truth + `exact`), `651dab2` (streaming per-read writer +
`common_samples`), `0cc5cc0` (both MAPQ denominators), `438a4c6` (checkpoint-10 warning fixes),
`874fc4c` (`.fai` temp name, `#`-line rule). Four findings verified empirically.

**No drift from TASK.md.** The gold-SAM/`exact` feature extends beyond §8.1 but is additive,
documented, and TASK.md was updated.

---

## CRITICAL

### C1. A gold-SAM truth whose contig is missing from `Contig2Genome` scores every read wrong

`truth.rs:274-277` falls back to **the contig name**; `metrics.rs:193-197` falls back to the literal
**`"NA"`**. The two must agree and do not. `build_truth.py:44` uses `"NA"`, so the claim that the
label is assigned "exactly as that script assigns it" is false in this case.

```
contig2genome omits 1001_geneA; gold.sam scored against ITSELF, scoring=full
  Truth = gold.sam           -> correct 0   correct_pct 0.0   reference_pct 0.0
  Truth = build_truth.py TSV -> correct 20  correct_pct 100.0 reference_pct 100.0
```

Reachable: an incomplete contig2genome is normal for a reference with plasmids/decoys/phiX. With a
*partially* missing map the numbers are quietly deflated with no message at all.

### C2. `scoring=species` with a gold-SAM truth reports 0% correct, silently

Same root cause. `species` does not require `Contig2Genome` (`meta.rs:364`), so `load_truth_sam`
labels the genome with the full contig name while `verdict()`'s species arm compares
`target.split(sep).next()`. `1001` never equals `1001_geneA`.

```
gold.sam as both Alignment and Truth, scoring=species, no Contig2Genome
  -> total 20  aligned 20  correct 0  correct_pct 0.0  recall_pct 0.0
```

This is the marker-DB configuration the module exists for. `docs/align.md:57` documents behaviour
that does not exist.

### C3. A `.paf` truth fabricates `exact_pct` — 100% or 0%, both invented

`sam.rs:597` gives every PAF record `cigar_fingerprint: 0`; `truth.rs:301` stores `Some(0)`, so a
missing CIGAR is treated as a real one.

```
Truth = mapper.paf, Alignment = mapper.paf -> exact_pct 100.0
Truth = mapper.paf, Alignment = good.sam   -> exact_pct   0.0
```

Exactly what TASK §9 forbids: a rendered zero for an unmeasured quantity.

---

## WARNING

- **W1** `load_truth_sam` drops unmapped and `CIGAR == "*"` gold records; `build_truth.py:41` skips
  only secondary/supplementary. Changes `total` and therefore every "of total" percentage.
  (`truth.rs:269` → `sam.rs:483-496`)
- **W2** `truths` is cached by path only, but its content now depends on `row.contig2genome`. Two
  rows sharing a truth SAM with different maps silently get the first row's labelling.
  (`mod.rs:145-158`)
- **W3** `--per-read` leaves a 0-byte file when no group yields rows — no header, where the other
  tables have one. `polars.read_csv` raises on it. (`report.rs:908-911`)
- **W4** A tab-separated commented-out header (`# read_id\tmate\t...`) now parses as five fields,
  fails the mate check and aborts the run; before `874fc4c` it was skipped. `load_contig2genome`
  instead inserts `"# contig" -> "genome"` into the map. (`truth.rs:146-179`, `360-370`)
- **W5** A partially written `align_reads.tsv` survives a mid-run error, beside no summary at all.
  `ensure_fai` already uses the write-then-rename pattern. (`mod.rs:99-101`)

---

## NOTE

- `mapq.rs:161-166`: the doc claims the base "makes the leading tool peak at 100% recall";
  `report.rs:264` sums per-sample bases, so on a multi-sample dataset no tool reaches 100% unless it
  wins every sample. Documentation, not behaviour.
- `mapq.rs:193-196`: `mappable == 0` with `total > 0` emits `recall_mappable_pct` as `0.0` rather
  than null.
- `metrics.rs:253`: `gold_alignments` uses `.any()`; `.all()` would be safer if mixed truth ever
  becomes reachable.
- `report.rs:223`: `summarize`'s doc references the deleted `restrict_to_common_samples`.
- `report.rs:913`: `ReadsWriter` writes through a bare `File`, not a `BufWriter`.
- `align_mapq.tsv` has no `sample` column although TASK §9 lists one. Pre-existing, cosmetic.

---

## Verified clean

Checkpoint-10 C1 (base denominators), C2/C3 (thread determinism — a 120 MB/500,000-record SAM with
every mate duplicated across batches gave byte-identical output at `-t 1/2/8`), W13 (stratum
denominators now "of total", matching `report.py:669-670`), W6 (query-underrun no longer charged,
matching `align_metrics.py:292`), W5 (`judged` excludes unknown, matching `clip_geometry.py:71-72`),
W4 (`contig_lengths_for` filters and shares), NOTE 11 (clippy clean for `src/align/`).

Also verified: `mapq::add_counts` index alignment; `best_by`'s lower-cutoff tie-break;
`Verdict::Exact = 5` ordering against `is_correct()` and the `>= Position` MAPQ rule; the interning
pool in `load_truth_sam`; `common_samples` equivalence to the old post-hoc restriction; the per-read
header written exactly once; `ReadsWriter`'s CSV settings matching `utils::write_df`.

## Not examined

- `base.rs` beyond the field additions; unchanged in this range.
- The `.fai` cross-process race (code read only). **Note: the added
  `concurrent_index_builds_do_not_collide` test may pass with the old code too, since threads that
  lose the race find the `.fai` already present and skip the build.**
- `clip_geometry.py` parity beyond the `judged` denominator.
- TASK §13 parity numbers against the real `SRR8797712` SAMs.
- `bioreader` internals.
- Polars memory on a tens-of-millions-row `--per-read` run (only 250,000 mates verified).
- The `docs/` changes in `ccf7dbb` beyond claims about the gold SAM.
