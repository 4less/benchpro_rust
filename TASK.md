# TASK: `benchpro align` — read-alignment benchmark module

Port the alignment benchmark from `~/git/4less/flexalign_benchmark/` (Python) into benchpro as a
new, self-contained Rust subcommand `align`, using the SAM reader from
[`4less/bioreader`](https://github.com/4less/bioreader).

Branch: `alignment`.

---

## 1. What this is

The taxonomic-profiling side of benchpro (`profile`, `strain`) asks *which taxa are in a sample*.
The alignment benchmark asks a different question about a different input: given a **read aligner's
SAM/PAF output** and **per-read ground truth**, *did each read go to the right place, and is the
alignment it reports any good?*

It is therefore a **separate module** that shares almost nothing with `profile.rs`/`common.rs`
(no `Taxonomy`, no `TaxonomicRank`, no `Lineage`). What it does share is benchpro's shape: a meta
samplesheet in, aggregated TSVs out, Polars for the tabular work, `log::` for output,
`thiserror` for errors.

### Source of truth for the semantics

| Python file | What to port |
|---|---|
| `scripts/align_metrics.py` | Base-level scoring: CIGAR parse, reference replay, identity, head-to-head, metric rows. **The single most important file.** |
| `scripts/report.py` (`parse_sam`, `parse_paf`, `score`, `mapq_counts`, `mappable_base`, `curve_from_counts`) | Mapping-level scoring and the MAPQ precision/recall curve |
| `scripts/bench.py` (`score_and_curve`, `species_score`, `full_score`, `SUMMARY_COLS`, `SAMPLE_COLS`, `restrict_to_common_samples`, `agg`/`pooled`/`weighted`) | Scoring entry points and cross-sample aggregation rules |
| `scripts/build_truth.py`, `build_truth_protal.py` | The truth TSV format this module consumes |
| `scripts/clip_geometry.py` | Optional (phase 6): dovetail-vs-contained clip diagnosis |
| `flexalign_benchmark/CLAUDE.md` §"Alignment (base-level) scoring", §"The protal marker-gene dataset" | The *reasoning* behind the metric choices — read these before implementing |

### Explicitly OUT of scope

Everything in `bench.py` that is **not scoring**: running aligners, `/usr/bin/time` parsing, index
builds, cold/hot cache control, load calibration, thread sweeps, HTML/plotly reports, slurm. That
orchestration stays in `flexalign_benchmark`; benchpro's job is to be the **scorer** it calls
instead of `report.py` + `align_metrics.py`. Keep the module free of any assumption that a tool was
run locally — it only ever sees files on disk.

---

## 2. CLI

```
benchpro align --meta <samplesheet> --outprefix <prefix> [options]
```

Add `Command::Align(AlignArgs)` to `src/options.rs` and dispatch in `src/main.rs`, matching the
existing arms (log the error and `exit(1)` on `Err`).

```rust
/// CLI arguments for the `align` subcommand.
#[derive(ClapArgs, Debug, Clone)]
pub struct AlignArgs {
    /// Samplesheet (.xlsx/.csv/.tsv) with columns:
    /// ID, Sample, Tool, Alignment, Truth (required); Reference, Contig2Genome, Peer, Scoring, Sep (optional).
    #[arg(short = 'm', long, required = true, value_name = "PATH")]
    pub meta: PathBuf,

    /// Output prefix; writes <prefix>.align_summary.tsv, <prefix>.align_samples.tsv,
    /// <prefix>.align_mapq.tsv and (with --per-read) <prefix>.align_reads.tsv.
    #[arg(short = 'o', long, required = true, value_name = "PREFIX")]
    pub outprefix: String,

    /// Scoring mode when the meta file does not set one per row.
    #[arg(long, value_enum, default_value_t = ScoringMode::Full)]
    pub scoring: ScoringMode,

    /// Separator splitting a marker contig into <species-prefix><sep><gene>, for --scoring species.
    #[arg(long, default_value = "_")]
    pub sep: String,

    /// bp of slack for the position stratum and for "same locus" in the head-to-head.
    #[arg(long, default_value_t = 100)]
    pub tolerance: u64,

    /// Replay this many alignments per file against the reference (0 = all).
    #[arg(long, default_value_t = 100_000)]
    pub verify_sample: usize,

    /// Skip the reference replay entirely, even where a Reference column is given.
    #[arg(long, default_value_t = false)]
    pub no_replay: bool,

    /// Worker threads for SAM parsing (0 = all available).
    #[arg(short = 't', long, default_value_t = 0)]
    pub threads: usize,

    /// Also write the per-read verdict table (large: one row per truth read per tool).
    #[arg(long, default_value_t = false)]
    pub per_read: bool,

    /// Seed for reservoir sampling and replay sampling, so runs are reproducible.
    #[arg(long, default_value_t = 0)]
    pub seed: u64,
}

#[derive(ValueEnum, Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScoringMode { Full, Species }
```

---

## 3. Meta samplesheet

Reuse benchpro's loader path (`utils::workbook_to_dataframe` for `.xlsx`, Polars CSV for
`.tsv`/`.csv`) but a **dedicated parser** — do not bend `meta.rs`'s `MetaEntry` (its required
columns are `Profile`/`GoldStd`, which do not apply here). Put it in `align/meta.rs` and mirror
`meta.rs`'s style: a `const REQUIRED_FIELDS`, case-insensitive column matching, a clear
`AlignMetaError` naming the missing column.

| Column | Req | Meaning |
|---|---|---|
| `ID` | ✓ | Row id / dataset label (results grouped by it) |
| `Sample` | ✓ | Sample id within the dataset; aggregation is across samples |
| `Tool` | ✓ | Contender name (the column head in the report) |
| `Alignment` | ✓ | Path to the tool's `.sam`/`.sam.gz`/`.paf`/`.paf.gz` |
| `Truth` | ✓ | Path to the per-read truth TSV (see §4) |
| `Reference` | – | Reference FASTA; present ⇒ replay alignments against it |
| `Contig2Genome` | – | `contig<TAB>genome` TSV; required for `scoring=full`'s genome stratum |
| `Peer` | – | Tool name to run the base-level head-to-head against |
| `Scoring` | – | `full` \| `species`, overriding `--scoring` for this row |
| `Sep` | – | Marker-contig separator, overriding `--sep` |

Validation to perform up front, before any parsing (a long run must not die on row 40):
every `Alignment`/`Truth`/`Reference`/`Contig2Genome` path exists; `Peer` names a `Tool` present
for the same `(ID, Sample)`; `scoring=full` rows have a `Contig2Genome`; `(ID, Sample, Tool)` is
unique. Report **all** problems at once, then bail.

---

## 4. Input formats

**Truth TSV** (headerless, produced by `build_truth.py` / `build_truth_protal.py`):

```
read_id <TAB> mate(1|2) <TAB> true_contig <TAB> true_pos(1-based) <TAB> true_genome
```

For a marker DB, `true_genome` is protal's `TIID` and `true_contig`/`true_pos` refer to the *source
genome*, not the marker — so only the genome field is meaningful there (this is exactly why
`scoring=species` exists).

**Join key is `(read_id, mate)`.** The FASTQ header convention is `<id>/<mate>`; some aligners keep
the `/1`,`/2` suffix and some strip it, so **strip a trailing `/1` or `/2` from QNAME** and take the
mate from the FLAG (`0x40` ⇒ 1, `0x80` ⇒ 2, neither ⇒ 1). PAF has no FLAG, so there the mate comes
from the `/N` suffix and a read without one is skipped. Getting this wrong silently scores zero
overlap — assert non-empty intersection and error out loudly if truth and alignment share no keys.

**Contig2Genome TSV**: `contig<TAB>genome`, headerless.

---

## 5. Module layout

Flat `src/align/` directory (the repo already has `src/taxonomic_profiling/`, so a directory module
is in keeping):

| File | Responsibility |
|---|---|
| `align/mod.rs` | `run_align(&AlignArgs) -> Result<(), AlignError>`; orchestration, per-`(ID, Sample)` grouping, Rayon over rows |
| `align/meta.rs` | Samplesheet parsing + validation (§3) |
| `align/truth.rs` | Truth TSV + contig2genome loading; `ReadKey`, `TruthEntry` |
| `align/sam.rs` | bioreader-backed SAM parsing → `AlnRecord`; PAF fallback parser |
| `align/cigar.rs` | CIGAR op iterator, `(match, ins, del, soft, hard)` counts, malformed check |
| `align/reference.rs` | `.fai` build/load + random-access `fetch`; the replay |
| `align/metrics.rs` | Mapping score, base-level summary, head-to-head |
| `align/mapq.rs` | Cumulative MAPQ counts, mappable base, precision/recall curve |
| `align/report.rs` | DataFrame assembly + TSV writing via `utils::write_df` |
| `align/error.rs` | `AlignError` (`thiserror`) |

Register `pub mod align;` in `src/main.rs`.

---

## 6. Core types

```rust
/// A read mate: the join key between an alignment and its truth.
/// QNAME with any trailing `/1`,`/2` removed, plus the mate number from the FLAG.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ReadKey { pub id: Box<str>, pub mate: u8 }

/// Where a read truly came from.
pub struct TruthEntry { pub contig: Box<str>, pub pos0: u64, pub genome: Box<str> }

/// One primary alignment, reduced to what scoring needs.
/// `cigar`/`seq` are only retained for the sampled subset that gets replayed — they are the bulk
/// of a record and a SAM can hold tens of millions.
pub struct AlnRecord {
    pub target: Box<str>,
    pub pos0: u64,
    pub mapq: u8,
    pub nm: Option<i64>,          // None = tool emits no NM tag (protal)
    pub matches: u64,             // M/=/X
    pub ins: u64,
    pub del: u64,                 // D and N
    pub clip: u64,                // S + H
    pub malformed: bool,          // SEQ length != CIGAR query consumption
    pub cigar: Option<Box<str>>,
    pub seq: Option<Box<str>>,
    pub vnm: Option<u64>,         // reference-replayed edit distance
}
```

Derived quantities (port verbatim from `Aln` in `align_metrics.py`):

```
read_len  = matches + ins + clip          # hard clips included; SEQ may be '*'
aln_len   = matches + ins + del           # identity denominator
coverage  = (matches + ins) / read_len
identity  = max(0, 1 - nm / aln_len)      # None when nm is None or aln_len == 0
has_indel = ins > 0 || del > 0
```

Counters collected during the parse: `records`, `unmapped`, `secondary`, `no_cigar`, `no_nm`,
`proper_pair`, `unknown_ref`.

---

## 7. Parsing with bioreader

Add the dependency:

```toml
bioreader = { git = "https://github.com/4less/bioreader" }
```

(A `path = "../bioreader"` override is fine while iterating; the committed manifest must use the
git form.) bioreader is `edition 2021` on **nightly** (`#![feature(slice_internals)]`) — benchpro is
already nightly (`#![feature(adt_const_params)]` in `main.rs`), so this is compatible. Verify
`cargo build --release` on the pinned toolchain as the very first step of phase 1; if the git
dependency does not resolve, that is a blocker to raise, not to work around silently.

API to use (from `~/git/4less/bioreader`):

- `bioreader::parallel::sam::read_sam_state_par_with_header(file, buffer_size, num_threads, state, f)`
  — the multithreaded driver. `State: Default + Clone + Send + Merge`; each worker gets its own
  state and they are merged at the end. Use a `ParseState { alns: HashMap<ReadKey, AlnRecord>,
  counters: Counters, reservoir: Vec<..> }` and implement
  `bioreader::parallel::fastq::Merge` for it.
- `RefSamRecord` accessors: `qname()`, `flag()`/`is_unmapped()`/`is_secondary()`/
  `is_supplementary()`/`is_primary()`/`is_proper_pair()`/`is_first_in_pair()`/`is_last_in_pair()`,
  `rname()`, `pos()` (**1-based**, subtract 1), `mapq()`, `cigar()`, `seq()`,
  `tag(b"NM").and_then(|t| t.as_int())`.
- The header comes back as `SamHeader`; `reference_sequences()` gives `(name, length)` pairs. Use it
  to warn when an alignment names a target the header does not declare, and as a cheap fallback
  source of contig lengths.
- `SamByteReader::new` takes any `std::io::Read`, so gzip is just
  `flate2::read::MultiGzDecoder` wrapped around the file (`flate2` is already a dependency). Use
  `MultiGzDecoder`, not `GzDecoder` — bgzipped SAM is a multi-member stream and `GzDecoder` stops
  after the first member, truncating the file without erroring.

Rules the parse must follow (all load-bearing, all from `align_metrics.py::parse_sam`):

1. Skip `@` headers, unmapped (`0x4`), secondary (`0x100`), supplementary (`0x800`), `CIGAR == "*"`.
2. **First record for a mate wins** — that is its primary. With threads this is order-dependent, so
   either merge deterministically (keep the record with the lower file offset, via
   `read_sam_state_par_ids` / `SamReadId.byte_offset`) or document that ties are broken arbitrarily.
   Prefer the deterministic route: a benchmark whose numbers wobble between runs is a bug generator.
3. `malformed = seq != "*" && seq.len() != matches + ins + soft` (hard clips are *not* in SEQ).
4. Retain `cigar`/`seq` only for a **reservoir sample** of `--verify-sample` records. Reservoir, not
   the first N: SAM order follows read order, which on simulated data is grouped by source genome,
   so the first N is a handful of species. Seed from `--seed`; with threads, sample per worker with
   a seed derived from the worker index and merge proportionally.

**PAF** (for `flexalign-paf`, `strobealign -x`): mapping-level only, no base-level metrics.
Columns used: `0` QNAME (mate from the `/N` suffix), `5` target, `7` target-start (**already
0-based**), `11` MAPQ. First line for a read is its primary. A small hand-rolled reader is fine —
bioreader has no PAF support.

Dispatch on extension after stripping `.gz`, and **honour the declared format** where the meta says
one: probing blindly is how a stale `<base>.paf` from an earlier contender gets scored as a new
SAM tool's output (`bench.py::find_output` documents exactly this trap).

---

## 8. Metrics to implement

### 8.1 Mapping score

`scoring=species` (marker DB) — correct iff the mapped contig's prefix before `sep` equals the
truth genome label:

```
total   = |truth|
aligned = #{ reads in truth that the tool placed }
correct = #{ aligned where target.split(sep).0 == truth.genome }
```

`scoring=full` (whole-genome reference) — nested strata, each implying the one before:

```
aligned   → placed anywhere
genome    → contig2genome[target] == truth.genome
reference → genome  AND  target == truth.contig
position  → reference AND |pos0 - truth.pos0| <= tolerance
```

Emit per-sample counts and pooled percentages: `align_pct = aligned/total`,
`correct_pct = correct/aligned` (of *aligned*, not of total — that is what `bench.py` reports).

### 8.2 Base-level metrics (needs SAM)

Port `summarize()` exactly. All shares are percentages of emitted primary alignments:

| Key | Definition |
|---|---|
| `aln_records` | primary alignments emitted |
| `aln_identity` / `aln_identity_high` | mean `1-NM/aln_len` over records carrying NM; share ≥ 0.95 |
| `aln_coverage` | mean `(M+I)/read_len` |
| `aln_full_length` | share with `clip == 0` |
| `aln_indel` | share with an indel |
| `aln_malformed` | share where SEQ length ≠ CIGAR query consumption |
| `aln_no_nm` | share missing the NM tag |
| `aln_proper_pair` | share flagged `0x2` |
| `aln_unknown_ref` | share whose target is absent from the reference (a *naming* problem, not a bad alignment — count it separately, never as a mismatch) |
| `aln_verified` | count replayed against the reference |
| `aln_identity_v` / `aln_identity_high_v` | the same identity from the **replayed** edit distance |
| `aln_nm_agree` | share where the tool's NM equals the replayed one |
| `aln_identity_correct` / `aln_identity_wrong` | mean identity split by right-vs-wrong target (replayed identity preferred, falling back to reported) |

`IDENTITY_HIGH = 0.95` as a module constant.

### 8.3 Reference replay — the load-bearing idea

Walk each record's SEQ against the **actual reference bases** through its CIGAR and recount the edit
distance. Necessary because (a) protal emits no NM at all, and (b) a tool's own NM is computed from
its own CIGAR, so it cannot expose a CIGAR that describes no real alignment.

```
M/=/X : mismatches in the zipped chunk, PLUS max(0, n - available_ref)   # ref running out counts
                                                                        # as mismatching: the record
                                                                        # claims bases the reference
                                                                        # does not have
I     : += n, advance query
D/N   : += n, advance reference
S     : advance query only
H/P   : no-op
```

Random access via a `.fai`: read it if present and newer than the FASTA, otherwise **build it**
(plain scan — do not shell out to samtools; benchpro should not acquire a runtime binary
dependency). Write to `<fasta>.fai.tmp` and `rename` so a killed run cannot leave a truncated index.
`.fai` columns: `name`, `length`, `offset`, `linebases`, `linewidth`; the byte offset of base `i` is
`offset + i/linebases*linewidth + i%linebases`.

Filter the loaded index to the contigs the alignments actually name — the marker reference is
14.5 M sequences and loading all of them costs gigabytes.

Fetch in **reference order** (sort the sampled keys by `(target, pos0)`): the marker FASTA is
16.6 GB and sorted seeks keep the reads local. This makes the replay single-threaded per file by
nature; parallelise across *files* instead, and keep the file handle per-thread.

### 8.4 Head-to-head vs `peer`

On a marker DB there is no base-level truth (truth knows the read's source species, not which marker
locus it should land on), so quality is judged against another aligner on the loci **both** found:

```
common      = mates both tools placed
same_locus  = same target AND |Δpos| <= tolerance     (elsewhere they solve different problems)
better/equal/worse, nm_delta = mean(self_nm - peer_nm)   # negative = this tool aligns better
```

Compare **replayed vs replayed** when both were verified, **reported vs reported** otherwise —
never a mix, or the comparison measures bookkeeping rather than alignment. Also report
`h2h_only_self` / `h2h_only_peer`.

### 8.5 MAPQ precision/recall curve

Computed from the same parse that produces the mapping score, so the curve and the headline numbers
cannot disagree — and from the *full* record set, not the sampled subset. Works for PAF contenders
too, so the whole field appears.

Accumulate raw `(mapq, correct, kept)` counts cumulatively from the top MAPQ downward. The recall
denominator is deferred to aggregation, because it depends on every contender:

> **`mappable_base` = the maximum `correct` count at MAPQ 0 across contenders** — not the truth size
> and not the union of what tools placed. On a marker reference only a few percent of whole-genome
> reads can land in a marker at all, so all-reads recall caps near 5%, dominates F1 and pins every
> tool's F1-optimal cutoff at 0. Using the union instead lets a tool that places everything at
> chance accuracy drag the denominator back up to the whole read set. See `report.py::mappable_base`
> for the full argument; carry that comment across.

Then `recall = correct/mappable_base`, `precision = correct/kept`, and report the F1-optimal cutoff
per tool.

### 8.6 Aggregation across samples

Port the rules from `bench.py` — they are decisions, not defaults:

- **Pooled rates, not means of rates**: sum numerators and denominators across samples, then divide.
- **`restrict_to_common_samples`**: intersect the sample set across contenders, drop the rest, and
  put the note in the output. A total over ten samples beside a total over five is nonsense that
  looks fine; means hide the same problem more politely.
- Report sd across samples alongside each aggregate where more than one sample exists.

---

## 9. Outputs

Written with Polars via `utils::write_df`, TSV, one file per concern. **Empty cell for a missing
value — never `0`, which reads as measured.**

`<prefix>.align_summary.tsv` — one row per `(ID, Tool)`:
```
dataset, tool, samples, total, aligned, correct, align_pct, correct_pct, recall_pct,
genome_pct, reference_pct, position_pct,          # scoring=full only
aln_records, aln_verified, aln_identity, aln_identity_high, aln_identity_v, aln_identity_high_v,
aln_nm_agree, aln_identity_correct, aln_identity_wrong, aln_coverage, aln_full_length, aln_indel,
aln_proper_pair, aln_malformed, aln_unknown_ref, aln_no_nm,
h2h_peer, h2h_common, h2h_only_self, h2h_only_peer, h2h_same_locus, h2h_better, h2h_equal,
h2h_worse, h2h_nm_delta,
mapq_best_cutoff, mapq_best_f1, note
```

`<prefix>.align_samples.tsv` — the per-sample rows behind every aggregate (same columns plus
`sample`), so a total can be checked by hand.

`<prefix>.align_mapq.tsv` — `dataset, tool, sample, mapq, correct, kept, recall_pct, precision_pct`.

`<prefix>.align_reads.tsv` (only with `--per-read`) — `dataset, sample, tool, read_id, mate,
truth_contig, truth_pos, truth_genome, aln_target, aln_pos, mapq, verdict, nm, vnm`, where `verdict`
is `unaligned|aligned|genome|reference|position` (`full`) or `unaligned|wrong|correct` (`species`).
This is the file that makes a surprising aggregate debuggable; it is also the biggest, hence opt-in.

Also log a human-readable summary table to stderr in the shape of `align_metrics.py`'s
`METRIC_ROWS`/`H2H_ROWS` (metric label down the left, one column per tool) — that is what the
Python tool prints and what makes the module usable standalone.

---

## 10. Implementation phases

- [x] **Phase 0 — wiring.** Branch `alignment` (done). Add `bioreader` to `Cargo.toml`, `pub mod
      align;`, `Command::Align`, `AlignArgs`, the `main.rs` arm, and an `align/mod.rs` whose
      `run_align` only parses the meta and logs it. `cargo build --release` must pass **before any
      logic is written** — the git dependency on a nightly crate is the one real integration risk.
- [x] **Phase 1 — inputs.** `align/meta.rs`, `align/truth.rs`, `align/error.rs`. Full up-front
      validation with all errors reported at once. Unit tests on small fixtures.
- [x] **Phase 2 — SAM/PAF parsing.** `align/sam.rs` + `align/cigar.rs` on bioreader; counters;
      deterministic primary selection; reservoir sampling; gzip via `MultiGzDecoder`. Unit-test the
      CIGAR arithmetic and the malformed check against hand-written records.
- [x] **Phase 3 — mapping score + MAPQ curve.** `align/metrics.rs`, `align/mapq.rs`, `align/report.rs`.
      At the end of this phase `benchpro align` reproduces `bench.py`'s `align_pct`/`correct_pct`
      and the MAPQ curve for both scoring modes. **Cross-check against the Python output on a real
      SAM before moving on.**
- [x] **Phase 4 — base-level metrics + replay.** `align/reference.rs`, `.fai` build, `summarize`,
      head-to-head. Cross-check `aln_identity_v` and `h2h_nm_delta` against `align_metrics.py`.
- [x] **Phase 5 — aggregation & polish.** Pooled aggregation, `restrict_to_common_samples`, sd,
      `--per-read`, the stderr table, doc comments with examples on every public item, `cargo clippy`
      clean, `cargo fmt`.
- [x] **Phase 6 — clip geometry.** Port `clip_geometry.py`: classify a clipped alignment
      as DOVETAIL (clip hangs off the contig end, so clipping is legitimate) vs CONTAINED (clip sits
      strictly inside, so those bases *could* have aligned and did not). Needs contig lengths, which
      the `.fai` and the SAM header both already provide. Gate behind `--clip-geometry`.

---

## 11. Testing

- **Unit** (co-located `#[cfg(test)]`, per the repo convention): CIGAR counting incl. `=`/`X`/`N`/`P`
  and hard clips; malformed detection; `ReadKey` derivation across the `/1`-kept and `/1`-stripped
  conventions; `.fai` offset arithmetic on a wrapped FASTA (and an unwrapped one); replay against a
  known reference incl. the run-off-the-end case; the cumulative MAPQ counts; pooled aggregation.
- **Integration** — `tests/integration_align.rs`, alongside `integration_profile.rs`. Build a tiny
  fixture under `data/test_data/align/`: a ~200 bp two-contig reference FASTA, a truth TSV of ~20
  mates, and two hand-written SAMs (one "good" tool, one that misplaces reads, lacks NM tags and
  contains a deliberately malformed record) plus one PAF. Assert exact counts on every headline
  metric. Keep it small enough to commit and to read by eye — the fixture doubles as the spec.
- **Parity check** (manual, recorded in the PR description, not committed as a test): run both
  implementations on `~/git/4less/flexalign_benchmark/SRR8797712.flexalign.sam` /
  `.minibwa.sam` and diff the metric tables. Sampled metrics (`aln_identity_v`) will differ within
  sampling noise; the mapping score, the counters and the MAPQ curve must match **exactly**.

## 12. Non-negotiables

- Polars for all tabular output; no manual CSV loops (`AGENTS.md`).
- `log::info!`/`warn!`/`error!` — never `println!` (the stderr metric table goes through `info!`).
- `thiserror` errors with contextual messages naming the offending file and row.
- Doc comments with examples on every public item; 100-column lines; 4-space indent.
- Rayon across files/rows; the replay stays sequential per file for seek locality.
- Determinism: identical inputs and `--seed` must give byte-identical outputs.

## 13. Outcome

All six phases are implemented on branch `alignment`, one commit each.

### Parity with the Python

Checked against `align_metrics.py` / `report.py` on real data from
`~/git/4less/flexalign_benchmark/SRR8797712.flexalign.sam`, not only on the fixture:

| what | records | result |
|---|---|---|
| mapping score + MAPQ curve | 128,222 | exact — `correct_pct` 69.88738282042084, `mappable_base` 89611, 61 cutoffs, best-F1 82.27495374897285 |
| base-level metrics, `--verify-sample 0` | 28,557 | all 16 metrics to 12+ significant digits, incl. `aln_nm_agree` 13.429281787302587 |
| head-to-head vs a perturbed peer | 25,707 common | exact — `same_locus` 88.38059672462754, `better` 61.443661971830984, `nm_delta` -0.9726232394366198 |
| clip geometry | 28,557 | exact — dovetail 0.46223342788%, contained 10.501803411%, mean clips 69.99242424 / 45.30076692 |

Residual differences are 1 ulp, from float summation order (a `HashMap` iterates in a different
order than a Python `dict`).

### Where the time goes, and why there is no Rayon

§12 lists "Rayon across files/rows". Measured on a 96 MB / 285,570-record SAM with a 1.28M-row
truth, the run divides as truth loading ~1.2 s, replay ~0.6 s, parse + score ~0.7 s. Parallelism was
therefore aimed at the largest cost rather than added where the list suggested:

- **Across `(dataset, sample)` groups — rejected.** Each concurrent group holds its own truth and
  record set resident, so it multiplies peak RSS by the concurrency. On a marker-DB run with
  ten-million-read truths that is the wrong trade; the module would get faster and then OOM.
- **Across files within a group — rejected.** The SAM parse is already internally multithreaded
  (bioreader, `--threads`, now defaulting to all cores), so a second layer only oversubscribes.
- **Inside the replay — rejected by design**, as the plan already said: records are fetched in
  reference order because sorted seeks over a 16.6 GB FASTA keep the reads local.
- **Truth loading — addressed directly instead.** Reading the file in one pass and interning the
  contig and genome columns took it from 1.41 s to 1.20 s and peak RSS from 314 MB to 292 MB, with
  byte-identical output. That helps every run rather than trading memory for wall time.

Parallelising the scoring passes (`metrics::score`, `mapq::counts`, `base::summarize`) over the
already-resident records would cost no extra memory and is the one place Rayon still fits; it is
unmeasured and therefore not done.

### Divergences from the Python, deliberate

1. **`recall_pct` uses `correct / total`**, matching `bench.py`'s summary. An earlier draft here
   divided by `mappable_base`; that base belongs only to the MAPQ curve's recall axis. Corrected in
   phase 5.
2. **PAF rows carry no base-level columns at all.** The Python computes them anyway and gets a row
   of zeroes; an empty cell is the honest rendering of "this format cannot answer that".
3. **The primary alignment is chosen by byte offset**, the replay sample is bottom-k on a hash of
   the read key, and `no_nm`/`proper_pair` are derived from the surviving record set — so none of
   them depends on the thread count or on how batches were split. The Python samples in file order,
   which is single threaded and so never had to be order free. Verified byte-identical across runs
   and across `-t 1/2/8` on a 76 MB SAM in which every mate is duplicated.
4. **No samtools dependency.** `.fai` is built by a plain scan; the Python shells out when samtools
   is on `PATH`.
5. **`genome_pct` is not emitted** although §9 lists it: of total it is `recall_pct`, of aligned it
   is `correct_pct`, so it would be a third name for a number already in the row.

### Answers to §13's open questions

1. **Truth format** — the existing headerless 5-column TSV is consumed as-is, so
   `build_truth.py` / `build_truth_protal.py` output works unchanged.
2. **Scope** — scoring only, as planned. `flexalign_benchmark` keeps orchestrating and can call
   `benchpro align` in place of `report.py` + `align_metrics.py`. Timing, index builds, cache
   control and the HTML report are untouched.
3. **Read pairs** — both mates are scored independently, as in the Python. Pair-level metrics
   beyond the `proper_pair` flag are not ported.

### Checkpoint-10 audit (`.claude/night/review-10.md`)

An adversarial review found three CRITICAL defects and ten warnings. All are fixed except one that
needs a decision (below). What it caught, and what it says about the earlier parity claims:

| # | Defect | Why the parity runs missed it |
|---|---|---|
| C1 | `aln_records`/`aln_no_nm`/`aln_proper_pair` divided by truth reads placed, not alignments emitted — percentages over 100% | the fixture and the parity SAM both have `aligned == records`, so the two denominators coincided |
| C2 | `no_nm`/`proper_pair` counted before the duplicate check, so they varied with `--threads` | the parity SAM has no duplicate primaries |
| C3 | a mate seen by two workers could occupy two sample slots or be stripped while still sampled | same |
| W13 | `reference_pct`/`position_pct` divided by `aligned`; `report.py` renders them "of total" | never compared — the parity run used `scoring=species`, which has no strata |

The lesson for any future parity work: **`aligned == records` on every test input hid two wrong
denominators.** Inputs where a tool places reads the truth does not cover are now part of the
fixture set.

### Open question — `mappable_base` is spec-vs-source conflicted

§8.5 of this document specifies `mappable_base` (the most any one contender got right at MAPQ 0) as
the MAPQ curve's recall denominator, quoting `report.py::mappable_base`. The implementation follows
it. But `bench.py:1697-1701` has since **reversed that decision** and uses the truth size, with the
comment that `mappable_base` is

> a tool-dependent, moving denominator: the field's own best result decides the scale every member
> of the field is then measured on. Adding or removing a contender silently rescales everyone's
> recall.

That is a real defect in the current behaviour: adding a contender to a samplesheet moves every
other contender's `mapq_best_cutoff` and `mapq_best_f1`. The code was left spec-conformant rather
than silently switched. **Decide which rule wins**; if the Python's current one does, the change is
`mapq::curve`'s denominator and the stale rationale at `mapq.rs:131-139`.

### Smaller divergences, recorded rather than fixed

- `--tolerance` defaults to 100 and feeds both the position stratum and the head-to-head's "same
  locus"; `align_metrics.py` defaults the latter to 50. A parity comparison must pass
  `--tolerance 50`. §2 asks for one knob, so this is spec-conformant.
- A mapped record with `POS 0` maps to `pos0 = 0`, where the Python gives `-1`. Only reachable for a
  malformed record; the position stratum would differ by one for it.
- The truth loader skips `#`-prefixed lines and rejects a mate that is not `1` or `2`; the Python
  does neither. A read id beginning with `#` is dropped silently.

### Known issues, pre-existing and unrelated

- `cargo test` reports 6 unit-test failures and 2 `integration_normalize` failures on a clean
  checkout, because `data/profile_examples/` is referenced by `include_str!` in
  `normalize_detect.rs` / `normalize_loader.rs` and by `tests/integration_normalize.rs` but was
  never committed. One of the six (`merge::tests::merge_promotes_unclassified_to_species_bucket`)
  fails on its own. None of this involves `align`.
- `Cargo.lock` needed `ethnum` bumped to 1.5.3 before the crate would build at all on the current
  nightly (committed separately, first commit on this branch).
