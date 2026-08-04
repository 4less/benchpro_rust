# Alignment benchmark (`align`)

`benchpro align` scores **read aligners**: given a tool's SAM/PAF output and per-read ground truth,
it asks where each read went and whether the alignment it reports is any good.

It is independent of the taxonomic-profiling side of Benchpro — no taxonomy, no ranks, no lineages.
What it shares is the shape: a samplesheet in, aggregated TSVs out.

!!! note "It scores files; it does not run tools"
    `align` never invokes an aligner, builds an index, or measures time. Give it output that already
    exists. Orchestration — running the tools, timing them, controlling the page cache — belongs to
    whatever harness calls this.

```bash
benchpro align --meta samplesheet.tsv --outprefix results/run1
```

## Samplesheet

`.xlsx`, `.csv` or `.tsv`. Column names are matched case-insensitively.

| Column | Required | Meaning |
|---|---|---|
| `ID` | ✓ | Dataset label; results are grouped by it |
| `Sample` | ✓ | Sample id; aggregation is across samples |
| `Tool` | ✓ | Contender name |
| `Alignment` | ✓ | `.sam`, `.sam.gz`, `.paf` or `.paf.gz` |
| `Truth` | ✓ | Per-read truth TSV (below) |
| `Reference` | – | Reference FASTA. Present ⇒ alignments are replayed against it |
| `Contig2Genome` | – | `contig<TAB>genome` TSV. Required for `full` scoring |
| `Peer` | – | Another `Tool` on the same sample, for the base-level head-to-head |
| `Scoring` | – | `full` or `species`, overriding `--scoring` for this row |
| `Sep` | – | Marker-contig separator, overriding `--sep` |

The whole file is validated before any scoring starts, and **every** problem is reported at once —
missing files, duplicate `(ID, Sample, Tool)`, a `Peer` that names no tool on that sample, `full`
scoring without a `Contig2Genome`.

```
ID	Sample	Tool	Alignment	Truth	Reference	Peer
gut	s1	flexalign	s1/flexalign.sam	s1/reads.truth.tsv	ref.fna	minibwa
gut	s1	minibwa	s1/minibwa.sam	s1/reads.truth.tsv	ref.fna	
```

## Truth file

Headerless, five tab-separated columns — the format written by the flexalign benchmark's
`build_truth.py`:

```
read_id	mate(1|2)	true_contig	true_pos(1-based)	true_genome
```

Reads join on `(read_id, mate)`. The FASTQ convention is `<id>/<mate>`; some aligners keep that
suffix and some strip it, so the suffix is removed from QNAME and the mate is taken from the FLAG.
If an alignment and its truth share **no** read at all, `align` fails loudly rather than reporting a
flawless run of zeroes — that is almost always a read-id convention mismatch.

On a marker database `true_genome` is the species id and `true_contig`/`true_pos` refer to the
read's *source genome*, not to a marker. Only the genome column is meaningful there, which is what
`species` scoring exists for.

## The two scoring modes

**`--scoring full`** — a whole-genome reference, so truth knows the exact locus. Correctness is
stratified, each level implying the one before:

```
aligned    placed anywhere
genome     contig2genome[target] == truth genome
reference  genome  AND target == truth contig
position   reference AND |pos - truth pos| <= --tolerance
```

**`--scoring species`** — a marker-gene reference, where a read may legitimately land on any of its
species' markers. Correct iff the target's prefix before `--sep` is the truth genome. The reference
and position strata are not defined and are omitted.

## Flags

| Flag | Default | Meaning |
|---|---|---|
| `--meta <PATH>` | – | Samplesheet |
| `--outprefix <PREFIX>` | – | Output prefix |
| `--scoring <full\|species>` | `full` | Default scoring mode |
| `--sep <STR>` | `_` | Marker-contig separator for `species` |
| `--tolerance <INT>` | `100` | bp slack for the position stratum and for "same locus" |
| `--verify-sample <INT>` | `100000` | Alignments replayed against the reference per file; `0` = all |
| `--no-replay` | off | Skip the replay even where a `Reference` is given |
| `--threads <INT>` | `0` | Parse workers; `0` = all cores |
| `--per-read` | off | Also write the per-read verdict table |
| `--clip-geometry` | off | Classify clipped alignments as dovetail or contained |
| `--seed <INT>` | `0` | Seed for the replay sample |
| `--validate-meta` | off | Check the samplesheet and exit |

## Outputs

| File | One row per |
|---|---|
| `<prefix>.align_summary.tsv` | `(dataset, tool)` — the headline table |
| `<prefix>.align_samples.tsv` | `(dataset, sample, tool)` — the rows behind every aggregate |
| `<prefix>.align_mapq.tsv` | `(dataset, tool, mapq)` — the precision/recall curve |
| `<prefix>.align_reads.tsv` | `(dataset, sample, tool, read)` — only with `--per-read` |

A missing value is an **empty cell, never `0`** — a tool that emits no `NM` tag has no reported
identity, and rendering that as zero would read as "every alignment is a complete mismatch". PAF
contenders therefore have every `aln_*` column empty: PAF carries no CIGAR, SEQ or NM.

### Mapping metrics

| Column | Definition |
|---|---|
| `align_pct` | placed / all truth reads |
| `correct_pct` | correct / **placed** — precision. Rewards a tool for placing less |
| `recall_pct` | correct / **all truth reads** — comparable however permissive the tool is |
| `reference_pct`, `position_pct` | strata, as shares of all truth reads (`full` only) |
| `position_precision_pct` | position / placed (`full` only) |

Both `correct_pct` and `recall_pct` are reported because neither alone is honest: a tool that places
its most confident 1% and gets them all right scores 100% precision.

### Base-level metrics (SAM only)

The load-bearing one is the **reference replay**: each record's SEQ is walked against the actual
reference bases through its CIGAR and the edit distance recomputed. That is what makes the identity
numbers tool-independent — some tools emit no `NM` tag at all, and a tool's own `NM` comes from its
own CIGAR, so it cannot reveal a CIGAR that describes no real alignment.

| Column | Definition |
|---|---|
| `aln_identity_v` | mean identity from the **replayed** edit distance |
| `aln_nm_agree` | share where the tool's `NM` matches the replay |
| `aln_identity` | mean identity as the tool reports it |
| `aln_identity_correct` / `aln_identity_wrong` | identity split by right vs wrong target — a healthy aligner's wrong-target hits are its worse ones |
| `aln_coverage`, `aln_full_length`, `aln_indel` | query coverage, unclipped share, indel share |
| `aln_malformed` | CIGAR does not match SEQ length: emitted but uninterpretable |
| `aln_unknown_ref` | target absent from the reference — a naming problem, not a bad alignment |

All are shares of **alignments emitted**, not of truth reads placed. The two differ whenever a tool
places a read the truth does not cover.

Replay needs random access, so `<fasta>.fai` is built on first use (a plain scan — no samtools
dependency) and reused. Only `--verify-sample` alignments per file are replayed; the metrics are
distributions, so a large sample answers them as well as the full set.

### Head-to-head (`Peer`)

On a marker database there is no base-level truth — truth knows the read's source species, not which
marker locus it should land on. So quality is also judged against another aligner on the loci both
found: `h2h_same_locus`, `h2h_better`/`equal`/`worse`, and `h2h_nm_delta` (negative = this tool
reaches the lower edit distance). Replayed distances are compared with replayed, reported with
reported, never a mix — a mixed comparison would measure bookkeeping rather than alignment.

### MAPQ curve

At each cutoff, only alignments with `MAPQ >= q` are kept. An aligner whose MAPQ is informative
trades a little recall for a lot of precision; one whose MAPQ is noise moves along a flat line.

Recall is reported against **two denominators**, because neither is right on its own and the choice
changes which cutoff looks best:

| Column | Divides by | Property |
|---|---|---|
| `recall_mappable_pct`, `f1_mappable` | the most any one contender got right | the only denominator under which the F1-optimal cutoff is meaningful on a marker reference — but **field-relative**, so adding or removing a contender rescales it |
| `recall_total_pct`, `f1_total` | every read in the truth | fixed and comparable across runs with different fields — but on a marker reference it caps near 5%, dominates F1, and pins the optimal cutoff at 0 |

The summary carries the F1-optimal cutoff under each: `mapq_best_cutoff_mappable` /
`mapq_best_f1_mappable` and `mapq_best_cutoff_total` / `mapq_best_f1_total`. They can disagree, and
when they do, the disagreement is the point — read the one whose denominator suits the question you
are asking.

The summary's own `recall_pct` is a separate, mapping-level number and always divides by the truth
size.

### Clip geometry (`--clip-geometry`)

A clipped alignment has two very different explanations:

- **dovetail** — the read runs off a contig edge, so the clipped bases have nowhere to align. The
  clipping is correct.
- **contained** — the alignment sits inside the contig with room on both sides, so those bases
  *could* have aligned and did not.

Contig lengths come from the reference's `.fai` where there is one and from the SAM's `@SQ` header
otherwise. Alignments on a contig of unknown length are counted in `clip_unknown_contig` and left
out of the percentages.

## Reproducibility

Identical inputs and `--seed` give **byte-identical** output, independent of `--threads`:

- the primary alignment of a mate is the record at the lowest byte offset, not whichever worker saw
  it first;
- the replay sample is the records with the lowest hash of `(read id, mate, seed)`, so it does not
  depend on how batches were split;
- means are accumulated in fixed point, so they do not depend on hash-map iteration order.

## Aggregating across samples

Rates are **pooled**, never averaged: numerators and denominators are summed across samples and
divided once, so a 1000-read sample does not weigh the same as a million-read one.

Samples that not every contender produced are **dropped**, and the note appears in the `note` column.
Aggregates are sums, and a sum over ten samples is not comparable with a sum over five.
