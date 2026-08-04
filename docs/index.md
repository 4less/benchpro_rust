# Benchpro

Benchpro is a CLI benchmarking suite for metagenomics tools. It answers two different questions,
with two independent halves.

## Taxonomic profiling — `profile`, `strain`, `merge`, `msa`, `normalize`

*Which taxa are in this sample, and did the profiler find them?* Profiles are compared against a
gold standard, taxa classified TP/FP/FN across ranks, and abundance correlations computed.

It supports two input styles:

- Per-sample profile files (one profile per sample)
- Abundance matrices (one matrix per tool, with multiple samples per file)

Key features:

- Binary classification of taxa into TP/FP/FN across ranks
- Summary metrics (F1, precision, sensitivity)
- Optional alternative-name matching (`--allow-alternatives`)
- Optional phylogeny-based adjustment (`--adjusted`)

See [Usage](usage.md) for CLI examples and [Input Formats](input-formats.md) for the two meta file
formats.

## Read alignment — `align`

*Where did each read go, and is the alignment any good?* A read aligner's SAM/PAF output is scored
against per-read ground truth: whether each read reached the right reference, whether its position
is right within a window, and — given a gold-standard SAM — whether the alignment is identical to
the one the simulator produced.

It also reports base-level alignment quality (identity replayed against the reference bases rather
than trusted from the tool's own `NM`), a head-to-head against a peer aligner, accuracy broken down
per source genome, and recovery of the reads that were hard — gapped or clipped.

This half shares none of the taxonomic machinery: no ranks, no lineages. See
[Alignment Benchmark](align.md).

!!! note "Benchpro scores; it does not run tools"
    Neither half runs a profiler or an aligner, builds an index, or measures time. Both take output
    that already exists. Orchestration belongs to whatever harness calls them.
