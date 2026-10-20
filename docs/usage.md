# Usage

## Run a benchmark

```bash
mkdir -p output
benchpro profile --meta path/to/meta.xlsx --outprefix output/run1
```

From a source checkout without installing, use `cargo run --release -- profile ...` instead.

The output folder must exist; Benchpro checks this before it starts.

## Flags

- `--meta <PATH>`: Meta file (`.xlsx`, `.csv`, `.tsv`). See [Input Formats](input-formats.md).
- `--outprefix <PREFIX>`: Output prefix for TSV outputs.
- `--validate-meta`: Check columns and file paths in the meta file, then exit.
- `--adjusted`: Add phylogeny-adjusted results. Needs `Dataset` and `GoldStdTree` columns.
- `--allow-alternatives`: Add results with alternative-name matching.
- `--ignore-abundance-error`: Warn instead of stopping when abundances cannot be normalized.
- `--cami-ignore-lineage-error`: Accept CAMI rows whose lineage length does not match the rank.
- `--log-level <level>`: `debug`, `info` (default), `warning`, `error`. Put it before the
  subcommand or after it: `benchpro --log-level error profile ...`.

The flags are switches: write `--adjusted`, not `--adjusted true`.

## Outputs

Benchpro writes two files based on `--outprefix`:

- `<outprefix>_detailed.tsv`: Per-taxon binary classification details.
- `<outprefix>.tsv`: Summary metrics grouped by rank, sample, and flags.

The summary file includes TP/FP/FN counts and derived metrics:

- `F1 = 2*TP / (2*TP + FP + FN)`
- `Sensitivity = TP / (TP + FN)`
- `Precision = TP / (TP + FP)`

See [Metrics](metrics.md) for the abundance metrics and adjustments.

## Errors

Problems with the meta file or output folder are reported as a single line naming the file,
column or row, and Benchpro exits with a non-zero code (2 for usage errors, 1 for input errors).
