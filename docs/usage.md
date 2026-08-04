# Usage

Benchpro is a set of subcommands. Every invocation names one:

| Subcommand | Purpose |
|---|---|
| `profile` | Benchmark taxonomic profiles against a gold standard |
| `strain` | Benchmark at the strain level |
| `merge` | Merge profiles into an abundance matrix |
| `msa` | Compute MSA summary statistics |
| `normalize` | Convert a profile to a canonical format |
| `align` | Benchmark read aligners — see [Alignment Benchmark](align.md) |

## Build and run

```bash
cargo run -- profile --meta path/to/meta.xlsx --outprefix output/path/prefix
```

## Flags

These are the `profile` subcommand's flags; `--log-level` is global. Run
`benchpro <subcommand> --help` for the others.

- `--meta <PATH>`: Meta file (`.xlsx`, `.csv`, `.tsv`).
- `--outprefix <PREFIX>`: Output prefix for TSV outputs.
- `--adjusted <true|false>`: Enable phylogeny-based adjustments. Default: `false`.
- `--allow-alternatives <true|false>`: Enable alternative-name matching. Default: `false`.
- `--log-level <level>`: Log level (`debug`, `info`, `warning`, `error`).

## Outputs

Benchpro writes two main files based on `--outprefix`:

- `<outprefix>_detailed.tsv`: Per-taxon binary classification details.
- `<outprefix>.tsv`: Summary metrics grouped by rank, sample, and flags.

The summary file includes TP/FP/FN counts and derived metrics:

- `F1 = 2*TP / (2*TP + FP + FN)`
- `Sensitivity = TP / (TP + FN)`
- `Precision = TP / (TP + FP)`
