# Benchpro

Benchmarking tool for taxonomic profilers (Rust rewrite of 4less/benchpro). Compares predicted
profiles against gold standards and reports TP/FP/FN, F1, precision, sensitivity and abundance
similarity per sample and rank.

Full documentation: `docs/` (MkDocs site).

## Requirements

Rust **nightly** (Benchpro uses unstable features). The included `rust-toolchain.toml` selects
nightly automatically when you build with `rustup`.

## Install

```bash
git clone git@github.com:4less/benchpro_rust.git
cd benchpro_rust
cargo install --locked --path . --root ~/.local   # binary: ~/.local/bin/benchpro
```

## Quick start with the test data

```bash
mkdir -p data/test_data/output
benchpro profile --meta data/test_data/meta/meta.xlsx --outprefix data/test_data/output/test
```

This writes:

- `data/test_data/output/test.tsv`: summary metrics per sample and rank
- `data/test_data/output/test_detailed.tsv`: one row per taxon (TP/FP/FN)

Check a meta file without running the benchmark:

```bash
benchpro profile --meta my_meta.tsv --validate-meta
```

Other subcommands: `strain`, `merge`, `msa`, `normalize`. Run `benchpro <subcommand> --help`.
