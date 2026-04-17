# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**Benchpro** is a Rust benchmarking suite for taxonomic profilers. It evaluates metagenomics profiling tools by comparing predictions against gold standard datasets, computing binary classification metrics (TP/FP/FN) and abundance correlation statistics across taxonomic ranks.

## Commands

```bash
# Build
cargo build --release

# Run tests
cargo test
cargo test --test integration_profile   # specific integration test

# Lint / format
cargo clippy
cargo fmt

# Install locally (via Justfile)
just install-local
```

## Architecture

The tool is structured as a pipeline with five subcommands: `profile`, `strain`, `merge`, `msa`, and `normalize`. The main entry point is [src/main.rs](src/main.rs), which dispatches to [src/benchpro.rs](src/benchpro.rs) (core logic) or the dedicated subcommand modules.

### Core Data Flow (profile/strain)

```
Meta file (.xlsx/.csv/.tsv)
  → meta.rs: parse MetaEntry[]
  → profile_handler.rs: load ProfileWrapper (prediction + gold standard)
  → format.rs: detect format, map columns
  → profile.rs: parse into Profile<Entry<T>>  (T = taxonomy type)
  → benchpro.rs: match taxa, classify TP/FP/FN, aggregate metrics
  → two TSV output files (per-taxon + summary)
```

### Key Abstractions

- **`Entry<T: Taxonomy>`** ([src/profile.rs](src/profile.rs)) — a single taxonomic entry; generic over taxonomy type to allow compile-time dispatch for GTDB, NCBI, ChocoPhlAn, Custom.
- **`TaxonomicRank`** ([src/common.rs](src/common.rs)) — 9-variant enum (Superkingdom → Strain); the central unit of analysis.
- **`ProfileFormat`** ([src/normalize_detect.rs](src/normalize_detect.rs)) — auto-detected format enum: CAMI, MetaPhlAn, Sylph, Bracken, Custom.
- **`BC` enum** ([src/profile.rs](src/profile.rs)) — binary classification: TP/FP/FN/TN (plus FFP/FFN for phylogeny-adjusted mode).

### Module Responsibilities

| Module | Role |
|--------|------|
| `benchpro.rs` | Core benchmarking logic, TP/FP/FN classification, metric aggregation |
| `profile.rs` | Profile data structures, `Entry<T>`, BC enum |
| `profile_handler.rs` | Loads & caches profiles from meta; parses `matrix::path::column` refs |
| `format.rs` | Format detection heuristics, column mapping |
| `normalize.rs` / `normalize_detect.rs` / `normalize_loader.rs` | Format auto-detection and conversion to canonical (CAMI/plain) form |
| `meta.rs` | Parses meta files, validates structure |
| `common.rs` | `TaxonomicRank`, `Lineage<T>`, `Taxonomy` trait + GTDB/NCBI/ChocoPhlAn impls |
| `merge.rs` | Merges multiple profiles into abundance matrices; sample regex extraction |
| `tree_handler.rs` / `tree_adjusted_benchmarks.rs` | Phylogeny-aware matching (alternative names via tree) |
| `utils.rs` | Stats helpers: Pearson/Spearman, F1, Bray-Curtis, L2 |
| `ncbi_taxdump.rs` | Loads lineages from NCBI taxdump files |

## Code Style (from AGENTS.md)

- Use Polars exclusively for tabular data processing — no manual CSV loops.
- Use `log::error!`/`log::info!`/etc., never `println!`.
- Error types use `thiserror`; provide contextual messages.
- All public items require doc comments with examples.
- Max 100 character line length, 4-space indentation, snake_case functions, PascalCase types.
- Prefer algorithmic efficiency and Rayon parallelism for bulk operations.

## Testing

Integration tests live in [tests/](tests/) (`integration_profile.rs`, `integration_matrix.rs`, `integration_normalize.rs`). Test data is under [data/test_data/](data/test_data/). Unit tests are co-located in each module.

## Documentation

MkDocs site source is in [docs/](docs/); deployed to GitHub Pages via `.github/workflows/docs.yml` on push to `main`.
