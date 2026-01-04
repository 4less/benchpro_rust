# Benchpro

Benchpro is a CLI tool for benchmarking taxonomic profiles against gold standards.
It supports two input styles:

- Per-sample profile files (one profile per sample)
- Abundance matrices (one matrix per tool, with multiple samples per file)

Key features:

- Binary classification of taxa into TP/FP/FN across ranks
- Summary metrics (F1, precision, sensitivity)
- Optional alternative-name matching (`--allow-alternatives`)
- Optional phylogeny-based adjustment (`--adjusted`)

See `Usage` for CLI examples and `Input Formats` for the two meta file formats.
