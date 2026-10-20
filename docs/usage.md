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
- `--taxonomy-index <DIR>`: `taxdet` index folder used to confirm each profile's taxonomy and
  find its release. Default: `~/.taxdet/index` if it exists. See
  [Taxonomy detection](taxonomy-detection.md).
- `--no-taxonomy-index`: Do not use an index, even if `~/.taxdet/index` exists.
- `--allow-taxonomy-mismatch`: Continue when a prediction and its gold standard look like
  different taxonomies (normally an error).
- `--log-level <level>`: `debug`, `info` (default), `warning`, `error`. Put it before the
  subcommand or after it: `benchpro --log-level error profile ...`.

The flags are switches: write `--adjusted`, not `--adjusted true`.

## Outputs

Benchpro writes three files based on `--outprefix`:

- `<outprefix>.tsv`: Summary metrics grouped by rank, sample, and flags. The `TaxonomySource`
  column says whether each row's `Taxonomy` was `given` or `inferred`.
- `<outprefix>_detailed.tsv`: Per-taxon binary classification details.
- `<outprefix>_inputs.tsv`: One row per meta row describing what Benchpro used:
    - the final `Taxonomy` and whether it was given or inferred;
    - for the prediction and the gold standard: inferred taxonomy, confidence, best-matching
      release (with an index) and the evidence behind it;
    - whether each profile loaded and how many species it has;
    - `GoldStdTree` and `AvailableSpecies` (or `NA`);
    - whether `--adjusted` ran, ran leniently, or skipped the row;
    - all warnings for the row.

The run ends with a short report: rows evaluated, taxonomies given or inferred, rows skipped by
`--adjusted`, the number of warnings, and the paths of all three files.

The summary file includes TP/FP/FN counts and derived metrics:

- `F1 = 2*TP / (2*TP + FP + FN)`
- `Sensitivity = TP / (TP + FN)`
- `Precision = TP / (TP + FP)`

See [Metrics](metrics.md) for the abundance metrics and adjustments.

## Errors and warnings

Problems with the meta file or output folder are reported as a single line naming the file,
column or row, and Benchpro exits with a non-zero code (2 for usage errors, 1 for input errors).

Benchpro stops (exit 1) when:

- a row has no `Taxonomy` and its taxonomy cannot be inferred with confidence;
- a prediction and its gold standard look like different taxonomies (GTDB vs NCBI), since
  their names cannot match; override with `--allow-taxonomy-mismatch`;
- rows of one `Tool` use different taxonomies, or different `Taxonomy` values (e.g.
  `GTDB r207` and `GTDB r214`). One tool with one database has one taxonomy; give different
  databases different `Tool` names, such as `MetaPhlAn4` and `MetaPhlAn4_GTDB`. The family mix
  can be overridden with `--allow-taxonomy-mismatch`;
- no row could be evaluated.

Benchpro warns, and keeps going, when:

- a given `Taxonomy` contradicts a file's contents (the given value is used);
- a taxonomy was inferred with only medium confidence;
- with an index, a prediction uses names missing from the gold standard's release;
- with an index, the rows of one tool share no best-matching release;
- a row's prediction had too little evidence and took the taxonomy of its tool's other rows;
- `--adjusted` skips rows without `GoldStdTree` or runs without `AvailableSpecies`;
- a row has 0 true positives at species level (usually a naming mismatch);
- a row could not be evaluated because a profile failed to load.

Repeated warnings are grouped into one line per kind; per-row details are in the `Warnings`
column of `<outprefix>_inputs.tsv`.
