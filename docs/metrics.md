# Metrics and Adjustments

This section describes how Benchpro computes TP/FP/FN and how optional adjustments work.

## Base matching (no alternatives)

For each rank, Benchpro builds sets of taxon names from prediction and gold standard profiles. It classifies:

- **TP**: taxon present in both prediction and gold standard
- **FP**: taxon present only in prediction
- **FN**: taxon present only in gold standard

Taxon equality is strict by default:

- Same rank
- Same name
- No alternative names involved (`Taxon::match_exact`)

## Alternative-name matching (`--allow-alternatives`)

Alternative names are represented as `Taxon.alternative_names` in prediction taxa. The gold standard must *not* contain alternative names; if it does, Benchpro panics.

When `--allow-alternatives` is enabled, Benchpro applies the following logic:

1. **Exact matches first** (`match_exact`).
2. **Alternative matches next** (`match_any`):
   - A prediction taxon can match a gold taxon if any of its alternative names matches the gold name at the same rank.

Ambiguity handling:

- If a prediction taxon matches *multiple* gold taxa (same alternative name matches many), it is considered ambiguous.
- Ambiguous matches are **not counted** toward TP abundance or TP counts, but they are still considered “matched” for FP purposes (they do not become FP).
- Only alternative matches that map uniquely to one gold taxon contribute to TP abundance/counts.

This behavior is implemented in `binary_classification_helper` and `binary_classification_alternatives_df` in `src/profile.rs`.

## Adjusted benchmarks (`--adjusted`)

When `--adjusted` is enabled, Benchpro re-labels some FP/FN using phylogeny:

- Only **Species** rank and **AllowAlternatives = false** are considered.
- Each sample must have a `GoldStdTree` path in meta.
- The sample’s taxa are mapped onto that tree and a subtree is extracted.

For each taxon in the sample:

1. Find its closest neighbor on the tree (`closest_neighbor`).
2. If a **FP** taxon has a closest neighbor that is **FN**, and the distance is below a threshold, the FP is converted into **FFP**.
3. The matched FN is converted into **TP**, with prediction abundance taken from the FP.

Threshold:

- The current distance threshold is hard-coded as `0.04`.

Detectability:

- `AvailableSpecies` (if provided) marks detectable taxa.
- Only pairs where the closest neighbor is **not** detectable are eligible for adjustment.
- Detectable is derived from `DetectableTaxon == "True"`; any other value (including `Unknown`) is treated as not detectable.

Adjusted output columns:

- `ClosestNeighbor`
- `ClosestNeighborType`
- `ClosestNeighborDistance`
- `ClosestNeighborAbundance`
- `Adjusted = True`

## Summary metrics

For each `(Rank, ID, AllowAlternatives, Adjusted)` group, Benchpro computes:

- `TP`, `FP`, `FN` counts
- `F1`, `Sensitivity`, `Precision`

Formulas:

- `F1 = 2*TP / (2*TP + FP + FN)`
- `Sensitivity = TP / (TP + FN)`
- `Precision = TP / (TP + FP)`
