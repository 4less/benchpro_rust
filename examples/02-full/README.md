# 02: All optional columns and benchmark modes

Every optional meta column is filled in, and both optional modes are on:

- `--adjusted`: an FP whose closest relative on the gold standard tree is an FN the tool
  cannot detect (not in `AvailableSpecies`) counts as a near miss (`FFP`), and the FN as TP.
- `--allow-alternatives`: predictions with alternative names (e.g. mOTUs' `(A/B)` names) can
  match any of them.

| Column | Value here |
|---|---|
| `Taxonomy` | `GTDB r207` / `GTDB r214` (given, so not inferred) |
| `Dataset` | `Gastrointestinal` |
| `ProfileColumns`, `GoldStdColumns` | `X|0|1` and `0|1|2`: column indices of name, lineage, abundance (see [Input Formats](../../docs/input-formats.md)) |
| `GoldStdTree` | GTDB species tree of the matching release |
| `AvailableSpecies` | Species in the tool's database |

## Run

```bash
just example-full
```

```
Rescaled 2 profile(s) from percentages to fractions
Done: 2 rows, taxonomy given for 2, inferred for 0
Adjusted: 2 rows, skipped 0 (no GoldStdTree)
Warnings: 0
```

## Look at the results

Species level from `output/run.tsv`; each row appears once per mode:

```
ID                 AllowAlternatives  Adjusted  TP  FP  FN  F1
MetaPhlAn4_GTDB_1  false              False     36  3   3   0.923
MetaPhlAn4_GTDB_1  false              True      39  0   0   1.0
MetaPhlAn4_GTDB_1  true               False     36  3   3   0.923
sylph_4            false              False     45  2   0   0.978
sylph_4            false              True      45  2   0   0.978
sylph_4            true               False     45  2   0   0.978
```

The three near misses behind `Adjusted = True`, from `output/run_detailed.tsv`:

```
Name                           Type  ClosestNeighbor                 ClosestNeighborType  ClosestNeighborDistance
s__Achromobacter xylosoxidans  FFP   s__Achromobacter ruhlandii      FN                   0.0077
s__Clostridium_AM drakei       FFP   s__Clostridium_AM scatologenes  FN                   0.0133
s__Clostridium_F sporogenes    FFP   s__Clostridium_F botulinum_A    FN                   0.0152
```

sylph has no FNs, so adjustment changes nothing for it.
