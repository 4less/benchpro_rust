# 05: Abundance matrices

Instead of one profile per sample, each tool has one matrix: taxa as rows, samples as columns.
One meta row then covers many samples. Benchpro pairs prediction and gold standard columns with
regexes whose capture group must match:

| Meta row | ProfileRegex | GoldStdRegex | Pairs, for example |
|---|---|---|---|
| `protal_Gastrointestinal` | `^GAST(\d+)$` | `^gastrointestinal_(\d+)$` | `GAST3` with `gastrointestinal_3` |
| `protal_Oral` | `^ORAL(\d+)$` | `^oral_(\d+)$` | `ORAL3` with `oral_3` |

Each pair becomes its own row, with ID `<meta ID>_<capture>` (e.g. `protal_Gastrointestinal_3`).
Matrix meta files still need a `Taxonomy` column.

## Run

```bash
just example-matrix
```

```
Done: 20 rows, taxonomy given for 20, inferred for 0
Warnings: 0
```

The two meta rows expanded to 20 sample rows (10 gastrointestinal, 10 oral); see
`output/run_inputs.tsv` for the list.

To build such a matrix from per-sample profiles, use `benchpro merge`
(see [06-other-commands](../06-other-commands/)).
