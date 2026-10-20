# 01: Minimal meta file

The smallest meta file Benchpro accepts: five columns, no `Taxonomy`.

[`meta.tsv`](meta.tsv) (paths shortened):

```
ID        Sample              Tool    Profile                                GoldStd
sylph_3   gastrointestinal_3  sylph   .../sylph/Gastrointestinal_3.profile   .../GTDB/r214/abundance_3_gtdb.profile
sylph_4   gastrointestinal_4  sylph   .../sylph/Gastrointestinal_4.profile   .../GTDB/r214/abundance_4_gtdb.profile
protal_3  gastrointestinal_3  protal  .../protal/3.profile                   .../GTDB/r214/abundance_3_gtdb.profile
protal_4  gastrointestinal_4  protal  .../protal/4.profile                   .../GTDB/r214/abundance_4_gtdb.profile
```

## Run

```bash
just example-minimal
```

The recipe first checks the meta file, then runs the benchmark:

```
Row 'sylph_3': Taxonomy GTDB (inferred)
...
Meta validation succeeded

Taxonomy inferred as GTDB for 4 row(s): protal_3, protal_4, sylph_3, sylph_4
Rescaled 2 profile(s) from percentages to fractions
Done: 4 rows, taxonomy given for 0, inferred for 4
Warnings: 0
Summary:  examples/01-minimal/output/run.tsv
Per taxon: examples/01-minimal/output/run_detailed.tsv
Inputs (what was given, inferred and skipped): examples/01-minimal/output/run_inputs.tsv
```

## Look at the results

Species-level metrics from `output/run.tsv`:

```
ID        Rank     TP  FP  FN  F1     BrayCurtisSimilarity
protal_3  Species  45  0   0   1.0    0.960
protal_4  Species  45  0   0   1.0    0.962
sylph_3   Species  45  1   0   0.989  0.989
sylph_4   Species  45  2   0   0.978  0.990
```

What Benchpro inferred, from `output/run_inputs.tsv`:

```bash
cut -f1,3,6-9 examples/01-minimal/output/run_inputs.tsv | column -t -s$'\t'
```

```
ID        Tool    Taxonomy  TaxonomySource  ProfileTaxonomy  ProfileConfidence
sylph_3   sylph   GTDB      inferred        GTDB             high
...
```

Column `ProfileEvidence` (11) shows why, for example
`109 lines with _A-style suffixes, 8 lines with placeholder names, ...`.
