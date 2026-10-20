# 03: GTDB and NCBI tools together

Four tools on the same sample, two per taxonomy, without a `Taxonomy` column:

| Tool | Database | Gold standard |
|---|---|---|
| `sylph` | GTDB r214 | GTDB r214 |
| `MetaPhlAn4_GTDB` | GTDB r207 | GTDB r207 |
| `MetaPhlAn4` | NCBI | NCBI (CAMI) |
| `mOTUs3` | NCBI | NCBI (CAMI) |

`Tool` means "tool + database": `MetaPhlAn4` and `MetaPhlAn4_GTDB` are different tools to
Benchpro, and each must use one taxonomy across all its rows (see [04-mistakes](../04-mistakes/),
case c).

The NCBI files need `ProfileColumns`/`GoldStdColumns` (`0|2|4|1`) because auto-detection does
not read these CAMI gold standards correctly yet.

## Run

```bash
just example-gtdb-and-ncbi
```

```
Taxonomy inferred as GTDB for 2 row(s): MetaPhlAn4_GTDB_4, sylph_4
Taxonomy inferred as NCBI for 2 row(s): MetaPhlAn4_4, mOTUs3_4
Rescaled 5 profile(s) from percentages to fractions
Done: 4 rows, taxonomy given for 0, inferred for 4
Warnings: 0
```

```bash
cut -f1,3,6-9,12,13 examples/03-gtdb-and-ncbi/output/run_inputs.tsv | column -t -s$'\t'
```

```
ID                 Tool             Taxonomy  TaxonomySource  ProfileTaxonomy  ProfileConfidence  GoldStdTaxonomy  GoldStdConfidence
sylph_4            sylph            GTDB      inferred        GTDB             high               GTDB             high
MetaPhlAn4_GTDB_4  MetaPhlAn4_GTDB  GTDB      inferred        GTDB             high               GTDB             high
MetaPhlAn4_4       MetaPhlAn4       NCBI      inferred        NCBI             high               NCBI             high
mOTUs3_4           mOTUs3           NCBI      inferred        NCBI             high               NCBI             high
```

With a taxonomy index (`~/.taxdet/index`), Benchpro also warns that the MetaPhlAn4 and mOTUs3
predictions use names that do not exist in the 2018 taxdump of the NCBI gold standard.
