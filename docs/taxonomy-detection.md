# Taxonomy detection

Benchpro needs to know whether a profile uses GTDB or NCBI names, because the two are parsed
differently and their names cannot be compared with each other. If the meta file has no
`Taxonomy` value for a row, Benchpro infers it from the prediction and gold standard files. If
it has one, Benchpro still checks the files and warns when they disagree.

## Text evidence (always used)

Benchpro counts lines with features that only one taxonomy uses:

| GTDB | NCBI |
|---|---|
| `_A`-style suffixes (`g__Prevotella_A`, `p__Firmicutes_A`) | numeric TAXPATHs (`2\|1239\|186801`) |
| placeholder names (`g__UBA1234`, `CAG-110`, `sp000432435`) | bracketed names (`[Clostridium] scindens`) |
| GTDB accessions (`RS_GCF_000145275.1`) | `sp.`, `uncultured`, `Candidatus` |
| `d__Bacteria` domain | `k__Bacteria` kingdom |

Phylum names are not used: NCBI adopted the GTDB-style `-ota` names (for example `Bacillota`)
in 2021, so they occur in both.

Confidence:

- **high**: at least 3 lines of evidence for one taxonomy and none for the other;
- **medium**: one taxonomy has at least 10 times the evidence of the other;
- **low** or **none**: too little or mixed evidence.

A taxonomy is only inferred with high or medium confidence. Otherwise Benchpro stops and asks
for a `Taxonomy` value for that row, showing the evidence it found.

## One taxonomy per tool

A tool run with one database writes every profile in that database's taxonomy. Benchpro
therefore treats the `Tool` column as "tool + database" and checks all its rows together:

- rows of one tool that use different taxonomies (given or inferred) are an **error**, as are
  different `Taxonomy` values within one tool;
- a row whose prediction has too little evidence takes the taxonomy of its tool's other rows,
  with a warning naming the tool;
- with an index, Benchpro warns when the rows of one tool share no best-matching release. This
  is only a warning because one sample's best release is estimated from few names.

Use different `Tool` names for different databases, e.g. `MetaPhlAn4` and `MetaPhlAn4_GTDB`.

## Index evidence (optional)

With a [taxdet](https://github.com/4less/ncbi-taxdump-analyzer) index, Benchpro also looks up
every name of each profile in all GTDB releases and NCBI taxdumps. This:

- confirms the taxonomy (the index whose releases contain at least 80% of the names, and at
  least 10 percentage points more than the other index, wins);
- identifies the release (for example `GTDB r214`); an inferred `Taxonomy` then includes it
  when prediction and gold standard agree;
- finds release mismatches: if at least 1% of a prediction's names (and at least 2) do not
  exist in the gold standard's release, Benchpro warns. Those taxa can only ever count as
  FP/FN.

Set up the index once, then Benchpro finds it in `~/.taxdet/index` automatically:

```bash
git clone https://github.com/4less/ncbi-taxdump-analyzer
mkdir -p ~/.taxdet/index
tar -xzf ncbi-taxdump-analyzer/data/index.tar.gz -C ~/.taxdet/index
```

The unpacked index takes about 310 MB (GTDB about 60 MB, NCBI about 250 MB) and adds well under
a second to a run. Use `--taxonomy-index <DIR>` for another location, or `--no-taxonomy-index`
to turn it off.

## What gets reported

Every decision is written to `<outprefix>_inputs.tsv` with its evidence, for example:

```
ProfileTaxonomy    GTDB
ProfileConfidence  high
ProfileRelease     GTDB r207
ProfileEvidence    text: GTDB features [28 lines with _A-style suffixes, 2 lines with
                   placeholder names, 43 lines with d__ domain]; NCBI features [none];
                   GTDB index: 128/128 names (100%), best GTDB r207 (128 names);
                   NCBI index: 93/128 names (73%), best NCBI 2025-02-01 .. NCBI 2026-03-01
                   (28 releases) (84 names)
Warnings           row 'MG-TK-v2_gastrointestinal_4': prediction best matches GTDB r207,
                   but 7 of its names (5.5%) do not exist in the gold standard's release
                   GTDB r214; those taxa can only count as FP/FN.
```
