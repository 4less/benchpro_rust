# 04: Common mistakes

Each meta file here contains one mistake. The recipe runs them all and shows the message
and exit code (1 = input problem, 2 = usage problem). None of them gets as far as computing
results.

```bash
just example-mistakes
```

| Case | Meta file | Mistake |
|---|---|---|
| a | `a-missing-column.tsv` | Required column `GoldStd` is missing |
| b | `b-wrong-gold-standard.tsv` | GTDB prediction paired with an NCBI gold standard |
| c | `c-tool-mixes-taxonomies.tsv` | An NCBI profile listed under the GTDB tool `MetaPhlAn4_GTDB` |
| d | `d-empty-cell.tsv` | Empty `Tool` cell |
| e | `e-adjusted-without-tree.tsv` | `--adjusted` without any `GoldStdTree` |
| f | (same as e) | Output folder does not exist |

Expected output (long evidence lines shortened):

```
=== a: required column GoldStd is missing
Failed to load meta file 'examples/04-mistakes/a-missing-column.tsv': Meta file is missing
required column(s): GoldStd. Found columns: ID, Sample, Tool, Profile
--> exit code 1

=== b: GTDB prediction paired with an NCBI gold standard
Taxonomy problems in the meta file:
  row 'sylph_4': prediction looks like GTDB but gold standard looks like NCBI; taxon names
  cannot match. Prediction: text: GTDB features [109 lines with _A-style suffixes, ...];
  NCBI features [none]. Gold standard: text: GTDB features [none]; NCBI features [...]
--> exit code 1

=== c: one tool, two taxonomies (an NCBI file in a GTDB tool)
Taxonomy problems in the meta file:
  tool 'MetaPhlAn4_GTDB' mixes taxonomies: GTDB in 1 row(s) (MetaPhlAn4_GTDB_3); NCBI in 1
  row(s) (MetaPhlAn4_GTDB_4). All profiles of one tool and database must use one taxonomy.
  Check these rows, or give different databases different Tool names (e.g. MetaPhlAn4 and
  MetaPhlAn4_GTDB).
--> exit code 1

=== d: empty required cell (Tool)
Failed to load meta file 'examples/04-mistakes/d-empty-cell.tsv': Meta file row 2: column
'Tool' is empty. This column is required for every row.
--> exit code 1

=== e: --adjusted without any GoldStdTree
--adjusted needs a GoldStdTree for at least one row; none has one
--> exit code 2

=== f: output folder does not exist
Output folder 'examples/04-mistakes/no-such-folder' does not exist. Create it with:
mkdir -p examples/04-mistakes/no-such-folder
--> exit code 2
```

Things that are **not** mistakes: leaving out `Taxonomy`, `GoldStdTree`, `AvailableSpecies` or
`Dataset`, or leaving those cells empty. Benchpro fills them in or skips the feature, and says so
(see [01-minimal](../01-minimal/)).
