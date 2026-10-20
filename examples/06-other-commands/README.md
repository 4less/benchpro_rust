# 06: Other commands

```bash
just example-other-commands
```

## `normalize`: identify and convert profiles

```bash
benchpro normalize -i data/profile_examples/metaphlan402.profile --detect
```

The recipe does this for three profiles from `data/profile_examples/`. Each prints the format,
the tool and its version (`X` = unknown):

```
tool	sylph	X
tool	metaphlan	4
tool	motus	4.0.4
```

`-o FILE --output-format cami` converts a profile to CAMI format.

Converting a **Bracken** profile needs NCBI lineages for its taxids. The first time, Benchpro
downloads the NCBI taxdump (about 700 MB on disk) to `~/.benchpro/ncbi_taxdump` and logs that it
does so. Set `BENCHPRO_NCBI_DIR` to use an existing copy of `nodes.dmp` and `names.dmp`.

## `merge`: per-sample profiles into one matrix

```bash
benchpro merge --input data/test_data/profiles/predictions/sylph/Gastrointestinal_{0,1,2}.profile \
    --output sylph_species_matrix.tsv --rank species \
    --sample-regex '.*/(.*?)/Gastrointestinal_(\d+)\.profile'
```

The capture groups build the column names (`sylph_0`, `sylph_1`, `sylph_2`). Add
`--test-sample-regex` to print the names without merging.

```
Taxon                                                                 sylph_0  sylph_1  sylph_2
d__Bacteria;...;g__Bifidobacterium;s__Bifidobacterium animalis        0.0      0.4649   0.0
```

## `strain`: strain-level trees against a gold standard

[`strain-samplesheet.tsv`](strain-samplesheet.tsv) lists two protal trees (partitioned and
unpartitioned MSA) for *Bacteroides ovatus*. Benchpro writes monophyly scores per gold genome
and MSA error rates:

```
output/strain.monophyly.tsv        output/strain.msa_error.tsv
output/strain.monophyly_tips.tsv   output/strain.msa_gene_error.tsv
output/strain.sample_error.tsv
```

The optional `GoldTree` column (gold standard tree) adds a genome proximity output; the test
data does not include gold trees.
