# Examples

Small, runnable demos built on the repository's test data. Each folder has the input files and
a README with the expected output; the commands are `just` recipes in the repository's
`Justfile`. Run them from the repository root:

```bash
just --list                  # the [examples] group lists one recipe per demo
just --show example-minimal  # print a demo's commands without running it
just example-minimal         # build Benchpro and run one demo
just examples                # run all demos
```

The recipes build `target/release/benchpro`. To use an installed copy instead:
`BENCHPRO=benchpro just example-minimal`. Without `just`, copy the commands that
`just --show` prints.

| Demo | Recipe | Shows | Time |
|---|---|---|---|
| [01-minimal](01-minimal/) | `just example-minimal` | The smallest meta file (5 columns); taxonomy inferred from the files | < 1 s |
| [02-full](02-full/) | `just example-full` | Every optional column, `--adjusted` and `--allow-alternatives` | < 1 s |
| [03-gtdb-and-ncbi](03-gtdb-and-ncbi/) | `just example-gtdb-and-ncbi` | GTDB and NCBI tools in one run; `Tool` = tool + database | < 1 s |
| [04-mistakes](04-mistakes/) | `just example-mistakes` | Six common mistakes and the message each one produces | < 1 s |
| [05-matrix](05-matrix/) | `just example-matrix` | Abundance matrices (samples as columns) instead of one file per sample | < 1 s |
| [06-other-commands](06-other-commands/) | `just example-other-commands` | `normalize`, `merge` and `strain` | about 3 s |

Output goes to `examples/<demo>/output/` (ignored by git).

If `~/.taxdet/index` exists, Benchpro also reports the release of each profile (see
[Taxonomy detection](../docs/taxonomy-detection.md)); the expected output below was made without
an index.
