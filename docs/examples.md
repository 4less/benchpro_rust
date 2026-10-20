# Examples

The repository's [`examples/`](https://github.com/4less/benchpro_rust/tree/main/examples) folder
has small, runnable demos built on the test data:

| Demo | Recipe | Shows |
|---|---|---|
| `01-minimal` | `just example-minimal` | The smallest meta file (5 columns); taxonomy inferred from the files |
| `02-full` | `just example-full` | Every optional column, `--adjusted` and `--allow-alternatives` |
| `03-gtdb-and-ncbi` | `just example-gtdb-and-ncbi` | GTDB and NCBI tools in one run; `Tool` = tool + database |
| `04-mistakes` | `just example-mistakes` | Six common mistakes and the message each one produces |
| `05-matrix` | `just example-matrix` | Abundance matrices (samples as columns) |
| `06-other-commands` | `just example-other-commands` | `normalize`, `merge` and `strain` |

From a checkout:

```bash
just example-minimal   # one demo
just examples          # all demos
```

Each demo folder has a README with the expected output.
