# Input Formats

Benchpro uses a meta file to describe which profiles to compare. The meta file can be an `.xlsx`, `.csv`, or `.tsv`.

There are two supported formats:

- **Per-sample profile meta** (one row per sample)
- **Abundance matrix meta** (one row per tool/dataset, expanded into samples via regex)

## Per-sample profile meta

This is the classic format where each row corresponds to one prediction profile and one gold standard profile.

### Required columns

- `ID`: Unique row identifier.
- `Sample`: Sample name.
- `Tool`: Tool name.
- `Profile`: Path to the prediction profile file.
- `GoldStd`: Path to the gold standard profile file.

### Optional columns

Leaving a column out, or leaving a cell empty, never fails silently: Benchpro reports what it
inferred or skipped in `<outprefix>_inputs.tsv` (see [Usage](usage.md#outputs)).

| Column | If missing or empty |
|---|---|
| `Taxonomy` | Inferred per file from its contents; see [Taxonomy detection](taxonomy-detection.md). Values starting with `GTDB` or `NCBI` choose the parsing rules. |
| `GoldStdTree` | `NA`. `--adjusted` skips the row and says so. |
| `AvailableSpecies` | `NA`. Detectability is `Unknown`; `--adjusted` warns that its results are more lenient. |
| `Dataset` | `NA`. |
| `ProfileColumns`, `GoldStdColumns` | Columns are auto-detected. |

### Rules

- Column names are case-sensitive: `GoldStd`, not `goldstd`.
- Required cells must not be empty; empty required cells are rejected with the row number.
- Paths are relative to the folder you run Benchpro from, not to the meta file.

### Minimal example (`meta.tsv`)

```
ID	Sample	Tool	Profile	GoldStd
mpa_s1	s1	MetaPhlAn4	profiles/mpa/s1.tsv	gold/s1.profile
mpa_s2	s2	MetaPhlAn4	profiles/mpa/s2.tsv	gold/s2.profile
```

Run `benchpro profile --meta meta.tsv --validate-meta` to check the file and see the inferred
taxonomy of each row before a full run.

### Column format strings

`ProfileColumns` and `GoldStdColumns` are optional format strings for custom parsing. The format is:

```
<Name>|<Lineage>|<Abundance>|<Rank>
```

Each field is a zero-based column index. `Name` can be `X` if missing. The `Rank` field is optional.

Example from the test data:

```
0|2|4|1
```

## Abundance matrix meta

This format is for a single matrix per tool, where each matrix contains multiple samples as columns. The first column must be a taxon lineage column. Each sample is a column header.

### Columns

- `ID`
- `Dataset`
- `Tool`
- `Taxonomy`
- `Profile`: Path to a prediction abundance matrix.
- `ProfileRegex`: Regex used to match prediction sample column headers.
- `GoldStd`: Path to a gold standard abundance matrix.
- `GoldStdRegex`: Regex used to match gold standard sample column headers.
- `GoldStdTree`
- `AvailableSpecies`

### Regex with capture group

Both `ProfileRegex` and `GoldStdRegex` must contain a capture group. Benchpro uses the capture group value as the sample key and pairs prediction/gold columns that share the same key.

Example:

- Prediction column header: `GAST10`
- Gold standard column header: `gastrointestinal_10`

Regexes:

```
ProfileRegex: ^GAST(\d+)$
GoldStdRegex: ^gastrointestinal_(\d+)$
```

This pairs `GAST10` with `gastrointestinal_10` because both capture `10`.

### How rows are expanded

Benchpro scans the headers of both matrices, matches columns by regex, then expands the single meta row into multiple per-sample rows. Each expanded row receives a unique `ID` of the form:

```
<ID>_<capture>
```

The internal profile references are represented as:

```
matrix:<matrix_path>::<column_name>
```

## Taxonomy parsing

The `Taxonomy` column influences how lineages are parsed:

- Values starting with `GTDB` use GTDB parsing (lineage split on `;` or `|`).
- Values starting with `NCBI` use NCBI parsing (lineage split on `|`).
- Other values are parsed as custom (lineage split on `;`).
