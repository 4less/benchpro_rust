# Audit: `benchpro`, `648761f..HEAD` (checkpoint 30)

Auditor: general-purpose agent (`night-reviewer` no longer available; same adversarial brief).
Five commits, all outside `src/align/`. Build green, fmt clean, 223 unit + 37 integration pass.

## CRITICAL

### C1. `normalize` writes fractions and percentages into the same column, unscaled

`normalize_loader.rs` takes each tool's native value verbatim; `normalize.rs:115` writes it to the
CAMI `PERCENTAGE` column. The tools disagree on units: metaphlan/sylph are percent (0-100), while
mOTUs, bracken and protal are fractions (0-1). A normalized mOTUs profile sums to 1.0 in a column
CAMI defines as 0-100. Scale-invariant metrics hide it; Bray-Curtis and L1/L2 do not.

**Pre-existing** — the loaders predate this range. In scope only because `11d6bc5` commits a fixture
set that encodes both conventions side by side without flagging the conflict.

## WARNING

- **W1** The rewritten merge assertion (`merge.rs:1001,1006`) is a prefix match:
  `contains("\t2.5")` also passes for `2.55`, `2.5000001`. The code is right and the old assertion
  was stale (confirmed), but the replacement is weaker than needed. `\t12.5` does not match, and
  this test has one sample column, so wrong-column/wrong-line are not reachable *here*.
- **W2** `motus404.relab`'s `# tool_version=4.0.4` line is invented — mOTUs writes its version in
  the `# git tag version …` banner. `detect_tool`'s motus branch only greps `tool_version=` and has
  no fallback, so a *genuine* mOTUs relab profile reports version `X`. The synthetic fixture is the
  only reason the test is green: it conceals a real detector gap.
- **W3** A real Bracken report panics `benchpro merge` (`format.rs:594` `.expect()` on a missing
  lineage column) instead of returning `MergeError::ProfileParse`. Pre-existing; newly reachable
  because this range commits the first real-Bracken-shaped file.
- **W4** `merge.rs:320-327` and `:330-339` are both dead (`lineage.is_none()` never holds); the
  good error message written for exactly W3's case never fires.
- **W5** `just check` says "as CI would" but there is no Rust CI (only `docs.yml`), and
  `cargo clippy --all-targets` lacks `-D warnings`, so the lint step cannot fail — ~150 warnings,
  exit 0.
- **W6** `just uninstall` uses the default prefix, so it cannot undo `just install-local`; it exits
  0 having removed nothing.

## NOTE

`metaphlan402.profile` pairs a 2-column GTDB-mode header with pipe-delimited SGB content — a
chimera; the 4-column index path is unexercised. `motus404.relab`'s body could not come from mOTUs
(header shape and GTDB-style lineages). `bracken.profile`'s fractions sum to 0.55, not 1.0.
`keeps_metaphlan_t_rank_when_present`'s second assertion is vacuous (metadata is always empty).
`is_metaphlan_tool`'s second clause is dead. mOTUs' `unassigned` is not canonicalized like
MetaPhlAn's `UNCLASSIFIED`. `prefix` can evaluate to `""` — destructive under sudo.

## Verified clean

**`0af7e30` is provably inert**: rustfmt applied to the parent's `strain.rs` is byte-identical to
HEAD's, and an independent token diff gives exactly four deleted `,` tokens, all trailing commas in
parameter lists. `main.rs` is an alphabetical reorder, both positions after `#[macro_use]`.
Justfile quoting is safe for paths with spaces; `install-local` still targets `~/.local`; `default`
is an improvement (a bare `just` previously ran `install-local`). Every flag in `docs/align.md`
exists in `--help`; every column named there is emitted by `report.rs`; all doc links resolve and
are in the nav. TASK.md §13 divergences 1, 4 and 6 verified against code. The merge *code* was
right and the old test assertion was stale.

## Not examined

The TASK.md parity numbers (needs the Python and SRR8797712). `mkdocs build --strict` (mkdocs not
installed). `src/align/` internals — unchanged in this range, and prior CRITICALs were not
re-checked for regression.
