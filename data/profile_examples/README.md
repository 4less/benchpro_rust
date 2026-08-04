# Profile format examples

Fixtures for the format-detection and loader tests in `normalize_detect.rs`,
`normalize_loader.rs` and `tests/integration_normalize.rs`. They are referenced by
`include_str!`, so the crate does not compile without them.

Provenance differs, and it changes what the tests using them prove:

| File | Provenance |
|---|---|
| `metaphlan3013_camioutput.profile` | real MetaPhlAn 3.0.13 CAMI output, from `data/test_data/` |
| `motus301_parenthesis.profile` | real mOTUs 3.0.1 output (`-C parenthesis`), from `data/test_data/` |
| `sylph090.profile` | real sylph output, from `data/test_data/` |
| `protal.profile` | real protal output, from `data/test_data/` |
| `bracken.profile` | **synthetic** — no Bracken report exists in this repository |
| `motus404.relab` | **synthetic** — no mOTUs 4.x output exists in this repository |
| `metaphlan402.profile` | **synthetic** — no MetaPhlAn SGB profile with `t__` ranks exists here |

The three synthetic files reproduce the column layout and header lines each detector keys on.
They are honest format exemplars, but a test that reads one is a **regression guard on the
detector**, not evidence that the detector matches what the real tool emits — the file was written
to the format the detector expects, which is circular. Replace them with captured output when any
becomes available, and the tests become validation rather than regression cover.
