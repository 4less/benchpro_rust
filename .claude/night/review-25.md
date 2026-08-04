# Audit: `benchpro align`, `eed0233..HEAD` (checkpoint 25)

Auditor: `night-reviewer` (agent `a89bcef9adab2d0e9`). Verdict: **STOP**.

No drift from TASK.md. `dd6c439`/`9fca287` are pure memory work with no behavioural change.
Everything below is in `fe76565` (the result cache). Both CRITICALs demonstrated on the committed
fixtures with the release binary.

---

## CRITICAL

### C1. The fingerprint stamps the peer's *alignment* but not the peer's *reference*

`cache.rs:179-192` stamps the peer's alignment file only. But `base::head_to_head` chooses
replayed-vs-replayed or reported-vs-reported from `(record.vnm, other.vnm)`, and `other.vnm` exists
only if the *peer's row* carries a `Reference` (`mod.rs:296`). So editing only the peer's row
changes the depending row's h2h while leaving it a cache hit.

```
run 1: no Reference on the sloppy row
run 2: Reference added to the sloppy row ONLY
cached run 2:  good  h2h_better=''   h2h_equal=''     h2h_worse=''   h2h_nm_delta=''
forced run 2:  good  h2h_better=0.0  h2h_equal=100.0  h2h_worse=0.0  h2h_nm_delta=0.0
```

Four summary columns differ between a cached and a forced run — exactly the property the commit
message claims is guaranteed. Also fires when the peer's reference FASTA changes, or when
`Reference` is *removed* from the peer row.

### C2. `--force` neither refreshes nor invalidates the cache

`mod.rs:102` makes every fingerprint `None` under `--force`, so `cache.store` is never called, and
`mod.rs:183-185` skips `cache.save`. The existing cache file is left untouched and the **next**
plain run serves the same stale entries.

```
md5 cache  b45b5b6d...      # run 1
--force                     # correct numbers in the tables
md5 cache  b45b5b6d...      # UNCHANGED
plain run                   # "unchanged, reusing" -- wrong numbers back
```

`cache.rs:93` and `docs/align.md:305` both advertise `--force` as *the* escape hatch for the
size+mtime false hit. It does not stick, which also makes C1 unrecoverable without deleting a file
the docs never mention.

---

## WARNING

- **W1** `FORMAT` is hand-maintained, and the change this module keeps making — a fixed denominator,
  no struct field touched — is invisible to it. A warm run after such a fix serves the old wrong
  number. `Cache`/`CacheEntry`/`SampleResult` also lack `#[serde(deny_unknown_fields)]`: a stored
  result with an obsolete field was reused without complaint. (`cache.rs:86`)
- **W2** mtime is truncated to whole seconds, making the documented false-hit window one second wide
  where nanoseconds are free. (`cache.rs:122`)
- **W3** `Cache::save` is non-atomic (no scratch+rename, unlike `ReadsWriter` and the `.fai`
  builder), and its failure propagates with `?` *after* all four tables were written — an unwritable
  cache turns a successful benchmark into `exit(1)`. (`cache.rs:377-387`, `mod.rs:184`)
- **W4** `vocabulary_mismatch` lives inside the per-row scoring block, which reused rows skip: run 1
  warns that the truth is labelled for the other scoring mode, run 2 on identical inputs is silent.
  (`mod.rs:347-365`)
- **W5** The 250k chunking path is never exercised — the largest fixture yields 20 rows per
  contender, so `chunks()` always returns one. The memory claim rests on untested code.
  (`mod.rs:63,377`)

---

## NOTE

- Peer propagation is one pass in sheet order; safe only because a reused row never recomputes its
  own h2h, but the invariant is unstated. (`mod.rs:235-244`)
- Entries are never pruned; a removed row's entry is written back out forever. (`cache.rs:344-362`)
- `FileStamp.path` is uncanonicalized, so two runs from different working directories with relative
  paths compare equal. (`cache.rs:117`)
- `ReadsWriter`'s docs still say "one group at a time"; it is now one contender, in chunks.
- `the_fixture_set_can_distinguish_every_denominator` guards its own samplesheet, not "the fixture
  set the other tests depend on" as documented; and `num()` maps an empty cell to 0.0, so a
  discriminator can be satisfied vacuously by a missing column.
- `f64_bits` is exact for NaN, ±inf and −0.0. No finding; the ULP rationale checks out.

---

## Verified clean

Cache key uniqueness; `row.scoring`/`row.sep` being the *effective* values, so a changed CLI default
invalidates only rows that did not override it; adding/removing a contender invalidating rows that
name it as a peer; `mappable` genuinely recomputed on every path including the fully-reused early
return; every other field of `SampleResult`/`BaseMetrics`/`ClipGeometry`/`PerGenome` walked to the
input that determines it — the peer's `vnm` (C1) was the only gap; `parsed[index]` alignment;
`truths[&key]` cannot panic; no underflow; cache saved only after all four tables succeed; per-read
row order unchanged by `dd6c439`/`9fca287` and byte-identical across `-t 1`/`-t 8`; `87a6bfa`'s
guard asserts what it claims. 29/29 integration tests pass; clippy clean in `src/align/`.

## Not examined

Parity against the real `SRR8797712` SAMs; `bioreader` internals; RSS for the chunking claim;
concurrent runs sharing an `--outprefix` (the four tables already clobber each other).
