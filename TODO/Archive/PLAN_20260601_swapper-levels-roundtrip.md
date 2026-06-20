# PLAN 2026-06-01 — protein_swapper ⇄ APB: multi-level capture, MuData, and the vendor round-trip

> **Status: DESIGN / DISCUSSION DRAFT.** This is an analysis document, not an
> approved implementation plan. It deliberately leaves many decisions open (see
> §8). Nothing here should be coded until the open questions are resolved and a
> follow-up implementation `PLAN_*.md` is written per workstream.

## 1. Why this document exists

A new sibling, `protein_swapper` (`../protein_swapper/`), was added to the
workspace. It generates an in-silico differential-abundance benchmark by
*swapping* the measured values of rank-aligned high/low protein pairs in one run
group (G2), leaving identities fixed. To do that it has independently grown two
capabilities that APB does **not** have today:

1. **Multi-level capture in one pass.** It decomposes a single Spectronaut
   report into *all five* hierarchy levels at once —
   `protein → peptide → peptidoform → precursor → fragment` — and holds them
   together (`protein_swapper/src/protein_swapper/datamodel.py`,
   `LEVELS`).
2. **A lossless vendor round-trip.** vendor TSV → internal model → vendor TSV,
   byte-faithful enough that `protein-swapper check` asserts equality
   (`protein_swapper/src/protein_swapper/io_spectronaut.py` + `anndata_io.py`).

protein_swapper *also already serialises each level to its own AnnData*
(`anndata_io.level_to_anndata`: `X` = ranking matrix, `layers` = other
measurement matrices, `obs` = identity, `var` = samples). In other words it has
re-implemented, for Spectronaut only, a chunk of exactly what APB is supposed to
own. Per the APB Coding Rules ("APB owns reusable proteomics parsing
infrastructure"; "Reuse before duplicate"), the right end-state is:

> **protein_swapper consumes APB for vendor ⇄ AnnData I/O, and APB grows the
> multi-level + round-trip capabilities to support it.**

This document maps the gap and sketches the workstreams to close it.

## 2. Where the two projects stand today

| Capability | APB (`anndata_proteomics`) | protein_swapper |
|---|---|---|
| Direction | vendor → AnnData **only** (RESTART_PLAN explicitly stops here) | vendor ⇄ AnnData ⇄ vendor (lossless round-trip) |
| Levels per run | **one** level per file, selected by the rule's `quantification_level` | **all five** levels simultaneously, one AnnData each |
| Container | single `AnnData` | a *bundle*: `level_<name>.h5ad` ×N + sidecar parquets |
| Vendors | DIA-NN, Spectronaut, MaxQuant, FragPipe, PEAKS, WOMBAT (TOML-driven) | Spectronaut only (hand-written `io_spectronaut.py`) |
| Rule system | pydantic `ParseRule` + packaged TOMLs + header recognition | none — vendor knowledge hard-coded |
| Sparsity/round-trip metadata | none | `observed.parquet` (row spine) + `report_columns.parquet` (column order) |
| Axis orientation | `obs` = feature identity, `var` = samples (pivot index/cols) | **same**: `obs` = ids, `var` = samples |

The orientation match is the key enabler: both put feature identity on `obs` and
samples on `var`, so an APB-produced `AnnData` is shape-compatible with what
protein_swapper's swap algorithm expects to operate on.

Note APB has already grown past the original RESTART_PLAN scope (`modifications/`,
`params/` subsystems, `quantification_level` is already a `ParseRule` field). So
expanding scope again is consistent with how the package has actually evolved —
but it **does** contradict RESTART_PLAN §"Non-Goals" ("vendor file + parsing TOML
→ AnnData" was meant to be the stopping point). That tension is a decision for
§8, not something to silently override.

## 3. Target picture

```
                 ┌─────────────────────── APB ───────────────────────┐
 vendor file ──► readers ──► converters ──► multi-level container ──► writers ──► vendor file
 (Spectronaut,                (TOML rules)     (MuData? bundle?)        (NEW)
  DIA-NN, …)                                        │
                                                     ▼
                                          protein_swapper.swap()
                                          (operates on the container,
                                           emits ground-truth tables)
```

protein_swapper keeps *only* the swap algorithm, pairing logic, ground-truth
generation, and CLI. Everything format-facing (read, decompose into levels,
hold the levels, write back) becomes APB.

## 4. The four capability gaps to close in APB

### 4A. Multi-level capture ("one file → many levels")

Today one rule produces one level. The hierarchy levels are *not independent
files* — they are all derivable from the same long vendor table by choosing
different `obs` identity keys (a fragment row already carries its precursor /
peptidoform / peptide / protein ids; cf. protein_swapper's `row_meta`).

Design options:

- **(a) Multiple rules, one per level**, sharing a vendor — reuses the existing
  `parse_<software>_<level>_<version>.toml` convention 1:1, and the converter
  runs each rule against the same loaded DataFrame. Simple, no schema change;
  cost is N re-pivots of the same table and no single object that knows the
  levels belong together.
- **(b) One rule that declares a level hierarchy** — a new `[[levels]]` array in
  the TOML, each entry naming its `obs` identity keys and ranking layer. One
  parse, one assemble pass, levels emitted together. Bigger schema change, but
  matches the "they belong together" reality and is the natural input to a
  MuData container.

Recommendation to debate: **(b)**, because the round-trip and the swap both need
the levels as a coherent set, and because parent-id columns must be consistent
across levels (a fragment's `protein` id must equal that protein level's id).

### 4B. A multi-level container — MuData, or APB's own bundle?

protein_swapper today uses a **directory bundle** (`anndata_io.write_bundle`):
one `.h5ad` per level plus parquet sidecars. The user raised **MuData** as the
candidate first-class container. The orientation question is decisive here:

- MuData (muon/scverse) conventionally holds modalities that **share `obs`**
  (the same cells across RNA/ATAC/protein). Here the levels share **samples**,
  which both projects currently place on `var`. So levels-as-modalities would
  share `var`, not `obs` — the *inverted* MuData layout.
- Two ways out: **(i)** accept a shared-`var` MuData (works, less idiomatic), or
  **(ii)** transpose the whole convention so samples live on `obs` and features
  on `var` (scverse-idiomatic, MuData-native — but flips APB's current
  orientation and every existing test/rule, and contradicts how proteomics
  people read these matrices).

This is the single biggest open question (see §8 Q1/Q2). It affects the rule
schema, every converter, protein_swapper's swap code, and downstream readers
(prolfqua/MSstats etc. expect features × samples).

A pragmatic middle path: keep per-level `AnnData` with the current orientation
**inside** a MuData, with samples on a shared dimension, and document the
non-idiomatic axis explicitly in `uns`. Defer the full transpose unless a
concrete downstream consumer demands it.

### 4C. The reverse direction: AnnData/MuData → vendor file

This is brand-new for APB and the hardest part to get *lossless*. A TOML rule
today is a one-way mapping (vendor columns → obs/var/layers). To write back you
need the inverse plus everything the forward direction discards:

- **Original column order and full column set.** protein_swapper stores
  `report_columns` so it can re-emit columns in the exact source order, including
  vendor columns that were never mapped into a layer. APB currently keeps the
  rule JSON in `uns` but not the untouched passthrough columns.
- **The sparsity spine.** A wide `obs × var` matrix is dense; the source long
  table had only the observed `(identity, sample)` rows. Writing back a dense
  matrix would invent rows that never existed. protein_swapper solves this with
  `observed.parquet`. APB has no equivalent.
- **Exact dtypes / null-vs-NaN.** protein_swapper stores a per-matrix
  `matrix_schema` in `uns` and maps float `NaN` back to `null` on the way out
  (`anndata_io._array_to_matrix`). APB would need the same to round-trip
  faithfully.

Design sketch: a `writers/` subsystem mirroring `readers/`, driven by the *same*
`ParseRule` read in reverse, plus a small amount of round-trip metadata persisted
in `uns` (column order, passthrough columns, dtype map) and an observed-row spine
carried alongside (its own `uns` table or a sidecar). A `check`-style
round-trip test per vendor, exactly like `protein-swapper check`, gates it.

### 4D. Round-trip provenance metadata

Independently of the container choice, lossless round-trip needs APB to persist,
per converted object: original column order, unmapped passthrough columns, the
observed `(identity, sample)` spine, and the exact per-layer dtype map. Decide
*where* (single `uns` blob vs. sidecar parquets vs. MuData `.uns`) once the
container in 4B is settled.

## 5. protein_swapper migration (the consuming side)

Once 4A–4D exist:

1. Replace `io_spectronaut.to_internal` / `from_internal` with calls into APB's
   multi-level converter + writer.
2. Replace `anndata_io.write_bundle` / `read_bundle` with APB's container I/O.
3. Keep `datamodel.SwapData` only if its shape still buys the swap algorithm
   something APB's container doesn't; otherwise have `swap()` operate directly on
   the APB container. (Prefer deleting `datamodel.py` + `io_spectronaut.py` +
   `anndata_io.py` if APB fully subsumes them — that is the whole point.)
4. protein_swapper's CLI `import` / `export` / `check` become thin wrappers over
   APB; `swap`, `subset`, ground-truth, and pairing stay protein_swapper-only.

A migration is only worthwhile if APB ends up **lossless for Spectronaut** —
protein_swapper's existing `check` test is the acceptance bar. If APB can't meet
it, protein_swapper keeps its own I/O and we reconsider.

## 6. Suggested phasing (each its own future PLAN)

1. **Spike: APB Spectronaut round-trip at the ion level.** Add a `writers/`
   prototype + round-trip metadata for the *one* vendor/level protein_swapper
   needs most, and port protein_swapper's `check` test. Proves 4C/4D before any
   schema churn. Lowest-risk first step.
2. **Multi-level rule schema** (4A) — pick option (a) or (b); implement and test
   against Spectronaut all five levels.
3. **Container decision + implementation** (4B) — MuData vs. bundle, orientation.
4. **Generalise the round-trip** to the other vendors APB already reads.
5. **protein_swapper migration** (§5) — flip it to consume APB, delete the
   duplicated I/O, keep `check` green.

Phase 1 is deliberately a throwaway-friendly spike: it answers "can APB be
lossless at all?" before we pay for schema and container changes.

## 7. Risks & non-goals

- **Scope creep vs. RESTART_PLAN.** This expands APB past its stated stopping
  point. If the team would rather keep APB strictly one-directional, the
  alternative is a *new* shared package (e.g. `anndata_proteomics_roundtrip`)
  that depends on APB for the forward direction — keeping APB's core minimal.
  (See §8 Q4.)
- **Orientation churn.** A transpose to scverse-idiomatic axes (4B-ii) would
  touch every rule, converter, and test. Treat as a separate, explicit decision.
- **Lossless is expensive.** Byte-exact round-trip across six vendors is a large
  surface. Scoping the guarantee to "round-trips the columns the rule maps +
  faithfully passes through the rest" may be the pragmatic contract.
- **Not in scope here:** the swap algorithm itself, ground-truth semantics,
  pairing/stratification — those stay in protein_swapper and are out of scope for
  APB regardless.

## 8. Open questions (must resolve before any implementation plan)

1. **Container:** MuData, or keep protein_swapper's directory-bundle pattern in
   APB? If MuData, do modalities share `var` (samples) under the current
   orientation, or do we transpose?
2. **Orientation:** keep `obs` = features / `var` = samples (proteomics-readable,
   matches both projects today, non-idiomatic for MuData), or transpose to
   scverse-native (`obs` = samples)? This blocks 4B.
3. **Rule shape for levels:** one rule per level (4A-a) or one rule with a
   `[[levels]]` hierarchy (4A-b)?
4. **Home for the round-trip:** grow APB to be bidirectional, or add a separate
   round-trip package that builds on APB's forward direction? (Honours
   RESTART_PLAN's non-goals either way.)
5. **Losslessness contract:** byte-exact (protein_swapper's current `check`
   bar), or "mapped columns exact + passthrough preserved"? Defines the test.
6. **Levels scope:** must APB capture all five levels, or is
   ion + peptidoform + protein enough for the benchmark protein_swapper builds?
7. **Migration appetite:** is protein_swapper willing to depend on APB now
   (tighter coupling, shared release cadence), or should APB stabilise the
   round-trip first behind its own tests and protein_swapper migrate later?

## 9. Pointers

- protein_swapper internals: `../protein_swapper/AGENTS.md`,
  `../protein_swapper/docs/algorithm.md`,
  `../protein_swapper/src/protein_swapper/{datamodel,anndata_io,io_spectronaut}.py`.
- APB architecture: `docs/RESTART_PLAN.md`, `docs/toml_schema.md`,
  `src/anndata_proteomics/{rules,readers,converters}/`.
- Why AnnData for proteomics:
  `../anndata_omics_bridge/docs/proteomics_rationale.md`.
