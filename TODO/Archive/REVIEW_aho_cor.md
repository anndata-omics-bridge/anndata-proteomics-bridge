# Review — `apb fasta` validation + `prozor` (Aho–Corasick) feature

**Date:** 2026-07-21
**Scope:** uncommitted working-tree changes in `apb/` and `apb_studio/`
**Companion doc:** [TODO/TODO_aho_cor.md](TODO/TODO_aho_cor.md)

> **Resolution:** The vendored algorithm code discussed below was subsequently
> extracted into the standalone Python `prozor` package. APB now depends on
> Prozor and retains only FASTA orchestration and AnnData/MuData storage.

## What was reviewed

Uncommitted changes implementing FASTA peptide-identification validation and a
vendored `prozor` peptide→protein matching package.

**`apb/` — new files**

- `src/anndata_proteomics/annotation/validate_fasta.py` (783 lines) — peptide validation
- `src/anndata_proteomics/fasta/config.py` — typed decoy/contaminant configuration
- `src/anndata_proteomics/fasta/anndata_io.py` — round-trip config through `uns`
- `src/anndata_proteomics/prozor/` — `ahocorasick.py`, `annotate.py`, `greedy.py`, `sparse_matrix.py`
- `tests/test_validate_fasta.py`, `tests/test_fasta_config.py`, `tests/test_prozor.py`

**`apb/` — modified**

- `annotation/var_fasta.py`, `fasta/annotation.py`, `scripts/cli.py` (unified `apb fasta`)
- `pyproject.toml` (+`scipy`, `ahocorapy`, `ahocorasick-rs`), `README.md`, `docs/ARCHITECTURE.md`
- `tests/test_fasta_annotation.py`, `tests/test_annotation_var_fasta.py`, `tests/test_cli_integration.py`

**`apb_studio/` — modified:** doc/TODO relinking only (`AGENTS.md`, `README.md`,
`pipeline.py` comment, `workflow/Snakefile` comment) — fine, no code behaviour change.

## Bottom line

A **well-built feature**. Tests are thorough (126 passing, including h5ad/h5mu
round-trips, backend equivalence, and the tricky MuLink idempotency semantics),
`ruff` is clean, `__pycache__` is correctly gitignored, and the "classify, never
filter" invariant is applied consistently across both FASTA paths.

One consistency issue is worth fixing before commit; the rest are conscious
decisions to confirm. None block the commit.

---

## Findings (ranked)

### 1. Two divergent copies of the same protein-field list — should unify

Directly hits the project rule (*"ensure column naming, parameter naming, and API
surface are CONSISTENT"*).

- [var_fasta.py:38-48](apb/src/anndata_proteomics/annotation/var_fasta.py#L38-L48) — `_PROTEIN_MATCH_FIELDS`
- [validate_fasta.py:61-71](apb/src/anndata_proteomics/annotation/validate_fasta.py#L61-L71) — `_LEADING_PROTEIN_FIELDS`

Same nine elements, **different order**. They diverge at index 1: protein
annotation prefers `PG_ProteinAccessions`, peptide validation prefers
`Leading_Razor_Protein`. Verified — for a var table carrying both columns the two
paths auto-select different columns.

In the normal MuData flow they act on different modalities, so there is no crash;
but the divergence looks unintentional (identical sets) and makes
`peptide_in_leading_protein` reason about a different column than the protein
modality's `match_on`.

**Fix:** define one shared tuple and import it in both (per "reuse before
duplicate"). If the different precedence is *deliberate* (razor-protein-first for
peptides), keep two but add a comment on each explaining why — nothing currently
signals intent.

### 2. `_uniprot_proteinname` can mangle non-decoy IDs — minor robustness

[fasta/annotation.py:143-152](apb/src/anndata_proteomics/fasta/annotation.py#L143-L152).
The prefix-preserving regex `^(.*?)(?:sp|tr)\|([^|]+)\|` correctly keeps decoy
prefixes (`REV_sp|P12345|X → REV_P12345`, `CON__sp|Q9|X → CON__Q9`). But any ID
literally containing `sp|`/`tr|` after some text is rewritten:
`contains_sp|P9|N → contains_P9`. All real decoy/contaminant prefixes in
[config.py:18-32](apb/src/anndata_proteomics/fasta/config.py#L18-L32) are safe
(none end in `sp`/`tr`), so this is edge-case only — note, not a fix.

### 3. `prozor` ships ~440 LOC of unwired inference code — confirm intent

Only `annotate_peptides_streaming` is used by the product (grep-confirmed).
[greedy.py](apb/src/anndata_proteomics/prozor/greedy.py) (241 lines), most of
[sparse_matrix.py](apb/src/anndata_proteomics/prozor/sparse_matrix.py) (201 lines),
and `AnnotationResult.{filter_tryptic,to_dataframe,to_sparse_matrix}` are exercised
only by `test_prozor.py`. [ARCHITECTURE.md](apb/docs/ARCHITECTURE.md) documents this
("currently unwired protein-inference primitives") and `__init__.py` says it is a
deliberate vendor from `diann_runner`, so this is a **decision, not a defect** —
just confirm you want to carry+maintain it now versus adding it when a consumer
exists.

### 4. MuLink ownership machinery persists a hidden second matrix — confirm it's warranted

[validate_fasta.py:678-730](apb/src/anndata_proteomics/annotation/validate_fasta.py#L678-L730)
tracks APB's own edges in a second sparse matrix
`varp['_apb_fasta_feature_mapping_contribution']` so re-runs replace only APB's
edges. It is correct and well-tested (idempotency, cross-producer preservation,
uint64/overflow-safe coords). Two things to weigh:

- That `_apb_*` key is written into every `.h5mu` and will be visible to any
  downstream MuLink reader enumerating `varp`.
- It is substantial complexity for an idempotent-rerun guarantee.

Both may be justified — just make sure the rerun requirement is real enough to pay
for a persisted side-channel matrix.

Minor: `_replace_owned_feature_mapping` relies on canonical (duplicate-free) CSR for
`np.isin(..., assume_unique=True)` at
[validate_fasta.py:713](apb/src/anndata_proteomics/annotation/validate_fasta.py#L713);
safe for scipy-produced / round-tripped matrices, fragile if a non-canonical matrix
is ever injected.

### 5. `prozor` diverges from APB house style

[annotate.py](apb/src/anndata_proteomics/prozor/annotate.py) /
[ahocorasick.py](apb/src/anndata_proteomics/prozor/ahocorasick.py) use
`from typing import Iterable, Iterator` (rest of APB uses `collections.abc`),
Google-style `Args:/Returns:` docstrings, inline `import pandas as pd`, and no
`from __future__ import annotations`. Lint passes, so low priority — but since APB
now *owns* this code, aligning it with the terse modern style in
`validate_fasta.py`/`config.py` would be worthwhile.

### 6. Public `fasta_to_dataframe` behaviour changed (internal-only impact)

[fasta/annotation.py:156](apb/src/anndata_proteomics/fasta/annotation.py#L156):
default `decoy_pattern` went `"^REV_|^rev_"` → `None`, and it now **retains +
classifies** decoys/contaminants instead of dropping them (row count changes, two
new columns `is_decoy`/`is_contaminant`). No external importer found in
`apb_studio`/`ProteoBench`, so risk is contained — noting it for API awareness. Side
effect: the in-silico digest now runs over decoy records too (~2× work on a
target+decoy DB), inherent to the "classify, not filter" design.

---

## What was verified

- `pytest` on all 6 affected test files → **126 passed**; `ruff check` on all changed
  files → **clean**.
- `_iter_sources` handles single inline-FASTA strings correctly (no char iteration),
  so `materialize_sources`/`describe_sources` provenance is sound.
- Decoy/contaminant `is_*` flags, nullable-boolean + categorical dtypes round-trip
  through h5ad/h5mu (covered by tests).
- The CLI reuses the resolved config from protein annotation for peptide validation
  ([cli.py:255-268](apb/src/anndata_proteomics/scripts/cli.py#L255-L268)) so inference
  is consistent across both passes — nicely done.

## Suggested next step

Fix #1 (unify the field list into a single shared constant) — small, safe change.
Everything else is a decision to confirm rather than a defect.
