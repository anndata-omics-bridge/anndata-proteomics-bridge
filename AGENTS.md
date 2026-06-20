# anndata_proteomics

Convert proteomics software output to AnnData format.

## Terminology

- **APB** means this project, `anndata_proteomics_bridge`.
- The Python package remains `anndata_proteomics`; use **APB** only as the project-level
  shorthand in plans, architecture notes, and cross-repo migration discussions.

Design lives in the sibling docs repo [anndata_omics_bridge](../anndata_omics_bridge/):
- **[conventions.md](../anndata_omics_bridge/docs/conventions.md)** — column / layer name sanitisation rules (apply on `obs.columns`, `var.columns`, layer names; **not** on `obs_names`/`var_names`/`uns` keys)
- **[adr_tool_specific_views.md](../anndata_omics_bridge/docs/adr_tool_specific_views.md)** — per-tool `uns['<app_name>']['column_roles']` schema (authoritative ADR)
- **[proteomics_rationale.md](../anndata_omics_bridge/docs/proteomics_rationale.md)** — why AnnData for proteomics; ProteoBench / prolfquapp synergies

In-repo docs: [docs/toml_schema.md](docs/toml_schema.md), [docs/RESTART_PLAN.md](docs/RESTART_PLAN.md).

TOML parsing-rule edits must follow [docs/toml_schema.md](docs/toml_schema.md).
In particular, `[columns.*.select]` is for original input-table columns only
(`"<sample>"` is the wide-file exception); APB-derived values such as
`proforma_sequence` and `stripped_sequence` must be declared via
`[[columns.var.compute]]`.

## Current Scope

**Ion/precursor level quantification only:**
- DIA-NN (`report.tsv`)
- MaxQuant (`evidence.txt`)
- Spectronaut (precursor exports)

## Status

The pre-restart `src/` was deleted on 2026-05-01. The package is being rebuilt against [docs/RESTART_PLAN.md](docs/RESTART_PLAN.md) — that doc is the authoritative target architecture (`rules/`, `readers/`, `converters/`, `parsing_rules/<vendor>/`) and implementation order. The TOML rule schema lives in [docs/toml_schema.md](docs/toml_schema.md). Old code is recoverable from git history (last full commit before deletion: `f6bffda`).

## Test Data

ProteoBench test data:
- `/Users/wolski/projects/ProteoBench/test/data/quant/quant_lfq_ion_DIA_AIF/`
- `/Users/wolski/projects/ProteoBench/test/data/quant/quant_lfq_ion_DDA_QExactive/`

## Coding Rules

- **Keep `__init__.py` files empty** (a single module docstring is acceptable). Put classes/functions in separate modules and import them directly from those modules.
- **APB owns reusable proteomics parsing infrastructure.** Modification cleanup/mapping
  rules currently duplicated in ProteoBench per-tool TOMLs should migrate into APB parsing
  TOMLs/schema instead of being reimplemented downstream.
- **Parameter parsing belongs in APB.** ProteoBench parameter parsers should move into APB
  as shared code; ProteoBench should consume APB rather than remain the upstream owner of
  generic vendor parameter parsing.

## Development

```bash
uv venv && source .venv/bin/activate
uv pip install -e .
pytest tests/
```
