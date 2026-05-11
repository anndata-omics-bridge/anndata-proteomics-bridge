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

## Workflow Rules

- **Plans live in the project `TODO/` folder, and that file is the source of truth.**
  - **Writing.** When entering plan mode (or when asked to "plan" / "analyze"), write the plan to `TODO/PLAN_YYYYMMDD_<short-slug>.md` at the project root, using today's date and a 2–4 word kebab-case slug (e.g. `TODO/PLAN_20260501_consolidate-agents-claude.md`). Create `TODO/` if it doesn't exist. Claude Code will additionally write its own copy under `~/.claude/plans/`; treat that copy as ephemeral.
  - **Reading / implementing.** When the user asks to implement a plan, read `TODO/PLAN_*.md` from the project root — that is the authoritative file, because the user may have edited it after planning. Do **not** read or rely on `~/.claude/plans/`. If the user references "the plan" without naming a file, use the most recent `TODO/PLAN_*.md` by date in the filename, unless they say otherwise.
  - **If they differ, the project file wins.** Never overwrite the project file from the `~/.claude/plans/` copy — the user's edits to `TODO/PLAN_*.md` are intentional. If you make material plan revisions during a session, update the project file before exiting plan mode.

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
