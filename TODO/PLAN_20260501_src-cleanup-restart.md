# PLAN 2026-05-01 — Clean out `src/` and reset to the RESTART_PLAN skeleton

## Context

The current `src/anndata_proteomics/` is the **pre-restart** code:
- Two parallel systems for vendor handling (`strategies/` *and* `params/`) — the `params/` files are mostly placeholders that say "TODO: implement, see ProteoBench".
- `builder.py` / `core.py` / `proforma.py` were built around the strategy-class abstraction we are walking away from.
- `benchmarks/` and `cli.py` are tied to that old surface.
- Tests, examples, and TOML configs were written for the old shape too.

[docs/RESTART_PLAN.md](../docs/RESTART_PLAN.md) already specifies the **new** target: TOML rules under `parsing_rules/<vendor>/`, validated by pydantic in `rules/`, generic readers in `readers/`, conversion in `converters/`, plus a slim `cli.py`. None of that exists yet.

User intent (this session): **"after looking at the software-specific files I'm inclined to delete all of them"** — i.e. don't try to evolve the old code into the new shape, just delete it and start clean against [RESTART_PLAN.md](../docs/RESTART_PLAN.md).

This is safe because:
1. The repo is a git repo at `main` with a clean working tree (last commit `f6bffda`). Anything we delete is recoverable via `git show <sha>:<path>`.
2. The 52 GB `test_data_download/` cache is gitignored and lives outside `src/` — untouched by this cleanup.
3. The RESTART_PLAN.md and toml_schema.md docs explicitly cover what to rebuild, so we are not throwing away design knowledge by deleting code.

## Recommendation

**Delete all current Python code under `src/anndata_proteomics/` (and the tests/examples that target it). Keep nothing.** Re-create only an empty package skeleton. Then implement against [RESTART_PLAN.md](../docs/RESTART_PLAN.md) in a follow-up session.

Why not "salvage" anything:
- `proforma.py` (262 lines) is the only file that *might* contain non-trivial logic worth keeping. But the RESTART_PLAN explicitly scopes ProForma handling out of the first restart pass (out-of-scope: "second-stage obs annotation, conditions, factors"). When ProForma is needed again, the file is one `git show f6bffda:src/anndata_proteomics/proforma.py` away.
- The 3 TOML configs (`configs/diann.toml`, `maxquant.toml`, `spectronaut.toml`) were written against the *old* schema and are pre-pydantic. The new pydantic schema in [docs/toml_schema.md](../docs/toml_schema.md) deserves fresh TOMLs written against it, not retro-fitted ones.
- The `params/` files are placeholders saying "TODO: implement"; they have no salvage value.
- The strategy classes (`strategies/diann.py` etc.) implement the abstraction the restart is walking away from; copying their logic across would just drag the old shape forward.

In other words: a clean delete + git history is strictly better than a partial rewrite.

## Concrete changes

### 1. Delete

```
src/anndata_proteomics/annotation.py
src/anndata_proteomics/benchmarks/         (whole dir: __init__.py, cli.py, loader.py)
src/anndata_proteomics/builder.py
src/anndata_proteomics/cli.py
src/anndata_proteomics/configs/            (whole dir: diann.toml, maxquant.toml, spectronaut.toml)
src/anndata_proteomics/core.py
src/anndata_proteomics/params/             (whole dir: __init__.py, diann.py, maxquant.py, parameters.py, spectronaut.py)
src/anndata_proteomics/proforma.py
src/anndata_proteomics/strategies/         (whole dir: __init__.py, diann.py, maxquant.py, spectronaut.py)
src/anndata_proteomics/utils.py
tests/conftest.py
tests/test_benchmark_data.py
tests/test_converter.py
tests/test_strategies.py
examples/convert_diann_example.py
examples/convert_spectronaut_example.py
```

Keep:
- `src/anndata_proteomics/__init__.py` — already empty per Coding Rules; leave as the package marker.
- `tests/` directory — keep the empty folder, just remove the `.py` files. Pytest will pass with no tests collected.
- `examples/` directory — same; keep the folder for the new examples to land in.
- `pyproject.toml`, `Snakefile`, `README.md`, `CLAUDE.md`, `AGENTS.md`, `.gitignore`, `docs/`, `TODO/`, `test_data_download/` — all untouched.

### 2. Update [pyproject.toml](../pyproject.toml)

Open it and check:
- Any `[project.scripts]` entries pointing into the deleted modules (`cli.py`, `benchmarks.cli`) — remove them. New entrypoints come back when [RESTART_PLAN.md](../docs/RESTART_PLAN.md) §`cli.py` is implemented.
- Any package-data declarations for `configs/*.toml` — remove. The new home is `parsing_rules/<vendor>/` (RESTART_PLAN §"Proposed Package Layout").

### 3. Update [CLAUDE.md](../CLAUDE.md)

The current "Project Structure" and "Strategy Interface" sections describe the deleted shape. Either:
- (a) replace them with a one-liner pointing to [docs/RESTART_PLAN.md](../docs/RESTART_PLAN.md) as the live source of truth until the new code lands, or
- (b) delete those two sections outright and let RESTART_PLAN.md carry the design.

Recommendation: (a) — leaves a breadcrumb so future Claude sessions don't re-read stale guidance. Concrete wording can be decided when executing.

### 4. Commit

Single commit: `chore: delete pre-restart src; reset to RESTART_PLAN skeleton`.

Body should reference `f6bffda` as the recovery point and link [docs/RESTART_PLAN.md](../docs/RESTART_PLAN.md) as the next-step doc.

## Non-goals (do NOT do in this plan)

- Do **not** start implementing the new `rules/`, `readers/`, `converters/`, or `parsing_rules/` packages yet. That's the next plan, driven by [RESTART_PLAN.md](../docs/RESTART_PLAN.md) §"First Implementation Order".
- Do **not** delete `test_data_download/`. It's the test data tool, gitignored, separately committed.
- Do **not** touch `docs/`. RESTART_PLAN.md and toml_schema.md are the inputs to the next plan.

## Critical files

- [src/anndata_proteomics/](../src/anndata_proteomics/) — most files deleted (see §1)
- [tests/](../tests/) — `.py` files deleted; folder kept
- [examples/](../examples/) — `.py` files deleted; folder kept
- [pyproject.toml](../pyproject.toml) — prune dead entrypoints / package-data refs
- [CLAUDE.md](../CLAUDE.md) — prune stale "Project Structure" + "Strategy Interface" sections

## Verification

```bash
cd /Users/wolski/projects/anndata_bridge/anndata_proteomics_bridge

# 1. Only the empty package marker remains under src/
find src -type f
# expected: just src/anndata_proteomics/__init__.py

# 2. tests/ and examples/ still exist as empty (or near-empty) directories
ls tests examples

# 3. Package still installs (pyproject.toml is consistent)
uv pip install -e .

# 4. Pytest passes with no tests collected (no import errors from stale modules)
pytest tests/ -q

# 5. The package imports cleanly
python -c "import anndata_proteomics; print(anndata_proteomics.__name__)"

# 6. Recovery is one git command away
git show f6bffda:src/anndata_proteomics/proforma.py | head -5
```

## After this plan lands

Open the next plan against [docs/RESTART_PLAN.md](../docs/RESTART_PLAN.md) §"First Implementation Order" — start with `rules/schema.py` (pydantic models) before any vendor-specific work.
