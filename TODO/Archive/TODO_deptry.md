# DONE: enable Deptry in APB pre-commit

> Completed 2026-07-23 with Deptry 0.25.1. Dependency declarations are clean
> and enforced for the complete APB package by pre-commit.
>
> Archived 2026-07-24 after the complete quality gate passed.

## Changes

- Added `deptry>=0.24,<1` to the locked development environment.
- Declared the directly imported `pyyaml` and `packaging` runtime dependencies.
- Configured `anndata_proteomics` as first-party and the `mkdocs-material`
  distribution-to-module mapping.
- Retained `pyarrow` as a runtime dependency for pandas Parquet I/O, with a
  package-specific `DEP002` exception because pandas loads the engine
  dynamically.
- Added `uv run --extra dev deptry .` as a full-scope local pre-commit hook with
  `pass_filenames: false`.

## Completion checks

```text
uv run --extra dev deptry .
Success! No dependency issues found.

uv run --extra dev pyright
0 errors, 0 warnings, 0 informations

ruff check .
All checks passed

ruff format --check .
119 files already formatted

pytest -q
471 passed, 4 skipped
```
