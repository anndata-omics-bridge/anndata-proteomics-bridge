# DONE: enable strict Pyright in APB

> Completed 2026-07-23 with Pyright 1.1.411. The configured `src/` and
> `tests/` scope is green and enforced by pre-commit.

## Enforced target

APB uses strict first-party checking while suppressing diagnostics whose only
signal is a propagated `Unknown` from incompletely typed scientific libraries:

```toml
[tool.pyright]
include = ["src", "tests"]
executionEnvironments = [
    { root = "tests", reportPrivateUsage = "none" },
    { root = "src" },
]
venvPath = "."
venv = ".venv"
pythonVersion = "3.13"
typeCheckingMode = "strict"
reportImportCycles = "error"
reportMissingTypeStubs = "none"
reportUnnecessaryTypeIgnoreComment = "error"
reportImplicitOverride = "error"
enableTypeIgnoreComments = false

reportUnknownMemberType = "none"
reportUnknownVariableType = "none"
reportUnknownArgumentType = "none"
reportUnknownParameterType = "none"
reportUnknownLambdaType = "none"
```

All other strict diagnostics remain enabled, including
`reportMissingTypeArgument`, missing parameter annotations, argument/return
compatibility, optional access, index/call issues, production private
cross-module usage, and import cycles. Tests may directly exercise private
helpers without weakening the production check.

The development environment includes the canonical `pandas-stubs` and
`scipy-stubs` packages. APB-owned library boundaries use precise annotations,
runtime narrowing, Pydantic validation, or small structural protocols rather
than blanket `Any`, casts, or checker ignores.

## Scope decisions

- `src/` and `tests/` are both checked.
- `src/anndata_proteomics/scripts/extract_raw_file_db.py` is included.
- Tests use typed fixtures and real guards where pandas/AnnData APIs return
  unions.
- Cross-module helpers that are intentionally consumed elsewhere have public
  names; genuinely internal helpers remain private.

## Enforcement

`.pre-commit-config.yaml` runs:

```text
uv run --extra dev pyright
```

with `pass_filenames: false`, so every commit checks the full configured scope.
A separate `pygrep` hook rejects both `# type: ignore` and
`# pyright: ignore` comments.

## Completion checks

```text
uv run --extra dev pyright
0 errors, 0 warnings, 0 informations

pytest -q
471 passed, 4 skipped

ruff check .
All checks passed

ruff format --check .
119 files already formatted
```

## Optional future work

Absolute strictness for every third-party `Unknown` remains an optional stretch
goal. It should start by improving or sourcing accurate AnnData/MuData stubs,
then re-enabling individual `reportUnknown*` rules only when they produce
maintainable signal.

Ruff's broader `ANN, B, C4, C90, I, PGH, PIE, RUF, SIM, UP` policy remains a
separate change because it is independent of the completed Pyright gate.
