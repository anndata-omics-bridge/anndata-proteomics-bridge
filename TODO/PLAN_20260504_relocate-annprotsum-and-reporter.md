# PLAN 2026-05-04 — Relocate `annProtSum`; build a single-file convert+report tool

## Context

Two coupled changes:

1. **Relocate the R package**: `anndata_proteomics_bridge/tools/annProtSum/` → `~/projects/anndata_bridge/annProtSum/` (sibling of `anndata_proteomics_bridge/`, alongside `ProteoBench/` and `anndata_omics_bridge/`). The R package is its own concern and shouldn't live inside the Python repo.

2. **Wire a single-file conversion + report pipeline**:
   - `python tools/generate_report.py <vendor file>` from inside `anndata_proteomics_bridge/`.
   - Python: read → recognize/load → convert → write `<software>_<sha8(input_path)>.h5ad` + sidecar `.meta.json` under `examples/results/`.
   - Python: shell out to a script in `annProtSum` (`inst/bin/render_report.R`) which renders an HTML report from the .h5ad via a parametrized Quarto vignette.
   - Python: rebuild `examples/results/index.html` from all sidecar `.meta.json` files: one row per conversion, columns **software | input | output (.h5ad) | layers (shape) | report**.

The "rich" report content (skim/ggpairs/DataExplorer panels) is **out of scope here**. The vignette ships minimal (shape + per-layer table + obs/var head) so the plumbing is verifiable end-to-end. Filling the vignette is a follow-up plan.

## Layout after the change

```
~/projects/anndata_bridge/
├── anndata_proteomics_bridge/
│   └── tools/
│       └── generate_report.py         orchestrator
└── annProtSum/                        sibling R package
    ├── DESCRIPTION
    ├── inst/
    │   ├── quarto/report.qmd          parametrized vignette (params.h5ad)
    │   └── bin/render_report.R        CLI wrapper
    └── tests/testthat/
```

`tools/install_r_deps.R` is **dropped** — `devtools::install(path, dependencies = TRUE)` already resolves Imports + Suggests. Don't reimplement what the language already does (this principle is now codified in the workspace-level CLAUDE.md).

## Decisions resolved

1. annProtSum target: `/Users/wolski/projects/anndata_bridge/annProtSum/`.
2. Filename stem: `<software>_<sha8(input_path)>` — deterministic, collision-safe.
3. `tools/install_r_deps.R`: removed. Use `devtools::install(..., dependencies = TRUE)`.

## Mechanics

### Move

```
mv anndata_proteomics_bridge/tools/annProtSum  ../annProtSum
rm  anndata_proteomics_bridge/tools/install_r_deps.R
rmdir anndata_proteomics_bridge/tools     # empty after move; tools/generate_report.py recreates it next
```

### Parametrized Quarto vignette `annProtSum/inst/quarto/report.qmd`

- YAML `params: { h5ad: "" }`.
- Setup chunk reads the file via `anndataR::read_h5ad(params$h5ad)`.
- Body: shape line, per-layer shape table (`gt`), `head(obs, 6)`, `head(var, 6)`, `uns$anndata_proteomics` summary.

### R-side wrapper `annProtSum/inst/bin/render_report.R`

```r
#!/usr/bin/env Rscript
# Usage: render_report.R <h5ad-path> <output-html>
```

- Locates the .qmd via `system.file("quarto/report.qmd", package = "annProtSum")`.
- Copies to a tempdir, calls `quarto::quarto_render(execute_params = list(h5ad = ...))`, copies the rendered `report.html` to the requested output path.

Adds `quarto` to `Suggests` in DESCRIPTION (only used by the inst/bin/ script, not in package code).

### Python orchestrator `anndata_proteomics_bridge/tools/generate_report.py`

- argparse CLI: `<data> [--rule-toml] [--output-dir]`.
- Pipeline: `read_table → recognize/load_rule → convert(df, rule) → write h5ad + meta.json → Rscript render_report.R → rebuild_index()`.
- Resolves `render_report.R` via `R -q -s -e 'cat(system.file(...))'`; falls back to dev path `<workspace>/annProtSum/inst/bin/render_report.R`.

`rebuild_index(output_dir)` globs `*.meta.json`, sorts, writes `index.html`. Stateless — every call regenerates from scratch, no append issues.

## Tests

- **R**: `annProtSum/tests/testthat/test-render.R` — happy path against real h5ad fixture if present (skip otherwise).
- **Python**: `anndata_proteomics_bridge/tests/test_generate_report.py` — `rebuild_index` with no metas; full pipeline against a synthetic DIA-NN-shaped TSV. Skips if Rscript / annProtSum unavailable.

## Verification

```bash
# Smoke test against real WOMBAT data
.venv/bin/python tools/generate_report.py \
    test_data_download/json_dir/Results_quant_peptidoform_DDA/.../input_file.csv
ls examples/results/    # wombat_<hash>.{h5ad, html, meta.json} + index.html

# Full suite
.venv/bin/python -m pytest tests/ -q     # 103 tests pass
R -e 'devtools::check("/Users/wolski/projects/anndata_bridge/annProtSum")'  # 0 errors, 0 warnings, 2 (expected) notes

# View
open examples/results/index.html
```

## Files touched

**Created (in annProtSum/)**
- `inst/quarto/report.qmd` (overwritten — placeholder → real parametrized vignette)
- `inst/bin/render_report.R`
- `tests/testthat/test-render.R`

**Modified (in annProtSum/)**
- `DESCRIPTION` — `quarto` in Suggests.
- `README.md` — points at the now-sibling layout.

**Created (in anndata_proteomics_bridge/)**
- `tools/generate_report.py`
- `tests/test_generate_report.py`

**Modified (in anndata_proteomics_bridge/)**
- `docs/ARCHITECTURE.md` — annProtSum is now a sibling, document orchestrator + index.html convention.
- `.gitignore` — drop now-irrelevant `tools/annProtSum/` entries.

**Deleted**
- `anndata_proteomics_bridge/tools/install_r_deps.R`.
- `anndata_proteomics_bridge/tools/annProtSum/` (moved out).

## Out of scope (follow-ups)

- Rich vignette content (skim, ggpairs, DataExplorer, gt panels).
- Multi-file batch mode.
- Initialising annProtSum as its own git repo.
- pkgdown site for annProtSum.
