# PLAN 2026-05-04 — Scaffold the `annProtSum` R package

## Context

The conversion pipeline produces `.h5ad` files; we want a polished HTML report with per-conversion summaries (skimr-style, DataExplorer-style auto-EDA, GGally `ggpairs` for `obs`). The user prefers R for this layer because of the tidyverse ecosystem (gt, skimr, GGally, DataExplorer).

This plan sets up the **scaffolding** for an R package living inside the Python repo at `tools/annProtSum/`. The actual `annprot_report()` function and the Quarto template land in a follow-up. Goal of this plan: a buildable, lintable, testable R package skeleton — so subsequent work has a clean home.

## Files / structure

```
tools/
└── annProtSum/                     R package root
    ├── DESCRIPTION                 (usethis-generated; deps populated)
    ├── NAMESPACE                   (roxygen-generated)
    ├── LICENSE / LICENSE.md        (MIT, via use_mit_license)
    ├── README.md                   (short — points back to project root)
    ├── R/
    │   └── annprot.R               package-level docstring + one stub function
    ├── tests/testthat/
    │   ├── helper-fixtures.R       (placeholder)
    │   └── test-annprot.R          smoke test
    ├── tests/testthat.R            (testthat boilerplate)
    └── inst/
        └── quarto/
            └── report.qmd          empty placeholder for the upcoming template

tools/install_r_deps.R              installs the not-yet-present CRAN deps
```

## R-side dependencies

Imports (DESCRIPTION):
- `anndataR` — native h5ad reader
- `skimr` — per-table descriptive stats
- `GGally` — `ggpairs()` for the obs panel
- `DataExplorer` — auto-EDA HTML report helpers
- `gt` — presentation tables for the manifest
- `here` — project-relative paths
- `readr` — manifest.csv ingest

Suggests:
- `testthat` (≥ 3.0.0)
- `quarto` — for rendering the .qmd template

Already installed (verified): `anndataR 1.1.0`, `GGally 2.4.0`, `here 1.0.2`, `readr 2.2.0`, `rmarkdown 2.31`, `knitr 1.51`, `devtools 2.4.6`, `testthat 3.3.2`, `usethis 3.2.1`.

Need install: `skimr`, `DataExplorer`, `gt`. `tools/install_r_deps.R` handles that idempotently.

## Build steps

1. `R -q -e 'usethis::create_package("tools/annProtSum", open = FALSE, rstudio = FALSE)'` — scaffold the package. Sets DESCRIPTION (placeholder), NAMESPACE, R/, .Rbuildignore.
2. `usethis::use_mit_license("Witold Wolski")` — license + LICENSE.md.
3. `usethis::use_package(pkg, "Imports")` for each import above.
4. `usethis::use_testthat(3)` — set up testthat.
5. Edit DESCRIPTION manually for Title, Description, Author/email (usethis won't infer those).
6. Write `R/annprot.R` with a package-level roxygen doc and one stub function (`annprot_summarize_skeleton(h5ad_path)` that returns a placeholder list — proves anndataR loads).
7. Write `tests/testthat/test-annprot.R` calling the stub against a known-good `.h5ad` in `examples/results/` (skip if absent, mirroring the Python skip-when-cache-missing pattern).
8. `R -q -e 'devtools::document(pkg = "tools/annProtSum")'` to generate NAMESPACE.
9. `R -q -e 'devtools::check(pkg = "tools/annProtSum", error_on = "error")'` to verify no errors before commit. Warnings about missing deps (skimr/DataExplorer/gt not yet installed) are expected and acceptable for the scaffold.
10. Commit.

## Out of scope (follow-up)

- The actual `annprot_report()` rendering function (uses skimr + GGally + DataExplorer + gt to assemble per-h5ad sections).
- The Quarto template `inst/quarto/report.qmd` — created empty here, populated next.
- The Python-side `tools/generate_h5ads.py` that writes `examples/results/*.h5ad` + `manifest.csv` (the input to `annprot_report`). Will be a sibling commit.
- pkgdown site / CI integration.

## Verification

```bash
# 1. Package scaffolds
ls tools/annProtSum/{DESCRIPTION,NAMESPACE,R}

# 2. R CMD check passes (or warnings only — missing-deps warnings are expected
#    until install_r_deps.R is run)
R -q -e 'devtools::check("tools/annProtSum", error_on = "error")'

# 3. After installing CRAN deps, R CMD check is clean
Rscript tools/install_r_deps.R
R -q -e 'devtools::check("tools/annProtSum", error_on = "warning")'

# 4. testthat smoke (skips if no h5ad fixtures present)
R -q -e 'devtools::test("tools/annProtSum")'
```

## Notes

- **Why inside the Python repo**: the R package's job is *visualizing the Python package's output*. Co-locating keeps the report regen step (`tools/generate_h5ads.py` → manifest → R report) close. The R package is small enough that publishing separately to CRAN/Bioconductor isn't on the cards; if it ever is, splitting to its own repo is mechanical.
- **`tools/` is gitignore-clean for outputs only**: the source under `tools/annProtSum/` IS tracked. `examples/results/` (output) IS gitignored.
