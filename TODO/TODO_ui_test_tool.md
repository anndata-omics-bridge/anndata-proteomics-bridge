# TODO: test-data browser → convert → inspect GUI

**Date:** 2026-06-22
**Status:** Design / spec — no code yet.
**Related:** [TODO_viewer.md](TODO_viewer.md) (the larger `APB_viewer` product + viz-scaling
research), [TODO_to_mu_data.md](TODO_to_mu_data.md) (the multi-level/MuData conversion this drives).

> **Note:** the existing single-file viewer was moved
> `anndataview.py` → `src/anndata_proteomics/scripts/anndataview.py` (per the same request).
> Its overview/stats/heatmap cells are the basis for the **summary panel** below.

## Does this make sense? — yes

It's a clean, useful internal tool: *browse the ProteoBench test corpus → pick a dataset and a
target → convert → inspect the result*. It also doubles as a live exercise of the registry +
converters on real, varied vendor files. Two refinements worth baking in from the start
(both grounded in what we just learned):

1. **Convertibility must be registry-driven, not a hardcoded vendor list** — so the tool stays
   correct as rules are added. The filter logic you described maps exactly onto
   `rules.registry.find_rule(software, level)` + `converters.recognize.matches(...)`.
2. **Convert is sometimes expensive** — a 740 MB / 175k-precursor DIA-NN file converted to
   **fragment** or **MuData** is the multi-GB / many-second case from the perf work. The tool
   needs a size guardrail + background conversion, not a naive synchronous call.

## Concept

Three-pane single-page app:

```
┌─────────────────────────────────────────────┬──────────────────────────────────┐
│ TARGET:  ( ) MuData (all levels)             │  CONVERTED (this session)        │
│          (•) AnnData level: [ion ▾]          │  ┌────────────────────────────┐  │
│ FILTER:  software [All ▾]  size ≤ [====o] MB │  │ dataset │ target │ shape │…│  │
│                                              │  │ DIANN…  │ ion    │ 6×73k │  │  │
│  LEFT — dataset catalog (filtered)           │  │ MaxQ…   │ ion    │ 6×41k │  │  │
│  ┌────────────────────────────────────────┐ │  └────────────────────────────┘  │
│  │ software │ version │ nr_prec │ size │ ✓ │ │             ▲ select a row        │
│  │ DIA-NN   │ 1.8.1   │ 175 476 │ 740M │ ✓ │ │   ┌──────────────────────────┐   │
│  │ DIA-NN   │ 2.0     │  69 629 │  95M │ ✓ │ │   │ SUMMARY PANEL            │   │
│  │ …        │         │         │      │   │ │   │ (anndataview.py views:   │   │
│  └────────────────────────────────────────┘ │   │  overview, obs/var/X,    │   │
│        ▲ select a row → [ Convert ]          │   │  layers, uns, heatmap;   │   │
│                                              │   │  per-modality for MuData)│   │
│                                              │   └──────────────────────────┘   │
└─────────────────────────────────────────────┴──────────────────────────────────┘
```

## Data source — the catalog (left table)

`test_data_download/raw_file_db_downloaded.csv` (83 rows today, all `status="ok"`). Useful columns:

| Column | Use |
|---|---|
| `software_name` | display + filter + registry slug lookup |
| `software_version` | display |
| `nr_prec` | "size" proxy (precursor count, 2 844 – 175 476) |
| `input_file_size_bytes` | "size" filter (2.3 MB – 740 MB) |
| `input_file_path` | path under `test_data_download/json_dir/…` to feed the converter |
| `status` | only show `ok` |

Software families present (count): DIA-NN 20, MaxQuant 11, i2MassChroQ 8, AlphaDIA 8,
Spectronaut 8, FragPipe (DIA-NN quant) 7, PEAKS 7, FragPipe 5, Sage 2, WOMBAT 2, AlphaPept 2,
quantms 1, ProlineStudio 1, MSAngel 1.

## Target selector + filter logic (the core interaction)

Choosing the **target** filters the left table to datasets whose software has the required
rule(s). Driven by the registry, not hardcoded:

- **MuData (all levels)** → software must have **all five** level rules → **DIA-NN only** (20).
- **AnnData level = ion** → DIA-NN, MaxQuant, Spectronaut, PEAKS, FragPipe (~51 rows).
- **= peptidoform** → DIA-NN, WOMBAT.
- **= peptide / protein / fragment** → DIA-NN only (today).

Current packaged coverage (the `✓` column / what's convertible):

| Software (catalog name) | slug | convertible targets |
|---|---|---|
| DIA-NN | `diann` | ion, peptidoform, peptide, protein, fragment, **MuData** |
| MaxQuant | `maxquant` | ion |
| Spectronaut | `spectronaut` | ion |
| PEAKS | `peaks` | ion |
| FragPipe | `fragpipe` | ion |
| WOMBAT | `wombat` | peptidoform |
| FragPipe (DIA-NN quant) | ? | **resolve via `recognize`** (column variant — TBD) |
| i2MassChroQ, AlphaDIA, AlphaPept, Sage, quantms, ProlineStudio, MSAngel | — | none yet |

**Implementation:** for each catalog row, map `software_name` → registry slug (reuse the
existing alias resolution, e.g. `params.registry` / `get_parser`; do **not** write a new map),
then a target is convertible iff `find_rule(slug, level)` resolves for the level(s). Optionally
confirm with `recognize.matches(headers, rule)` (a cheap header read) to catch DIA-NN column
variants — recall some cached DIA-NN files lack `Fragment.*` or `PG.Normalised`, so "DIA-NN"
does **not** guarantee fragment/MuData for *that* file.

## Convert action

Selected left row + target → **Convert** →
- read the file (`readers.dispatch.read_table`), `load_packaged_rule(slug, level)`,
  `converters.assemble.convert(...)`; for MuData, convert each level and assemble (the
  `tests/test_mudata_levels.py` assembly is the reference pattern to factor into a helper).
- append a row to the **right table**.

**Guardrails (required, not optional):**
- **Size/memory warning** before converting `fragment` or `MuData` above a threshold (the perf
  note: full 6-run fragment ≈ 6.5 GB). Show the estimate; offer **per-run** or **row-capped**
  conversion for big files.
- **Run conversion in the background** with a spinner/progress; never block the UI on a
  multi-second/GB convert.
- Decide **artifact storage**: in-memory dict keyed by `(dataset, target)` (simplest, RAM-bound)
  vs. write `.h5ad`/`.h5mu` to a gitignored cache (`test_data_download/_converted/`) so the
  summary panel reloads lazily. Recommend **disk cache** given the sizes.

## Right table — converted artifacts (this session)

Columns: dataset (software+version), target, shape (`n_obs × n_var`; per-modality for MuData),
layers, build time, peak memory (nice-to-have), artifact path. Select a row → drives the panel.

## Summary panel — reuse `anndataview.py`

Show the existing views for the selected artifact: overview table, obs/var/X/layers/uns tabs,
optional heatmap. For MuData, a modality selector + per-modality summary. **Refactor first:**
extract `anndataview.py`'s inline stats into importable helpers (`_matrix_stats`,
`_layer_stats_rows`, `_format_uns` — already recommended in the June code review) so both the
standalone viewer and this tool call the same code.

## Reuse before duplicate (AGENTS.md)

- catalog → `raw_file_db_downloaded.csv` (already the index; `test_data.py` already reads it).
- software→slug and rule lookup → `rules.registry` / `params.registry` (don't reinvent aliases).
- convertibility / column-variant check → `converters.recognize.matches`.
- conversion → `converters.assemble.convert` + `rules.loader.load_packaged_rule`.
- MuData assembly → factor the helper out of `tests/test_mudata_levels.py`.
- summaries → shared helpers extracted from `anndataview.py`.

## Framework

Recommend **marimo** for the prototype: reactive (selecting target re-filters the left table,
selecting a right row re-renders the panel — no manual callbacks), tables with selection
(`mo.ui.table(selection="single")`), `mo.ui.dropdown` / `mo.ui.range_slider` / `mo.ui.run_button`,
and it reuses the existing marimo viewer directly. This is also a low-risk testbed for the
framework question in [TODO_viewer.md](TODO_viewer.md); if `APB_viewer` later goes Dash/Panel,
the conversion + summary helpers carry over unchanged (they're framework-agnostic).

## Phased build

1. **Catalog + filters** — left table from the CSV, software + size filters (read-only).
2. **Convertibility column** — registry-driven `✓`/targets per row; target selector re-filters.
3. **Convert (single AnnData)** — ion-level convert for one row → right table; disk cache.
4. **Summary panel** — wire right-row selection to the refactored `anndataview` helpers.
5. **MuData + fragment** — add the heavy targets behind the size guardrail + background run.

## Open questions

- **Filename:** captured as `TODO_ui_test_tool.md` (spec). Want me to also scaffold the marimo
  app as `src/anndata_proteomics/scripts/ui_test_tool.py`?
- **Artifacts:** in-memory vs disk cache (recommend disk)?
- **Convertibility check:** registry-only (fast) or also `recognize.matches` per row (one header
  read each — catches column variants but costs I/O on first load)?
- **Big-file policy:** hard cap, warn-only, or auto-offer per-run/subset for fragment/MuData?
- **Same app as `APB_viewer` or separate internal tool?** (Recommend: prototype here, fold in later.)
- How does `FragPipe (DIA-NN quant)` map — `fragpipe`, `diann`, or its own rule?
