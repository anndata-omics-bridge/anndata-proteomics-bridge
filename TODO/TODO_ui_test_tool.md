# TODO: test-data browser -> convert -> inspect GUI

**Date:** 2026-06-22
**Status:** Implemented prototype; this note tracks current behavior and remaining gaps.
**Related:** [TODO_viewer.md](TODO_viewer.md), [TODO_to_mu_data.md](TODO_to_mu_data.md),
[HOWTO/test_gui.md](HOWTO/test_gui.md).

## Current Flow

The APB GUI is `src/anndata_proteomics/scripts/ui_test_tool.py` and is launched with `make ui`.
It is a compact marimo app:

1. Browse/filter the downloaded ProteoBench catalog.
2. Select a target and dataset.
3. Start a background conversion subprocess through `jobrunner.py`.
4. Write durable outputs under `logs/ui_converted/<run>/`.
5. Scan converted outputs into a compact table.
6. Select a converted run and inspect it with `_ui_support.summarize()`.

The standalone `src/anndata_proteomics/scripts/anndataview.py` is **not** the active GUI result
path. It remains useful future work for a richer `.h5ad` viewer, but urgent parameter visibility
belongs in `_ui_support.summarize()`.

## Conversion Contract

The active conversion path is:

`ui_test_tool.py -> convert_one.py -> _ui_support.convert_target() -> assemble.convert(..., params_path=...)`

Key behavior:

- Rule selection is version-aware: the co-located param/log file is parsed for software version.
- Rule selection is level-semantic: a standalone AnnData level is exposed only when the vendor
  output has real quantitative layers for that level. Peptide and peptidoform identifiers are
  mandatory `.var` metadata/link columns where available or computable, but DIA-NN precursor
  layers are not aggregated into peptide/peptidoform matrices during parsing.
- Finished single-level conversions write `result.h5ad`.
- Finished MuData conversions write `result.h5mu`.
- Parsed search parameters are stored under
  `uns["anndata_proteomics"]["search_parameters"]` with `search_parameters_path`.
- MuData stores shared search parameters at `mdata.uns["anndata_proteomics"]` and keeps
  modality-level AnnData provenance.

## Converted-Output Table

The converted-output table is durable, not session-local. It scans `logs/ui_converted/`.

Visible columns:

- `run_name`
- `software_name`
- `software_version`
- `target`
- `status`
- `result_type`
- `nr_prec`
- `size_mb`

Full filesystem paths stay internal for selection/loading. For finished artifacts, table metadata
comes from the stored result artifact first; log/catalog reconstruction is fallback only.

## Result Inspector

The active result inspector renders `_ui_support.summarize(...)` as JSON:

- AnnData: shape, obs/var columns, layers, uns keys, X stats, full non-empty search parameters.
- MuData: n_obs, shared search parameters once at the MuData level, then per-modality summaries.
- When MuData-level search parameters are available, modality summaries suppress duplicate
  `search_parameters`.

The summary uses `params.anndata_io.read_search_parameters()` and
`get_search_parameters_path()`, not an ad hoc JSON key whitelist.

## Remaining Gaps

- Add a visible error/status column for parameter-parse failures instead of silently making rows
  non-convertible.
- Add big-file guardrails for fragment/MuData: hard cap, per-run subset, or row-capped conversion.
- Later, factor common matrix/layer/uns rendering helpers between `_ui_support.summarize()` and
  `anndataview.py` if the standalone viewer becomes important again.
- Decide whether `FragPipe (DIA-NN quant)` should be its own rule family or map through existing
  FragPipe/DIA-NN rules.
