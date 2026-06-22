# HOWTO: the APB test-data browser GUI (`make ui`)

How the test-data browser → convert → inspect GUI is built, so it can be maintained and
extended. Spec/rationale live in [../TODO_ui_test_tool.md](../TODO_ui_test_tool.md);
the viz-scaling background is in [../TODO_viewer.md](../TODO_viewer.md).

## What it does

Browse the ProteoBench test corpus, filter by **target / software / size**, select a dataset,
**Convert** it to an AnnData level or a MuData, and inspect the result in a summary panel.
Picking a target re-filters the catalog to what's actually convertible (e.g. *MuData* → DIA-NN
only; *ion* → DIA-NN/MaxQuant/Spectronaut/PEAKS/FragPipe).

Launch: `make ui` (the Makefile installs the `gui` extra on first run via `uv`).

## File map

| File | Role |
|---|---|
| `Makefile` (repo root) | `make ui` / `ui-edit` / `viewer` / `test` / `lint` |
| `pyproject.toml` | added the `gui` optional-dependency extra (`marimo`, `plotly`) |
| `src/anndata_proteomics/scripts/_ui_support.py` | **all non-UI logic** — catalog, convertibility, conversion, summaries (no marimo import; unit-testable) |
| `src/anndata_proteomics/scripts/ui_test_tool.py` | the **marimo app** — thin reactive glue over `_ui_support` |
| `src/anndata_proteomics/scripts/anndataview.py` | the pre-existing single-`.h5ad` viewer (moved here from repo root) |

## Architecture: split testable logic from UI glue

The single most important decision: **marimo never appears in `_ui_support.py`.** Everything
that has behaviour — reading the index, deciding convertibility, running a conversion,
summarising — is a plain function there, so it can be imported and tested without a browser.
`ui_test_tool.py` only wires widgets to those functions.

Why: a marimo app can't be unit-tested headlessly (it needs a live reactive runtime + browser),
so anything important must live outside it. This also lines up with the reuse rule — the helpers
call the existing `rules.registry`, `rules.loader`, `readers.dispatch`, `converters.assemble`,
and `test_data` rather than reinventing them.

```
raw_file_db_downloaded.csv ──load_catalog()──▶ catalog DataFrame
                                                  │  (+ slug, targets, size_mb)
target/software/size widgets ──filter_catalog()──▶ left table
            select row + Convert ──convert_one subprocess──▶ logs/ui_converted/<run>/
                                                              │ result.h5ad | result.h5mu
list_converted_runs() ◀───────────────────────────────────────┘
            converted row selection ──load_converted_result() + summarize()──▶ inspector
```

## `_ui_support.py` walkthrough

- **`software_slug(name)`** — maps a catalog `software_name` ("DIA-NN") to a parsing-rule vendor
  slug ("diann") with `re.sub(r"[^a-z0-9]", "", name.lower())`. This deliberately yields no
  match for vendors without a rule dir (e.g. "FragPipe (DIA-NN quant)" → `fragpipediannquant`),
  which is exactly the "not convertible" signal we want.
- **`convertible_levels(slug)`** — registry-driven, *not* a hardcoded list: tries the known
  quantification levels and keeps only rules that resolve for the software version and match the
  file headers. A level TOML is valid only when the vendor output has real layers for that level;
  metadata/link columns such as `ProForma_peptide` do not create a standalone level.
- **`available_targets(slug)`** — the convertible levels, plus `"mudata"` when at least two
  report-backed levels resolve.
- **`load_catalog()`** — reads `test_data.DOWNLOADED_DB` (reusing that module's path constants),
  keeps `status == "ok"`, and adds `size_mb`, `slug`, `targets` (tuple), `targets_str`. Targets
  are cached per slug so `find_rule` isn't re-hit for every row. Returns an empty frame (not an
  error) when the gitignored cache index is absent.
- **`filter_catalog(catalog, target=, software=, max_size_mb=)`** — the three GUI filters; the
  target filter checks membership in each row's `targets` tuple.
- **`is_heavy(target, size_mb)`** — `True` for `fragment`/`mudata` above 100 MB; drives the UI
  warning (full-report fragment/MuData is the multi-GB case from the perf work).
- **`convert_target(input_file_path, slug, target)`** — `read_table` →
  - a level → `convert(df, load_packaged_rule(slug, level))`;
  - `"mudata"` → `_build_mudata`: convert the report-backed levels that resolve and
    `MuData(mods, axis=0)`, **prefixing each level's `var_names`** (`frg:/ion:/prt:` for current
    DIA-NN v1) so modalities do not collide on the global axis.
- **`list_converted_runs()`** — scans `logs/ui_converted/` for conversion run folders, detects
  `result.h5ad`, `result.h5mu`, and `console.log`. Finished results read software/version and
  parameter provenance from the stored artifact first; the conversion log/catalog are only
  fallbacks.
- **`converted_runs_table()`** — trims the internal scan result to user-facing columns (`run_name`,
  software, version, target, status, result type, size/precursor counts). Full paths stay internal
  for selection and loading.
- **`load_converted_result(path)`** — loads a selected `result.h5ad` or `result.h5mu`; unsupported
  suffixes fail clearly.
- **`summarize(obj)`** — a render-ready dict: shape / obs+var columns / layers / uns / X stats for
  an AnnData, or `{n_obs, modalities: {name: <per-modality summary>}}` for a MuData (detected via
  `hasattr(obj, "mod")`). When search parameters are present in `uns`, it surfaces the full
  non-empty parsed parameter payload plus `search_parameters_path`.

## `ui_test_tool.py` — the marimo specifics (the tricky part)

A marimo notebook is a **reactive DAG of cells**: each `@app.cell` function's *arguments* are its
dependencies and its *return tuple* publishes names other cells consume. When a value changes,
only dependent cells re-run. Three things needed care:

1. **Converted results are durable, not session state.** Each conversion writes a run directory
   under `logs/ui_converted/`; the right-side table is rebuilt by scanning those folders. Restarting
   marimo does not lose completed conversions.

2. **"Run once per click", not on every selection change.** The Convert control is an
   `mo.ui.run_button`, so changing the left table selection does not start a job. The run cell only
   launches a subprocess when the button value changes and no conversion is already running.

3. **The background job owns conversion state.** The app stores only the current/last
   `jobrunner.Job`, a static-server handle, and the finished `run_key`. The subprocess writes
   `console.log` and `result.*`; the GUI polls and renders snapshots without holding converted
   AnnData/MuData objects in memory.

Cell graph (top to bottom):
`imports → title → catalog → controls(target/software/size) → state(+button) → filtered →
left_table → convert/status → converted_runs → converted_table → result_viewer → layout`.
The final **layout cell** stacks the input dataset table above the converted-output table and
summary inspector, keeping each table to a small page size so the workflow stays compact.

The result-viewer cell maps the selected visible `run_name` back to the full internal run row,
loads `result.h5ad` or `result.h5mu`, and renders `summarize(...)` as JSON. For MuData, shared
search parameters are shown once at the MuData level and suppressed from modality summaries.
Log-only or incomplete folders stay visible in the table but render a concise "no result file"
message.

## Conversion outputs

Every Convert writes a per-run output directory and the UI surfaces the live status, output folder,
zip download, console log, and durable converted-output table row.

- **Where:** `logs/ui_converted/<timestamp>_<slug>_<target>/` under the repo root. `logs/` is
  gitignored.
- **Files:** success writes `result.h5ad` for a single level or `result.h5mu` for MuData, plus
  `console.log`; failures usually have only `console.log` and remain visible as incomplete rows.
  New conversions pass the selected parameter file through to `assemble.convert(...)`, so
  searchable parameter metadata is stored under `uns["anndata_proteomics"]["search_parameters"]`.
  MuData writes the same shared search-parameter container at `mdata.uns["anndata_proteomics"]`
  and keeps the modality-level AnnData provenance.
- **Subprocess:** `ui_test_tool.py` launches `python -m anndata_proteomics.scripts.convert_one`.
  `jobrunner.py` tails the log, tracks process state, zips output directories, and starts a small
  static HTTP server so browser links open reliably.

## Makefile & packaging

- A new **`gui`** extra (`marimo`, `plotly`) in `pyproject.toml`; they were used but never
  declared. GUI targets run `uv run --extra gui marimo run …`, so `uv` provisions them on first
  use — no manual install step.
- `make` defaults to `help`, which self-documents by `grep`-ing the `## ` comments on targets.
- The moved viewer is launchable via `make viewer FILE=path.h5ad`.

## Run / test / extend

```bash
make ui                       # launch the browser GUI
make ui-edit                  # marimo editor (live-edit the app)
make viewer FILE=foo.h5ad     # the single-file AnnData viewer

# test the logic (no GUI needed):
uv run --extra gui python -c "from anndata_proteomics.scripts import _ui_support as u; \
print(u.available_targets('diann'), len(u.load_catalog()))"
```

- **New convertible vendor/level:** just add the parsing-rule TOML. `convertible_levels` /
  `available_targets` pick it up via the registry; no GUI change. Only add a level TOML when the
  vendor file has real quantitative layers for that level. Capture identifiers and links as
  `.var` metadata, not as derived abundance matrices.
- **New summary field:** extend `summarize()`; the panel renders whatever it returns.
- **Editing the UI:** prefer `make ui-edit` (marimo's reactive editor surfaces DAG errors live).

## Validation performed

- `_ui_support` exercised against the real cache: 83 catalog rows; convertibility correct for
  report-backed DIA-NN levels plus MuData; filters correct; one real end-to-end convert +
  summarize.
- `ruff` clean on both modules; `gui` extra installs and imports (marimo 0.23.10, plotly 6.8.0);
  the marimo app imports (all cells register); `make help` works.
- **Not** click-tested: marimo's live reactive runtime can't be driven headlessly here, so the
  button → state → table wiring is validated structurally only. Run `make ui` to confirm.

## Known limitations / next steps

- Heavy conversions (`fragment`, `mudata` on big files) run in the background, but there is still
  no hard cap or per-run/subset option yet (see the guardrail item in the spec).
- The right-side viewer is a lightweight JSON summary/inspector, not the full scalable viewer from
  `TODO_viewer.md`.
- `anndataview.py` is deferred for the urgent parameter-provenance path. The active GUI result
  viewer is `_ui_support.summarize()`; fold both onto shared helpers when building the larger
  viewer.
