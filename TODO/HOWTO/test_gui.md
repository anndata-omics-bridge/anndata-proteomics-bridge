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
            select row + Convert ──convert_target()──▶ AnnData | MuData ──summarize()──▶ panel
```

## `_ui_support.py` walkthrough

- **`software_slug(name)`** — maps a catalog `software_name` ("DIA-NN") to a parsing-rule vendor
  slug ("diann") with `re.sub(r"[^a-z0-9]", "", name.lower())`. This deliberately yields no
  match for vendors without a rule dir (e.g. "FragPipe (DIA-NN quant)" → `fragpipediannquant`),
  which is exactly the "not convertible" signal we want.
- **`convertible_levels(slug)`** — registry-driven, *not* a hardcoded list: tries
  `rules.registry.find_rule(slug, level)` for each of the five `LEVELS` and keeps the ones that
  resolve. Add a TOML for a new (vendor, level) and it shows up automatically.
- **`available_targets(slug)`** — the convertible levels, plus `"mudata"` only when **all five**
  levels resolve (today: DIA-NN only).
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
  - `"mudata"` → `_build_mudata`: convert all five levels and `MuData(mods, axis=0)`, **prefixing
    each level's `var_names`** (`frg:/ion:/pfm:/pep:/prt:`) so the modalities don't collide on the
    global axis. This mirrors `tests/test_mudata_levels.py` (a candidate to factor into one shared
    helper later).
- **`summarize(obj)`** — a render-ready dict: shape / obs+var columns / layers / uns / X stats for
  an AnnData, or `{n_obs, modalities: {name: <per-modality summary>}}` for a MuData (detected via
  `hasattr(obj, "mod")`).

## `ui_test_tool.py` — the marimo specifics (the tricky part)

A marimo notebook is a **reactive DAG of cells**: each `@app.cell` function's *arguments* are its
dependencies and its *return tuple* publishes names other cells consume. When a value changes,
only dependent cells re-run. Three things needed care:

1. **Accumulating results across clicks → `mo.state`.** marimo is functional, so session memory
   (the list of converted artifacts) lives in `mo.state([])`. The convert cell appends with the
   functional setter form `set_artifacts(lambda items: items + [record])` so it never *reads*
   the state it writes (which would create a feedback loop).

2. **"Run once per click", not on every selection change.** The Convert button is a counter:
   `mo.ui.button(value=0, on_click=lambda c: c + 1)`. The convert cell compares
   `convert_button.value` to a `get_last_click()` state value; it only acts when they differ, then
   records the new count. Changing the table selection re-runs the cell but the counter is
   unchanged, so nothing converts.

3. **Never `mo.stop()` in a cell the layout depends on.** The first version used
   `mo.stop(...)` to short-circuit before a click — but the layout cell depends on the convert
   cell's output, so on first load (no click) the cell stopped *before defining its output* and
   the whole UI failed to render. Fix: the side-effecting convert cell **returns nothing** and
   writes a `get_status/set_status` state string instead; a separate tiny cell turns that string
   into `convert_msg = mo.md(...)`. The status then **persists** (it's state) instead of flashing,
   and the layout always renders.

Cell graph (top to bottom):
`imports → title → catalog → controls(target/software/size) → state(+button) → filtered →
left_table → convert(side effects) → convert_msg → right_table → summary_panel → layout`.
The final **layout cell** composes everything with `mo.vstack` / `mo.hstack(widths=[1,1])` into
the two-pane view (catalog+convert on the left, results+summary on the right). UI elements are
created in their own cells and *arranged* in the layout cell — the standard marimo pattern.

The summary cell reads the right table's selection (`getattr(right_table, "value", None)` — the
right "table" is an `mo.md` placeholder until something is converted), maps the selected `#`
back into the artifacts list, and renders `summarize(...)` as JSON.

## Conversion logs ("Open log ↗")

Every Convert — success **or** failure — writes a per-run log file and the UI surfaces both a
new-tab link and the on-disk path.

- **Where:** `logs/ui_test_tool/<timestamp>_<slug>_<target>.log` under the repo root
  (`LOG_DIR` in `_ui_support.py`; `logs/` is gitignored). The absolute path is shown in the UI
  under the status line ("Log on server: …").
- **What's in it:** a header (time, dataset, absolute input file, slug, target), the per-step
  progress (each level's `obs × var` for a MuData build), then either `RESULT: OK` + the
  `summarize()` JSON, or `RESULT: FAILED` + the full traceback. This is how the
  `KeyError: 'PG.Normalised'` case reads as "ion/peptidoform/peptide OK → protein FAILED",
  pinpointing the missing column and the level.
- **How the new-tab link works:** `log_open_link_html()` embeds the log text in a
  `data:text/plain;base64,…` URL inside an `<a target="_blank">`. Deliberately **not** a
  `file://` link (browsers block those from an http page) and **not** `data:text/html` (Chrome
  blocks top-level navigation to it). text/plain opens reliably in a new tab without the marimo
  server needing to serve files. It's rendered with `mo.Html(...)` (raw HTML; `mo.md` would
  strip the `target` attribute).
- **Implementation:** `convert_with_log()` wraps `convert_target()`, threading a `log()`
  callback through `_convert_level`/`_build_mudata`, and always writes the file in a `finally`-
  style path before returning a `ConvertResult(ok, obj, log_path, summary, message)`.

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
  `available_targets` pick it up via the registry; no GUI change.
- **New summary field:** extend `summarize()`; the panel renders whatever it returns.
- **Editing the UI:** prefer `make ui-edit` (marimo's reactive editor surfaces DAG errors live).

## Validation performed

- `_ui_support` exercised against the real cache: 83 catalog rows; convertibility correct
  (DIA-NN → 5 levels + mudata, others → ion/peptidoform/none); filters correct; one real
  end-to-end convert + summarize.
- `ruff` clean on both modules; `gui` extra installs and imports (marimo 0.23.10, plotly 6.8.0);
  the marimo app imports (all cells register); `make help` works.
- **Not** click-tested: marimo's live reactive runtime can't be driven headlessly here, so the
  button → state → table wiring is validated structurally only. Run `make ui` to confirm.

## Known limitations / next steps

- Heavy conversions (`fragment`, `mudata` on big files) run **synchronously** and only show a
  size *warning* — no background execution or per-run/subset option yet (see the guardrail item
  in the spec).
- Converted artifacts are held **in memory** for the session; no disk cache / reload.
- `anndataview.py` and this tool duplicate summary logic — fold both onto the `_ui_support`
  helpers (the `_matrix_stats`/`_format_uns` extraction the June review flagged).
