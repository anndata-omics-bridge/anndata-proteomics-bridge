# TODO / Research: APB_viewer — visualising large AnnData/MuData

**Date:** 2026-06-22
**Status:** Research / decision doc — no code yet.
**Question asked:** We have a marimo GUI for AnnData. We're worried about hitting a scale
limit and are considering a Dash app (`APB_viewer`) as a top-level project. **Is Dash the
right direction, and is Plotly enough or do we need more efficient viz libraries?**

---

## TL;DR

- **Don't frame it as "Dash vs marimo."** At our data sizes the binding constraint is the
  **rendering + data-access strategy**, not the app framework. Any framework — Dash, Panel,
  marimo, Streamlit — chokes the moment you push ~10⁵–10⁶ glyphs or a 10⁶-row table into a
  browser. Pick the *rendering engine* first, then the framework.
- **Plotly alone is not enough** for the big views. Plotly is excellent for *small/aggregated*
  charts (volcano, PCA, per-feature profiles) and we should keep it for those. For the large
  views it caps out: SVG traces die at ~10⁴ points, WebGL (`Scattergl`) at ~10⁵–2×10⁵.
- **The three things that actually scale** (all framework-independent):
  1. **Data access:** store/read AnnData/MuData as **Zarr**, read **lazily/backed** — never load
     a 10⁶-feature matrix into RAM.
  2. **Big matrices / dense scatter → Datashader** (server-side rasterisation to an image;
     handles 10⁷–10⁸ points).
  3. **Big tables (browse 827k-feature `var`) → Perspective** (FINOS; WebAssembly + Arrow,
     virtualised, streams only visible rows).
- **Is Dash "right"?** Dash is a *defensible, production-grade* choice for a deployable
  multi-user product — **but** its native big-data story is weaker than HoloViz's, so you'll
  hand-build the Datashader/Perspective plumbing. If you want the least plumbing for big
  scientific data, **Panel + HoloViews + Datashader** is the stronger technical fit. **marimo**
  stays great for the exploratory/dev viewer and is already in place.
- **Before building anything custom**, evaluate **Vitessce** — a web viewer that already does
  large **AnnData _and MuData_** via Zarr tiling. APB now emits MuData, so this is unusually
  well-aligned and might give the matrix/embedding views "for free."

---

## 1. Where we are today

`src/anndata_proteomics/scripts/anndataview.py` is a single-file **marimo** app: load one `.h5ad`, show overview +
obs/var/X/layers/uns tables + an optional Plotly `px.imshow` heatmap. It already **sidesteps
scale by truncating** everything:

- X preview: 10 × 20; heatmap: 50 × 50; tables via `mo.ui.table` (fine for obs, strained for var).

It works *because* it never shows the data at full size. The worry is correct: the moment a
view needs the real feature axis, the current approach breaks. (Plotly and marimo aren't even
pinned deps yet.)

## 2. The data regime (this is what makes the choice)

| Level | obs (samples/runs) | var (features) |
|---|---|---|
| protein | 6–100s | ~7k |
| peptide / peptidoform | 6–100s | ~60–66k |
| ion | 6–100s | ~73k |
| **fragment** | 6–100s | **~827k** (≈90 % dense) |
| **MuData (all 5)** | shared obs | **~1M+** features |

Shape is **wide and short**: a small obs axis, a *huge* var axis, dense-ish float matrices,
and now a **multi-modal MuData** with cross-level foreign-key links. Implications:

- A full **heatmap** (100 samples × 827k features ≈ 8×10⁷ cells) cannot be drawn as
  SVG/canvas — it must be **rasterised** (Datashader) or feature-selected/clustered first.
- A **volcano/MA/scatter** over 10⁵–10⁶ features needs WebGL or Datashader, not SVG.
- A **`var` table** of 10⁵–10⁶ rows needs server-side/virtualised paging, not a DOM table.
- The **cross-level linked navigation** (ion → peptidoform → peptide → protein via the MuData
  FK columns) is **APB's unique view** and exists in no off-the-shelf tool.

## 3. Reframe: three independent layers

Decide each layer separately; they compose.

### Layer A — Data access (highest leverage, framework-independent)
Write AnnData/MuData to **Zarr** (chunked, lazy/cloud-readable) and read **backed/lazy**
(`anndata` backed mode, `read_lazy`, dask) so the browser/app only ever materialises the
slice in view. This single change is what lets *any* of the frameworks below scale, and it's
the format the domain web-viewers (Vitessce, CELLxGENE) already expect.

### Layer B — Rendering engine (per view type)
| View | Small/aggregated | Large (10⁵–10⁸) |
|---|---|---|
| Heatmap / matrix | Plotly `imshow` | **Datashader** raster → image pane |
| Scatter (volcano, MA, PCA/UMAP) | Plotly `Scattergl` (≤~10⁵) | **Datashader** |
| Line / per-feature profile | Plotly | **plotly-resampler** (1-D/time-series only) |
| Table (browse obs/var) | `mo.ui.table` / Dash DataTable | **Perspective** (Arrow/WASM, virtualised) |

### Layer C — App framework
| Framework | Big-data fit | Deploy / multi-user | Notes |
|---|---|---|---|
| **Dash** | Plotly-native; WebGL ≤~2×10⁵; Datashader integration is **clunky** per community; resampler is 1-D only | **Strong** (WSGI scaling, auth, Enterprise) | Structured, callback-based; best for a real *product*; you build the big-data plumbing yourself |
| **Panel + HoloViews + Datashader** (HoloViz) | **First-class Datashader** (10⁷–10⁸, 2-D), Bokeh/WebGL | Good (deployable, Jupyter-native) | Least plumbing for big scientific data; the technical sweet spot here |
| **marimo** (current) | Reactive (only reruns dependents — no Streamlit rerun tax); embeds Plotly, and Perspective/Datashader via `anywidget` | OK (notebook *is* app) | Best for the exploratory/dev viewer; already in place |
| **Streamlit** | Full-rerun model; RAM grows per user; weak at scale | Session-affinity headaches | **Not recommended** at this scale |

## 4. Domain-specific prior art — evaluate before building

- **Vitessce** — web framework for **large AnnData / MuData / SpatialData via Zarr tiling**
  (browser pulls only visible tiles; WebGL/deck.gl). APB emits MuData, so the data model lines
  up unusually well. Config-driven (JSON view config). *Caveat:* it's single-cell/spatial-
  oriented (embeddings, spatial coords); our level-specific matrices + FK-linked navigation may
  not map cleanly to its built-in view types. But even if we don't adopt it wholesale, it sets
  the **target format (Zarr)** we'd want anyway and could cover the matrix/embedding views.
- **CZ CELLxGENE Discover** — scalable web viewer for huge AnnData (embeddings + expression,
  TileDB/Census backend, millions of cells). Engine reuse is hard (single-cell-shaped), but
  it's the reference point for "web viewer for huge AnnData."
- **Reality check:** neither does proteomics quant levels or APB's cross-level FK navigation.
  That **custom** view is the reason to build *something* of our own — but we can still lean on
  Zarr + Datashader/Perspective rather than reinventing rendering.

## 5. Recommendation

**Adopt a layered stack, not a single framework bet:**

1. **Now (cheap, decouples everything):** make APB write **Zarr** alongside `.h5ad`, and have
   the viewer read **lazily/backed**. Do this regardless of framework choice.
2. **Rendering:** keep **Plotly** for small/aggregated charts; add **Datashader** for
   matrices/large scatter and **Perspective** for large tables.
3. **Framework — two viable paths:**
   - **If `APB_viewer` is a real deployable, multi-user product:** **Dash** is a reasonable,
     mature pick — *provided* we budget for Datashader-as-image and Perspective plumbing. Its
     enterprise/deploy story is the best of the four.
   - **If we want the least big-data engineering:** **Panel + HoloViews + Datashader** is the
     stronger technical fit for wide/dense matrices and is also deployable. Strongly consider it
     as the actual recommendation if "handle the big views with minimal custom code" is the goal.
   - **Keep marimo** as the developer/exploratory viewer either way — it's already here, it's
     reactive, and it can embed the same engines.
4. **In parallel, spike Vitessce** on a real APB Zarr/MuData export. If it covers the
   matrix + embedding views acceptably, the custom app shrinks to just the **proteomics-specific
   cross-level navigation**, which is the only thing nobody else does.

**Direct answer to "is Dash right?"** Yes-ish: Dash is a sound choice for a *product*, but it is
**not** where the scaling comes from, and on its own it does **not** beat the limit you're
worried about. The scaling comes from Zarr + Datashader + Perspective. Choose Dash for
structure/deployment; choose **Panel/HoloViz** if you'd rather the framework hand you the
big-data rendering. **Don't** pick Streamlit. **Don't** expect Plotly alone to carry the large
views.

## 6. Suggested de-risking spikes (small, before committing)

- **S1 — Lazy Zarr read:** export one MuData to Zarr; open backed/lazy; confirm a single
  level's matrix slice loads without materialising the whole object. *(Validates Layer A.)*
- **S2 — Datashader heatmap:** rasterise a 100 × 827k fragment matrix to an image in <1 s and
  display it (Plotly image / Bokeh / Panel). *(Validates the matrix view at scale.)*
- **S3 — Perspective table:** load an 827k-row `var` table into a Perspective viewer; confirm
  instant scroll/sort/filter. *(Validates the table view at scale.)*
- **S4 — Vitessce config:** point Vitessce at the S1 Zarr; see how far a heatmap + embedding +
  obs/var view gets with zero custom code. *(Validates buy-vs-build.)*
- **S5 — Framework bake-off:** build the *same* one view (Datashader heatmap + a linked filter)
  in **Dash** and in **Panel**; compare lines-of-plumbing and responsiveness. *(Decides Layer C.)*

## 7. Open questions for the user (these change the recommendation)

- **What views does `APB_viewer` actually need?** (heatmap, volcano/MA, PCA/UMAP, per-feature
  profiles, QC/missingness, **cross-level FK navigation**, raw table browsing?) The view list,
  not the framework, should drive the choice.
- **Who is the audience / deployment model?** Single analyst on a laptop (marimo/Panel is plenty)
  vs. a hosted multi-user service (Dash/Panel server, auth) vs. shareable static export
  (Vitessce on object storage).
- **Is `APB_viewer` a separate top-level repo?** Reasonable — it consumes APB outputs, mirroring
  how `annProtSum` (R/Quarto reports) is already a separate sibling. Keep APB emitting the data
  contract (Zarr/h5ad + `uns` schema); keep the viewer downstream.
- **Buy vs build:** are we willing to adopt Vitessce/Zarr conventions, or is the cross-level
  proteomics navigation important enough to justify a custom app from day one?

## Vitessce — resources (downloads, repos, docs)

**Project & documentation**
- Website / live app: https://vitessce.io/
- Docs (getting started): https://vitessce.io/docs/
- Python docs: https://python-docs.vitessce.io/
- R docs: https://r-docs.vitessce.io/
- Examples gallery: https://vitessce.io/examples/
- Tutorials: https://vitessce.io/docs/tutorials/
- Data file types & formats (AnnData/MuData/SpatialData → Zarr): https://vitessce.io/docs/data-file-types/
- Data hosting (static server / object storage): https://vitessce.io/docs/data-hosting/
- View config (JSON): https://vitessce.io/docs/view-config-json/
- JavaScript/React API overview: https://vitessce.io/docs/js-overview/

**Repositories** (org: https://github.com/vitessce)
- Core JS/TS library: https://github.com/vitessce/vitessce
- Python API + Jupyter widget: https://github.com/vitessce/vitessce-python
- R API + htmlwidget: https://github.com/vitessce/vitessceR
- **easy_vitessce** (configure Vitessce from scverse/scanpy plotting APIs — most relevant for us): https://github.com/vitessce/easy_vitessce
- Data-loading utils (HuBMAP formats): https://github.com/vitessce/vitessce-data
- R analysis/conversion helpers: https://github.com/vitessce/vitessceAnalysisR
- Python tutorial (2026, latest): https://github.com/vitessce/vitessce-python-tutorial-2026
- Python tutorial (2023, HuBMAP): https://github.com/vitessce/vitessce-python-tutorial-2023
- Deploy-to-GitHub-Pages demo: https://github.com/vitessce/vitessce-demo-gh-pages
- Paper figures: https://github.com/vitessce/paper-figures

**Packages** (no conda-forge build as of 2026-06-22 — PyPI/npm only)
- PyPI (`pip install vitessce`, needs Python ≥3.9): https://pypi.org/project/vitessce/
- npm (`npm i vitessce`): https://www.npmjs.com/package/vitessce
- Python releases/changelog: https://github.com/vitessce/vitessce-python/releases

**Paper / citation**
- Vitessce, *Nature Methods* 2025: https://www.nature.com/articles/s41592-024-02436-x

## Sources

- [Dash: high-performance visualization (WebGL limits)](https://plotly.com/python/performance/)
- [Visualizing a billion points: Dash + plotly-resampler](https://medium.com/dbsql-sme-engineering/visualizing-a-billion-points-databricks-plotly-dash-and-the-plotly-resampler-45461bc3f466)
- [Plotting large datasets in Dash (Ploomber)](https://ploomber.io/blog/plotly-large-dataset/)
- [HoloViews: working with large data using Datashader](https://holoviews.org/user_guide/Large_Data.html)
- [Dash + HoloViews docs](https://dash.plotly.com/holoviews)
- [Datashader introduction](https://datashader.org/getting_started/Introduction.html)
- [Why Panel (vs Dash) for data apps](https://medium.com/@marcskovmadsen/i-prefer-to-use-panel-for-my-data-apps-here-is-why-1ff5d2b98e8f)
- [FINOS Perspective (WASM + Arrow, virtualised datagrid)](https://github.com/finos/perspective)
- [Vitessce: integrative visualization of multimodal single-cell data (Nature Methods 2025)](https://www.nature.com/articles/s41592-024-02436-x)
- [CZ CELLxGENE Discover: scalable single-cell data platform (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11701654/)
- [marimo vs Streamlit (reactive vs full-rerun)](https://marimo.io/features/vs-streamlit-alternative)
- [Streamlit vs Dash in 2025 (Squadbase)](https://www.squadbase.dev/en/blog/streamlit-vs-dash-in-2025-comparing-data-app-frameworks)
