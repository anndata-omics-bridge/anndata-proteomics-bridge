# TODO (handoff) — Column aliasing + required-column extraction (qpx-led, alphabase secondary)

Date: 2026-07-07
Status: HANDOFF — scope + references only. Design and implementation happen in a **separate session**.

## Task

Extend APB's **column contract** with two interlocking features:

1. **Column aliasing (review item 2).** Let a declared internal column accept **several candidate
   vendor names** (first-present-wins) instead of a strict 1:1 `Internal_Name = "Vendor col"`. This
   absorbs *minor* intra-version naming drift in one rule — complementary to, not a replacement for,
   APB's version-folder mechanism for *major* structural version differences.

2. **Required-column extraction.** Define the set of columns APB treats as **mandatory** and
   **errors** on when absent, using `qpx` (primary) and `alphabase` (secondary) as the reference for
   *which* columns are essential per level.

These are the **same schema surface** — "this internal column may come from any of these vendor
names, and it is mandatory (error if none present)" — so they are one design, not two.

Do **not** design or implement here. This note records scope + where to look.

## Starting point in APB

- APB **already** has a `required` flag — but **only on layers**:
  `apb/src/anndata_proteomics/rules/schema.py` (`required: bool`, `layer_required(...)`; the
  `x_layer` is always required). The task extends this notion to **identity / `var` / `obs`
  columns**, and adds the alias-list form to `[columns.*.select]`.
- The column-selection surface to change: `[columns.obs.select]` / `[columns.var.select]` in the
  parsing-rule TOMLs, the pydantic schema, and the loader/converter column resolution + `recognize`.

## References — qpx (primary)

qpx is itself a search-engine → MuData converter and already implements exactly this contract:

- **Per-field schemas:** `related_work/qpx/docs/spec/schemas/*.yaml` (`psm.yaml`, `pg.yaml`,
  `feature.yaml`, `run.yaml`, `sample.yaml`, …) — each field carries a `required` flag.
- **Validator:** `related_work/qpx/qpx/core/data/schema.py` — checks "Required columns are present"
  and emits `Missing required column: <name>`; field loading at
  `related_work/qpx/qpx/core/data/loader.py` (`required = fdef.get("required", False)`).
- Reconcile with the archived reviews: `apb/TODO/Archive/REVIEW_qpx.md`,
  `apb/TODO/Archive/REVIEW_qpx_vs_apb.md`.

## References — alphabase (secondary)

- **Alias lists (first-present-wins):** `column_mapping` values that are lists, e.g.
  `mobility: ['IM','IonMobility']`, in
  `related_work/alphabase/alphabase/constants/const_files/psm_reader.yaml`; resolution logic in
  `related_work/alphabase/alphabase/psm_reader/utils.py` (`get_column_mapping_for_df`) and
  `psm_reader.py` (`_get_actual_column`).
- **Implicit essentials:** alphabase's peptide-centric must-haves (`sequence`, `charge`,
  `raw_name`, `intensity`) — a sanity check for APB's required set, though APB's set is level-driven
  and richer.
- Its `pg_reader.yaml` `measurement_regex` groups (raw/lfq/ibaq/razor/unique/maxlfq) are a useful
  checklist that each measurement variant is captured as its own layer.

## Reference documents

- [REVIEW_what_we_learn_from_alphapep_base.md](../../REVIEW_what_we_learn_from_alphapep_base.md) — item 2 (alias lists) and the required-columns discussion.
- `related_work/qpx/docs/spec/schemas/` + `related_work/qpx/qpx/core/data/schema.py`, `loader.py`.
- `related_work/alphabase/alphabase/constants/const_files/psm_reader.yaml`, `psm_reader/utils.py`.
- APB: `src/anndata_proteomics/rules/schema.py` (existing `required`), `docs/toml_schema.md` (the `[columns.*.select]` contract).
