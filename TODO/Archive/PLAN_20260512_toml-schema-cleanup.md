# PLAN 2026-05-12 — TOML schema & docs cleanup

## Background

Two-agent cross-review of `docs/toml_schema.md` and the six packaged
`parsing_rules/<vendor>/parse_*.toml` files surfaced five prioritised fixes.
Both reviewers (doc-first and TOML-first) independently converged on the
same top-5 list. This plan executes those five fixes. No new public API.

## Fix 1 (HIGH, docs) — Document `[modifications]` + `[[modifications.map]]`

The `[modifications]` block appears in every shipped TOML and is the most
edited section per vendor, yet is undocumented. A new contributor must
reverse-engineer it from `src/anndata_proteomics/rules/schema.py:88-112`.

**Action.** Add a new section `## Modifications` to `docs/toml_schema.md`
between Wide and "Software Families Already Seen". Cover:

- `source_column` (string, required) — vendor column whose values contain
  embedded modification tokens (e.g. `Modified.Sequence`).
- `parser` (default `"token_regex"`) — one of `token_regex`,
  `already_proforma`, `separate_mod_column`. Document the
  parser-consistency rules from `schema.py:99-112` (which fields are
  required / forbidden per parser).
- `token_pattern` (required for `token_regex`) — regex whose first
  capture group is the vendor token, e.g. `"\\(([^()]*)\\)"`.
- `token_position` (default `"after_residue"`) — one of
  `before_residue`, `after_residue`, `n_term`, `c_term`, `embedded`,
  `unknown`.
- `case_sensitive` (default `false`).
- `unknown_policy` (default `"preserve"`) — one of `preserve`, `drop`,
  `error`.
- `sequence_column` (optional) — only used when `parser =
  "separate_mod_column"`.
- `output_column` (default `"proforma_sequence"`) — the name of the
  derived column that `how = "proforma_sequence"` exposes.
- `[[modifications.map]]` (list of `{ token, accession }`) — vendor
  token → UNIMOD/MOD accession. Canonical name / target / position /
  mass_delta are filled at rule-load time from
  `src/anndata_proteomics/modifications/unimod_registry.toml`; do not
  duplicate them per-tool.

Add a worked example using FragPipe's actual block (numeric mass-delta
tokens) and a second example for DIA-NN (`UniMod:N` style tokens).

## Fix 2 (HIGH, TOMLs) — De-leak dataset-specific `column_pattern` regexes

`parse_fragpipe_ion_1.toml:46` uses `^(?P<sample>LFQ_.+) Intensity$`;
`parse_peaks_ion_1.toml:46` uses `^(?P<sample>LFQ.+) Normalized Area$`;
`parse_wombat_peptidoform_1.toml:35` uses `^abundance_(?P<sample>[AB]_[123])$`.
The `LFQ_`, `LFQ`, `[AB]_[123]` fragments encode the *current ProteoBench
sample naming* into a *vendor parsing rule*. Copy-pasting these to a new
dataset silently matches zero columns.

**Action.**

- FragPipe: change `(?P<sample>LFQ_.+)` → `(?P<sample>.+)` on every layer.
- PEAKS: change `(?P<sample>LFQ.+)` → `(?P<sample>.+)` on every layer.
- WOMBAT: change `^abundance_(?P<sample>[AB]_[123])$` →
  `^abundance_(?P<sample>.+)$` and same for `number_of_psms_`.
- Run pytest. If a now-broader pattern picks up an unintended column,
  narrow only the *suffix* (the vendor-fixed part of the column name),
  never the sample-token prefix.
- In `docs/toml_schema.md`, add a Wide note: "the `(?P<sample>...)`
  group must describe the *vendor column shape*, not the *user's
  sample naming*. Use `.+` for the sample token; filter or rename via
  `[sample_name_cleanup]` if needed."

## Fix 3 (MED-HIGH, TOMLs + docs) — Resolve WOMBAT `ProForma` overload

In 5 of 6 TOMLs `ProForma` is the ion-level compute
`how = "proforma_ion"` (`<peptidoform>/<charge>`). In WOMBAT the same
name is used for the sequence-level compute `how = "proforma_sequence"`.
Two different concepts under one name across the corpus is a footgun.

**Action.**

- In `parse_wombat_peptidoform_1.toml`:
  - Rename `[[columns.var.compute]] name = "ProForma"` →
    `name = "Peptidoform"` (the compute is `how = "proforma_sequence"`,
    which is exactly what `Peptidoform` means in every other TOML).
  - Update `axis.var_keys = ["ProForma"]` → `["Peptidoform"]`.
- In `docs/toml_schema.md`, add a "Naming convention" note:
  - `Peptidoform` → result of `how = "proforma_sequence"`.
  - `ProForma` → result of `how = "proforma_ion"` (ion-level only).
  - Peptidoform-level rules use `Peptidoform` as the `var_keys` entry;
    ion-level rules use `ProForma`.

## Fix 4 (MED, docs) — Document `obs_keys` convention

Long rules currently use four different obs-keys (`Run`, `Raw_File`,
`R_FileName`, plus `sample` for wide). There is no documented
convention, so a 7th-vendor author will invent yet another spelling.
Renaming all six is invasive and breaks downstream `adata.obs` access;
documentation is the lower-risk path.

**Action.** Add a "Obs-axis conventions" subsection to
`docs/toml_schema.md`:

- Wide rules: `obs_keys = ["sample"]`, derived from the
  `(?P<sample>...)` capture group.
- Long rules: `obs_keys` is the single column that uniquely identifies a
  run within the vendor's export. Preserve the vendor's natural name on
  the LHS (DIA-NN `Run`, MaxQuant `Raw_File`, Spectronaut `R_FileName`)
  — these are vendor-specific identifiers, not synonyms for "sample".
- Additional obs-side annotations (`Experiment`, `Fraction`,
  `R_Condition`, …) may appear in `[columns.obs.select]` even when not
  in `axis.obs_keys`; they enrich `adata.obs` but do not participate in
  uniqueness.

No TOML edits in this fix — it is purely documentation of existing
practice.

## Fix 5 (MED, docs) — General doc cleanup pass

Targeted edits to `docs/toml_schema.md`:

1. Title and framing — replace `# TODO: AnnData Mapping For Quant Parsing
   Rules` with `# TOML Schema Reference for Parsing Rules`. Rewrite the
   Purpose paragraph: this is a contract reference, not a planning
   document. Drop the "this is a TODO meant to become source material…"
   paragraph and the `## Decisions` heading (keep the content under a
   `## Background concepts` heading).
2. **Remove** the "tool-specific annotation file" sentence
   (line 96-97) — no such mechanism exists in `schema.py` or any TOML.
3. **Fix** the `columns.obs.compute` claim at line 156: schema
   (`schema.py:208-210`) forbids computed obs columns; `axis.obs_keys`
   must reference `[columns.obs.select]` only.
4. **Align** the long example's modifications block (line 286-290):
   shipped DIA-NN uses `parser = "token_regex"` with a `map`, not
   `parser = "already_proforma"`. Use the real shipped block.
5. **Clarify** reserved literal `stripped_sequence`: it is both the
   value of `how` and the conventional output name of that compute. The
   exact derived name written into `adata.var.columns` is whatever the
   compute's `name = "..."` is.
6. **State** that `Layer.encoding_mode = "factor"` requires a non-empty
   `categories` table (enforced by `schema.py:62-67`).
7. **Hoist** `encoding_mode` / `categories` documentation out of the
   per-shape Long / Wide sections into the common section — they apply
   to both.
8. **State** that `[sample_name_cleanup]` is wide-only (forbidden on
   long by `schema.py:162-165`).
9. **Clarify** filename casing: the folder and filename use lowercase
   vendor name (`diann/parse_diann_…`) while `software_name` inside the
   TOML preserves the canonical spelling (`"DIA-NN"`). The
   filename-vs-`quantification_level` match is what the packaged-rules
   test enforces (`tests/test_packaged_rules.py`).
10. **Drop** `uns` from the orientation bullet list, or move it to a
    "what the parser writes (not the rule)" note — TOML rules cannot
    declare `uns` entries.

## Execution order

1. Write this plan to `TODO/PLAN_20260512_toml-schema-cleanup.md`. ✓
2. Apply Fix 2 (TOML regex de-leak — three files).
3. Apply Fix 3 (WOMBAT rename — one file).
4. Run pytest; verify all green.
5. Apply Fixes 1, 4, 5 (docs only — single file edit).
6. Commit as one or two commits (TOML edits → docs).

## Out of scope (deferred)

- Stripping default-valued keys from shipped TOMLs (`mode = "error"`,
  default mod-block flags). Reviewer A defended the explicit defaults as
  intentional documentation; reviewer B wanted them stripped. Leave as
  is until we add a "minimal vendor template" doc.
- `Stripped_Sequence_Normalized` orphan compute (declared in every TOML
  but referenced nowhere in `var_keys` / downstream computes). Decide
  separately whether to remove or surface.
- Renaming long-rule obs-keys to a canonical name. Documentation-only
  for now.
- Porting remaining 5 vendor parameter parsers (separate task).
