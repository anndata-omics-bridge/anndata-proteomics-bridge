# PLAN — MaxQuant, Sage, and AlphaPept ion coverage

> Implementation plan. Ground truth, diagnostics, and the evidence behind every decision live in
> [TODO_maxquant_sage_alphapept_dda_coverage.md](TODO_maxquant_sage_alphapept_dda_coverage.md) —
> this document is the file-by-file execution order and must not restate the analysis.

**Date:** 2026-07-30 · **Approved:** 2026-07-30 · **Repos touched:** `apb`, `apb_studio`.
**Status:** implemented and verified 2026-07-30 — see "Outcome" at the end.

**Goal:** the 9 MaxQuant, 1 Sage, and 2 AlphaPept submissions in `Results_quant_ion_DDA` report
`supported` with branches `('mudata', 'ion')`, plus the 3 MaxQuant submissions in `dda_astral`,
`dia_astral`, and `dia_singlecell`. 15 submissions total.

**Order rationale:** Phase 1 is the only schema change and Phase 2 depends on it. Phases 3 and 4 are
independent of both and of each other — Sage needs no code change at all and can be pulled forward
if MaxQuant stalls.

---

## Phase 1 — `optional_select` in the rule schema

New concept: a vendor column the rule captures **when the export has it** and skips when it does
not. `select` keeps its current meaning (required, gates recognition). Mirrors the existing
`Layer.required = false` semantics, so the vocabulary is one idea, not two.

### 1.1 `src/anndata_proteomics/rules/schema.py`

- `ColumnGroup`: add `optional_select: dict[str, str] = Field(default_factory=dict)`.
- `ColumnGroup.names`: include `optional_select` keys, after `select` and before `compute`, keeping
  the existing `dict.fromkeys` de-duplication.
- `ColumnGroup`: new `model_validator` rejecting a name present in both `select` and
  `optional_select`.
- `PartialColumnGroup`: add `optional_select: dict[str, str] | None = None` so a level can extend the
  base's optional set through `_merge_rule_dicts` (dict-on-dict deep merge already handles it).
- `ParseRule._computed_column_consistency`: seed `available_var_columns` from
  `columns.var.select | columns.var.optional_select`, so a compute may name an optional source.
- `ParseRule._derived_columns_are_not_selected`: extend `selected_sources` with both groups'
  `optional_select` values — an APB-derived column must not sneak in through the optional door
  either.
- `ParseRule._axis_keys_are_declared_columns`: **unchanged**. It reads `ColumnGroup.names`, which now
  includes optional names, so add an explicit check that no `axis.obs_keys` / `axis.var_keys` entry
  resolves to an `optional_select` name. An axis key is required by definition; a rule whose feature
  key can vanish is not a rule.
- New `model_validator`: reject `columns.obs.optional_select` when `input_shape == "wide"`.
  `converters/wide.py` accepts only the `<sample>` placeholder on the wide obs axis, so an optional
  wide obs column has no meaning. **`converters/wide.py` therefore needs no change.**

### 1.2 `src/anndata_proteomics/converters/recognize.py`

- `_expected_long_columns`: read only `columns.obs.select` / `columns.var.select`.
- `_required_var_columns`: read only `columns.var.select`.

Optional columns must never gate recognition — that is the whole point.

### 1.3 `src/anndata_proteomics/converters/assemble.py`

`_materialize_columns` runs before the long/wide dispatch, so var-side optional handling covers both
shapes from this one place.

- `_materialize_column_group`: after the `select` loop, iterate `optional_select`, skipping sources
  absent from the frame; return the set of skipped output names.
- `_materialize_columns`: collect the skipped names from both groups and thread them into the
  compute loop.
- `_compute_column`: take an `allow_missing: frozenset[str]` parameter. For `coalesce` and
  `join_nonempty`, drop `from_` entries that are in `allow_missing` before the existing
  missing-source check, then raise as today on anything still missing. Do **not** relax the check
  generally — a typo in a rule must still fail loudly. If every source of a coalesce was skipped,
  raise: that is a rule declaring a column it can never produce.
- `_columns_needed_for_long`: add `optional_select` values from both groups to `needed`. It already
  filters to columns present in the frame, so absent ones drop out naturally.

### 1.4 Schema artifacts and docs

- Regenerate both JSON Schemas: `uv run python -m anndata_proteomics.rules._export_schema`
  (rewrites `parsing_rules/_schema/parse_rule.schema.json` and `parse_rule_document.schema.json`).
- [docs/json_schema.md](../docs/json_schema.md): document `optional_select` beside the optional-layer
  paragraph, stating that `select` gates recognition and `optional_select` does not, that axis keys
  may not be optional, and that it is forbidden on the wide obs axis.

### 1.5 Tests (`tests/`)

1. A long level with an `optional_select` column converts identically with and without that column
   present, and the column appears in `var`/`obs` only in the first case.
2. A `coalesce` whose first source is a skipped optional select falls through to the next source.
3. A `coalesce` whose *every* source was skipped raises.
4. A required `select` source that is absent still raises `cannot select column …`.
5. A name in both `select` and `optional_select` is rejected at validation.
6. An `axis.var_keys` entry naming an `optional_select` column is rejected.
7. `columns.obs.optional_select` on a wide rule is rejected.
8. Recognition ignores optional columns: a rule matches a header set lacking them.

---

## Phase 2 — MaxQuant

`src/anndata_proteomics/parsing_rules/maxquant/rules.json` only. No code.

1. `software_version`: `^2\.6\.7\.0$` → `^(1\.5|1\.6|2\.)`.
2. `base.columns.obs`: `select` keeps `{"Raw_File": "Raw file"}`; move `Experiment` and `Fraction`
   into `optional_select`.
3. `levels.ion.columns.var`: move `Leading_Proteins` and `Leading_Razor_Protein` into
   `optional_select`, and add the 1.5.2.8 title-case spellings there as
   `"Leading_Proteins_Legacy": "Leading Proteins"` and
   `"Leading_Razor_Protein_Legacy": "Leading Razor Protein"`.
4. Extend the existing `Proteins` coalesce `from` to
   `["Proteins", "Leading_Proteins", "Leading_Proteins_Legacy"]`.

Integration test: convert one 1.x and one 2.x fixture, assert `obs_names` equal the six
`dda_qexactive.toml` `raw_file` values and that `ProForma_ion` var names carry `UNIMOD` accessions.
Assert 2.6.7.0 still converts unchanged — it is the only currently-working MaxQuant fixture and the
one thing this phase can regress.

---

## Phase 3 — Sage

New `src/anndata_proteomics/parsing_rules/sage/rules.json`. No code. Single `ion` level, wide.

- `software_name` `"Sage"`, `software_version` `"^0\\.15\\."`.
- `base`: `input_shape` `wide`; `axis.obs_keys` `["sample"]`; `axis.duplicates.mode` `"error"`;
  `columns.obs.select` `{"sample": "<sample>"}`.
- `base.modifications`: `token_regex` on `peptide`, `token_pattern` `"\\[([^\\]]+)\\]"`,
  `token_position` `after_residue`, `unknown_policy` `preserve`; map `+57.021465 → UNIMOD:4`,
  `+15.994915 → UNIMOD:35`, `+42.010567 → UNIMOD:1`.
- `levels.ion.axis`: `var_keys` `["ProForma_ion"]`, `x_layer` `"Intensity"`.
- `levels.ion.columns.var.select`: `Peptide ← peptide`, `Charge ← charge`, `Proteins ← proteins`,
  `Q_Value ← q_value`, `Score ← score`, `Spectral_Angle ← spectral_angle`.
- `levels.ion.columns.var.compute`: `ProForma_peptidoform` and `ProForma_peptide` from `Peptide`,
  `ProForma_ion` from `[ProForma_peptidoform, Charge]`.
- `levels.ion.column_roles.protein_accessions`: `"Proteins"`.
- `levels.ion.layers`: one entry, `Intensity`,
  `source` `"^(?P<sample>.+)\\.(?:mzML|mzml|mzML\\.gz|raw|d)$"`, `missing_values` `[0]`.

Integration test: assert `shape == (6, 86878)`, the six `obs_names`, and — guarding the PEAKS failure
mode recorded in [TODO_fragpipe_peaks_version_coverage.md](TODO_fragpipe_peaks_version_coverage.md) —
that no non-run column (`peptide`, `charge`, `proteins`, `q_value`, `score`, `spectral_angle`) is
captured as a sample.

**Superseded by Phase 5.** The `^0\.15\.` pin shipped first as a proxy for
`lfq_settings.combine_charge_states == false`, then was replaced by a real parameter gate at `^0\.`
once the upstream default (`true`) showed the proxy would mis-classify the common case.

---

## Phase 4 — AlphaPept

### 4.1 `apb`: new `src/anndata_proteomics/parsing_rules/alphapept/rules.json`

Single `ion` level, long. `software_version` `"^0\\.5\\."`.

- `base.axis`: `obs_keys` `["Raw_File"]`, `duplicates.mode` `"keep_first"`.
- `base.columns.obs.select`: `Raw_File ← shortname`, `Sample_Group ← sample_group`,
  `File_Path ← filename`.
- `base.modifications`: `token_regex` on `sequence`, `token_pattern` `"[a-z]+"`, `token_position`
  `before_residue`, `case_sensitive` `true`, `unknown_policy` `preserve`; map `c → UNIMOD:4`,
  `ox → UNIMOD:35`, `a → UNIMOD:1`.
- `levels.ion.axis`: `var_keys` `["ProForma_ion"]`, `x_layer` `"Intensity"`.
- `levels.ion.columns.var.select`: `Sequence ← sequence`, `Sequence_Naked ← sequence_naked`,
  `Charge ← charge`, `Precursor ← precursor`, `Protein ← protein`,
  `Protein_Group ← protein_group`, `Decoy ← decoy`.
- `levels.ion.columns.var.compute`: the three ProForma computes, from `Sequence`.
- `levels.ion.column_roles.protein_accessions`: `"Protein_Group"`.
- `levels.ion.layers`: `Intensity ← ms1_int_sum_apex_dn` (the x_layer), then `MS1_Int_Sum_Apex`,
  `MS1_Int_Sum_Area`, `MS1_Int_Max_Apex`, `MS1_Int_Max_Area`, `Score ← score`, `Q_Value ← q_value`,
  `Retention_Time ← rt`, `FWHM ← fwhm`, `N_Fragments_Matched ← n_fragments_matched`. All non-x
  layers stay optional, which is what lets one document serve 0.5.0 (85 columns) and 0.5.3 (81).

### 4.2 `apb_studio`: delete the duplicate header reader

`src/apb_studio/capabilities.py`:

- Delete `read_table_headers` (lines ~208–224) and its `pandas` / `pyarrow.parquet` imports if they
  become unused.
- Import `read_table_columns` from `anndata_proteomics.readers.dispatch` and call it at line ~113,
  wrapping the result in `tuple(...)` to keep `headers` hashable.

Test updates:

- `tests/test_capabilities.py:445` — monkeypatch `capabilities.read_table_columns` instead of
  `capabilities.read_table_headers`. The assertion (the header read must not happen for a vendor
  with no rule document) is unchanged.
- `tests/test_remaining_quality.py:156` — `read_table_columns` raises `UnknownFormat` (a `ValueError`
  subclass) with `unsupported extension '.xlsx' for …`, so change the target call and the
  `pytest.raises` match string from `unsupported table extension` to `unsupported extension`.
- New regression test: a comma-delimited `.txt` fixture yields its full header list, not one column.

Integration test in `apb`: convert both AlphaPept fixtures, assert shapes `(6, 49318)` and
`(6, 58668)` and that `Intensity` is sourced from `ms1_int_sum_apex_dn`.

---

## Phase 5 — Sage level gating (implemented 2026-07-30)

Verifying the Sage version pin against Sage's own docs and source turned this from a nice-to-have
into a defect fix. `DOCS.md` gives `lfq_settings.combine_charge_states` **default `true`**, and
`crates/sage-cli/src/runner.rs` writes a combined row as `charge.unwrap_or(-1)`. So the *default*
Sage configuration produces a peptidoform-level `lfq.tsv`, and the `^0\.15\.` pin was a proxy for a
parameter whose default is the opposite of what the proxy assumed. Reproduced: a charge-collapsed
table presented as version `0.15.0` resolved to the ion rule and died with
`ValueError: charge must be positive, got -1` — a `FAILED` cell, the very thing the pin existed to
prevent. It only worked because the one cached 0.15.0-beta.2 submission sets `false` explicitly.

1. `params/model.py` — `Parameters.combine_charge_states: bool | None`.
2. `params/parsers/sage.py` — read `quant.lfq_settings.combine_charge_states`, defaulting to Sage's
   own `true` when the LFQ block omits the key, and `None` when there is no LFQ block at all (no
   quantification, so no charge-state decision to report).
3. `rules/schema.py` — `LevelRuleFragment.requires_search_parameters` plus `is_available()`, and
   `ParseRuleDocument.level_is_available()` / `available_levels()`. Excluded from `as_merge_dict`,
   so the gate cannot leak into the merged `ParseRule` body. A gated level is **unavailable when no
   parameters could be parsed** — guessing would pick the wrong level for exactly the inputs the
   gate exists to separate.
4. `rules/loader.py` `resolve_rule_locator` and `converters/pipeline.py`
   `_column_matching_rule_variants` — both apply the gate. The second is not redundant: it iterates
   packaged rules directly, bypassing the resolver, and two levels of one document can share both
   the version regex and the header schema.
5. `parsing_rules/sage/rules.json` — rewritten: `software_version` `^0\.`, shared body in `base`,
   `ion` gated on `combine_charge_states: false` and `peptidoform` on `true`. The peptidoform level
   does not select `charge`, so the `-1` sentinel never reaches the output.

Result: **Sage 1/2 → 2/2.** 0.15.0-beta.2 → `('mudata', 'ion')` at `(6, 86878)`; 0.14.6 →
`('mudata', 'peptidoform')` at `(6, 31200)`. Corpus-wide 76 supported, 22 unsupported, 0 failed.
Verified that `ion` is unreachable for a combined export at every version tried (`0.14.6`, `0.15.0`,
`0.16.0`, `None`) and that missing parameters yield no level rather than a wrong one.

Fallout: level count 15 → 16; `test_converters_e2e` gained `_level_available_for`, because
`find_test_data`'s "first cached file for this vendor+version" heuristic cannot pick between two
parameter-gated levels — it now parses the sibling parameter file exactly as production does.

## Still out of scope

- `Layer.sources` first-present-wins, which would recover `MS_MS_Count` for MaxQuant 1.5.2.8.
- i2MassChroQ (5), MSAngel, ProlineStudio, quantms — 8 submissions, each needing a rule document and
  a parameter parser.

---

## Verification

Run after each phase, not only at the end:

1. `uv run pre-commit run --hook-stage pre-push --all-files` in `apb`, then in `apb_studio`.
2. Capability sweep over `dda_qexactive` via `apb_studio.capabilities.discover_capabilities`:
   MaxQuant 9, Sage 1, AlphaPept 2 report `supported` with `('mudata', 'ion')`; i2MassChroQ,
   MSAngel, ProlineStudio and quantms keep their existing "no packaged parsing-rule document"
   diagnostic verbatim.
3. No regression across the corpus: MaxQuant 2.6.7.0 / 2.6.3.0 / 2.7.5.0, FragPipe, PEAKS, WOMBAT,
   DIA-NN, Spectronaut, AlphaDIA.
4. `obs_names` of every newly converting fixture match its module annotation `raw_file` values, so
   the `annotate` and `proteobench` stages can run on the new branches.
5. Dated entry in `apb/CHANGES.md` and `apb_studio/CHANGES.md`.

---

## Outcome (2026-07-30)

`Results_quant_ion_DDA` went from 20 `UNSUPPORTED` submissions to 8 — the four out-of-scope
vendors. Corpus-wide: **75 supported, 23 unsupported, 0 failed** over 98 fixtures. All 12 MaxQuant,
both AlphaPept, and the charge-resolved Sage submission convert; DIA-NN (28), FragPipe (5),
FragPipe/DIA-NN quant (10), PEAKS (7), Spectronaut (8) and WOMBAT (2) are unchanged. Both repos pass
`pre-commit --hook-stage pre-push --all-files`, and `apb convert` writes a real `.h5mu` for one
fixture per vendor.

Two things differed from the plan as written:

- **`converters/_axis.py` also needed the change.** `build_axis_frame` is called with the rule's
  *declared* column set (`ColumnGroup.names`), so a skipped optional column raised `KeyError` in the
  axis frame builder even though materialization had correctly skipped it. It now filters to present
  columns, with a docstring recording why a declared-but-absent column can only be a skipped
  optional select. `converters/wide.py` needed no change, as planned.
- **`_compute_column` was split.** Adding the "every source skipped" branch pushed it past the C901
  complexity gate, so the `coalesce` / `join_nonempty` path moved into `_compute_generic_column`.

Test-suite fallout, all resolved: packaged-document counts moved 7 → 9 and levels 13 → 15 (five
tests); `test_data._PROTEOBENCH_PARAM_FIXTURES` gained the already-present `Sage` and `AlphaPept`
fixtures to keep the vendor-parity invariant; the MaxQuant coalesce assertion moved to
`optional_select`; and `find_test_data` became version-aware, which is what stops the end-to-end test
from handing the Sage ion rule the charge-collapsed 0.14.6 file. New tests:
`tests/test_optional_select.py` (11) and `tests/test_maxquant_sage_alphapept_coverage.py` (29).
