# MaxQuant, Sage, and AlphaPept ion coverage for the DDA modules

> Make the 12 cached MaxQuant, 2 cached Sage, and 2 cached AlphaPept submissions convert, not
> merely resolve: allow a rule to declare configuration-dependent vendor columns as optional, add
> two new parsing-rule documents, and stop `apb_studio` from reading comma-delimited `.txt` as TSV.

**Date:** 2026-07-30 · **Source:** APB capability discovery plus direct conversion runs over the
current `test_data_download` cache. **Owner:** apb (rule schema, recognizer, column materialization,
parsing rules) and apb_studio (capability header read).
**Status:** implemented and verified 2026-07-30. Execution record:
[PLAN_maxquant_sage_alpha.md](PLAN_maxquant_sage_alpha.md).

---

## Verified ground truth

`Results_quant_ion_DDA` (module `dda_qexactive`) holds 26 submissions. FragPipe (3), PEAKS (1) and
WOMBAT (1) convert today; the remaining 20 are `UNSUPPORTED`. Reproduced with
`apb_studio.capabilities.discover_capabilities` over every fixture in the module:

| Software | Cached | APB diagnostic today |
|---|---:|---|
| MaxQuant | 9 | `No APB parsing rule matches software 'maxquant', version '<v>', and the input headers.` |
| i2MassChroQ | 5 | `APB has no packaged parsing-rule document for software 'i2masschroq'` |
| AlphaPept | 2 | `APB has no packaged parsing-rule document for software 'alphapept'` |
| Sage | 1 | `APB has no packaged parsing-rule document for software 'sage'` |
| MSAngel | 1 | `APB has no packaged parsing-rule document for software 'msangel'` |
| ProlineStudio | 1 | `APB has no packaged parsing-rule document for software 'prolinestudio'` |
| quantms | 1 | `APB has no packaged parsing-rule document for software 'quantms'` |

This plan covers **MaxQuant, Sage, and AlphaPept** (12 submissions in this module, 15 across the
corpus). i2MassChroQ, MSAngel, ProlineStudio and quantms are out of scope and stay unsupported.

Parameter parsing is **not** the gap: `params/parsers/{maxquant,sage,alphapept}.py` are all
registered and already return a version (`0.15.0-beta.2`, `0.5.0`, `0.5.3`) plus enzyme and FDR.
The gap is the result-table side only.

### MaxQuant — one over-narrow version regex and two configuration-dependent columns

`parsing_rules/maxquant/rules.json` pins `software_version` to `^2\.6\.7\.0$`, which is the only
cached MaxQuant version whose `evidence.txt` carries all 12 columns the rule declares as required.
Diffing the rule's required set against every cached MaxQuant file:

| Module | Version | Missing required columns |
|---|---|---|
| dda_astral | 2.6.7.0 | — |
| dda_qexactive | 1.5.2.8 | `Experiment`, `Fraction`, `Leading proteins`, `Leading razor protein` |
| dda_qexactive | 1.5.3.30 · 1.5.8.2 · 1.6.3.3 · 2.1.3.0 · 2.1.4.0 · 2.3.1.0 · 2.4.13.0 · 2.5.1.0 | `Fraction` |
| dia_astral | 2.6.3.0 | `Experiment`, `Fraction` |
| dia_singlecell | 2.7.5.0 | `Fraction` |

Three distinct causes, and only the first is about versions:

1. **`Fraction` is absent from 11 of 12 files.** MaxQuant emits it only for fractionated designs.
2. **`Experiment` is absent from 2 files** — and those two straddle the version axis (1.5.2.8 and
   2.6.3.0, while 2.6.7.0 has it). It is design-dependent, not version-dependent, so no version
   split can express it.
3. **1.5.2.8 uses title-case `Leading Proteins` / `Leading Razor Protein` / `MS/MS Count`** where
   1.5.3.30 and later use `Leading proteins` / `Leading razor protein` / `MS/MS count`.

Because (1) and (2) are configuration-driven, the fix is not another version family — the schema
needs to express "capture this vendor column when the export has it". Optional **layers** already
work exactly this way (`Layer.required = false`); `columns.*.select` has no equivalent and
`recognize._expected_long_columns` therefore treats every selected column as mandatory.

**Verified:** a scratch document with the version regex broadened and those columns removed from the
required set matches all 12 cached MaxQuant files, and `converters.assemble.convert` produces
correct output for both extremes — 1.5.2.8 → `(6, 10538)` and 2.5.1.0 → `(6, 10566)` on 40k-row
truncations, with `obs_names` equal to the `dda_qexactive.toml` `raw_file` values and `ProForma_ion`
var names such as `[UNIMOD:1]-AAAAAAAAAAGAAGGR/2`. 1.5.2.8 yields layers
`Intensity, PEP, Retention_Time, Score`; 2.5.1.0 additionally yields `MS_MS_Count`.

### Sage — no rule document; `lfq.tsv` is a wide matrix

Both cached Sage submissions export the same 12-column `lfq.tsv`:
`peptide, charge, proteins, q_value, score, spectral_angle` plus one column per run named
`<run>.mzML`. `peptide` carries ProForma-shaped mass deltas already
(`LGMLSPEGTC[+57.021465]K`, N-terminal `[+42.010567]-M[+15.994915]ESQQ…`).

`peptide` + `charge` is unique across all 86 878 rows of the DDA submission, so
`duplicates.mode = "error"` is safe.

**The two submissions are not the same quantification level.** `lfq_settings.combine_charge_states`
in the Sage parameter JSON differs:

| Module | Version | `combine_charge_states` | `charge` column | Actual level |
|---|---|---|---|---|
| dda_qexactive | 0.15.0-beta.2 | `false` | 2 / 3 / 4 | ion |
| dda_astral | 0.14.6 | `true` | `-1` for all 31 200 rows | peptidoform |

`-1` is Sage's sentinel for charge-collapsed quantification. Converting the Astral file with an ion
rule fails in `_format_charge` with `charge must be positive, got -1` — correctly, because that
table has no ion axis.

**Verified:** a scratch Sage ion document converts the 0.15.0-beta.2 DDA submission to
`(6, 86878)`, `x_layer = Intensity`, `obs_names` exactly the six `dda_qexactive.toml` run names, and
var names `EAGELKPEEEITVGPVQK/3`, `LGMLSPEGTC[UNIMOD:4]K/2`. No code change was needed — the
existing `token_regex` parser maps all three mass deltas to `UNIMOD:4/35/1` via the bundled Unimod
registry, including the N-terminal form.

### AlphaPept — no rule document, and a delimiter bug in apb_studio

Both cached AlphaPept submissions ship `input_file.txt` that is **comma**-delimited: 85 columns
(0.5.0) and 81 columns (0.5.3), one row per PSM. 0.5.3's header is a strict subset of 0.5.0's — the
extras are `matching_p`, `mz_calib`, `rt_calib`, `type` — so one document covers both.

- `shortname` already equals the annotation `raw_file` values, so it is the obs key; `filename` is
  a Windows absolute path and is metadata only.
- Modified sequences use a lowercase token *before* the residue (`FESDTDTEcCIAK`, `HTGPIToxMLQFNPK`,
  N-terminal `aAGGKAGK`). Tokens observed: `c`, `ox`, `a`.
- `charge` is a float (`2.0`); `_format_charge` already normalizes it to `/2`.
- 795 of ~151 000 `(run, precursor)` pairs repeat (up to 5 rows), so a duplicate policy is required.

**Which intensity is the x_layer is settled empirically, not guessed.** Correlating each `ms1_int*`
column against ProteoBench's own per-run intensities in the fixture's `result_performance.csv`
(run `…Condition_A_Sample_Alpha_01`, 17 209 shared precursors):

| Candidate column | Fraction matching ProteoBench exactly |
|---|---:|
| `ms1_int_sum_apex_dn` | **0.980** |
| `ms1_int_sum_apex` | 0.000 |
| `ms1_int_max_apex` | 0.000 |
| `ms1_int_sum_area` | 0.000 |
| `ms1_int_max_area` | 0.000 |

`ms1_int_sum_apex_dn` is the x_layer. Re-running only that column under each duplicate policy:
`keep_first` 0.9948, `max` 0.9936, `sum` (APB's `aggregate`) 0.9797. **`keep_first`** both matches
ProteoBench best and is the honest semantics — it does not invent a sum across PSMs of one
precursor.

**146 rows are decoys** (`decoy = True`, `protein = REV__…`, `sequence = ELLPELR_decoy`). With
`unknown_policy = "preserve"` these render as `ELLPELR[decoy]` with `decoy` recorded in
`unknown_mod_tokens` and a correct stripped sequence — distinct from the target peptidoform and not
corrupting it. APB has no row-filter construct and does not filter MaxQuant `Reverse` rows either,
so this is left as-is and **not** a reason to add one.

**Verified:** a scratch AlphaPept ion document converts 0.5.0 to `(6, 49318)` and 0.5.3 to
`(6, 58668)` with 10 layers, `obs_names` equal to the annotation run names, and var names
`VILGGLK/2`, `FESDTDTEC[UNIMOD:4]IAK/2`.

**But a correct rule alone is not enough.** `apb_studio.capabilities.read_table_headers` hardcodes
`.txt → "\t"`:

```python
separators = {".csv": ",", ".tsv": "\t", ".txt": "\t"}
```

It returns **1 column** for AlphaPept's comma-delimited `.txt`, so header matching fails and the
fixture stays `UNSUPPORTED` regardless of the rule. This function duplicates
`anndata_proteomics.readers.dispatch.read_table_columns`, which handles `.parquet`, `.csv`, `.tsv`
and content-detects `.txt` via `readers.tabular.detect_text_delimiter` — and which the conversion
path already uses. This is the same class of reader bug recorded in
[TODO_fragpipe_peaks_version_coverage.md](TODO_fragpipe_peaks_version_coverage.md) for a
comma-delimited PEAKS `.txt`; Studio kept a private copy and so kept the bug.

---

## Plan

### Phase 1 — `optional_select` in the rule schema (apb)

Add a sibling of `select` for vendor columns that a rule should capture when present and skip when
absent, mirroring `Layer.required = false`. `select` keeps its current meaning: required.

1. `rules/schema.py`
   - `ColumnGroup.optional_select: dict[str, str] = {}`; include it in `ColumnGroup.names`.
   - `PartialColumnGroup.optional_select: dict[str, str] | None = None` so a level can extend the
     base's optional set through the existing deep merge.
   - `_computed_column_consistency`: treat optional-select names as available compute sources.
   - `_derived_columns_are_not_selected`: apply the same guard to `optional_select`.
   - Reject a name declared in both `select` and `optional_select` in the same group.
   - An `axis` key must not resolve to an `optional_select` name (`_axis_keys_are_declared_columns`
     stays restricted to `select` + `compute`) — an axis key is by definition required.
2. `converters/recognize.py` — `_expected_long_columns` and `_required_var_columns` read only
   `select`, so optional columns never gate recognition.
3. `converters/assemble.py`
   - `_materialize_column_group` — materialize `optional_select` after `select`, skipping absent
     sources; collect the skipped names.
   - `_columns_needed_for_long` — add `optional_select` sources (it already keeps only columns
     present in the frame).
   - `_compute_column` — for `coalesce` / `join_nonempty`, tolerate a `from` source only when it is
     a skipped optional select. Pass that set in explicitly; do **not** relax the existing
     missing-source error, which must still fire for a genuinely broken rule.
4. `converters/wide.py` — same optional handling in the obs select loop (line ~176).
5. Regenerate `parsing_rules/_schema/parse_rule_document.schema.json` via `rules/_export_schema.py`
   and document `optional_select` in [docs/json_schema.md](../docs/json_schema.md) next to the
   optional-layer paragraph.

Tests: a level declaring an optional select converts both with and without the column; a coalesce
over a skipped optional source still resolves; a required `select` column that is absent still
raises; `select`/`optional_select` name collision is rejected.

### Phase 2 — MaxQuant (apb)

Edit `parsing_rules/maxquant/rules.json` only:

- `software_version` → `^(1\.5|1\.6|2\.)` (covers 1.5.2.8 … 2.7.5.0; verified to still match
  2.6.7.0, so dda_astral does not regress).
- `base.columns.obs`: keep `select = {Raw_File: "Raw file"}`; move `Experiment` and `Fraction` to
  `optional_select`.
- `levels.ion.columns.var`: move `Leading_Proteins` / `Leading_Razor_Protein` to `optional_select`
  and add the 1.5.2.8 title-case spellings `Leading Proteins` / `Leading Razor Protein` as
  `Leading_Proteins_Legacy` / `Leading_Razor_Protein_Legacy`.
- Extend the existing `Proteins` coalesce to
  `[Proteins, Leading_Proteins, Leading_Proteins_Legacy]`.

Expected result: all 12 cached MaxQuant submissions convert, across four modules.

### Phase 3 — Sage (apb)

New `parsing_rules/sage/rules.json`: wide, single `ion` level, `software_version` `^0\.15\.`.

- obs `select = {sample: "<sample>"}`, `duplicates.mode = "error"`.
- var: `Peptide, Charge, Proteins, Q_Value, Score, Spectral_Angle`; computes
  `ProForma_peptidoform` / `ProForma_peptide` / `ProForma_ion`; `column_roles.protein_accessions =
  "Proteins"`.
- modifications: `token_regex` on `peptide`, `\[([^\]]+)\]`, `after_residue`, mapping
  `+57.021465 → UNIMOD:4`, `+15.994915 → UNIMOD:35`, `+42.010567 → UNIMOD:1`.
- one layer `Intensity`, `source = "^(?P<sample>.+)\.(?:mzML|mzml|mzML\.gz|raw|d)$"`,
  `missing_values = [0]`. Verified not to capture `peptide`/`charge`/`proteins`/`q_value`/`score`/
  `spectral_angle` as samples — the failure mode recorded for PEAKS in the FragPipe/PEAKS TODO.

The `^0\.15\.` pin is deliberately narrow and is a **proxy**: it excludes the 0.14.6 Astral
submission not because of its version but because `combine_charge_states = true` makes it
peptidoform-level. Recorded as Phase 5 rather than silently widened, so no fixture gets a rule that
will fail in `_format_charge`.

### Phase 4 — AlphaPept (apb + apb_studio)

1. New `parsing_rules/alphapept/rules.json`: long, single `ion` level, `software_version` `^0\.5\.`.
   - obs `select = {Raw_File: "shortname", Sample_Group: "sample_group", File_Path: "filename"}`,
     `duplicates.mode = "keep_first"`.
   - var: `Sequence, Sequence_Naked, Charge, Precursor, Protein, Protein_Group, Decoy` plus the
     three ProForma computes; `column_roles.protein_accessions = "Protein_Group"`.
   - modifications: `token_regex` on `sequence`, `[a-z]+`, `before_residue`,
     `case_sensitive = true`, mapping `c → UNIMOD:4`, `ox → UNIMOD:35`, `a → UNIMOD:1`.
   - layers: `Intensity ← ms1_int_sum_apex_dn` (x_layer) plus `ms1_int_sum_apex`,
     `ms1_int_sum_area`, `ms1_int_max_apex`, `ms1_int_max_area`, `score`, `q_value`, `rt`, `fwhm`,
     `n_fragments_matched`. The four 0.5.0-only columns stay unreferenced (or optional) so one
     document serves both versions.
2. **apb_studio:** delete `capabilities.read_table_headers` and call
   `anndata_proteomics.readers.dispatch.read_table_columns` instead, so the capability probe and the
   conversion path agree on every delimiter. Update `tests/test_capabilities.py:445` and
   `tests/test_remaining_quality.py:156`, which patch/assert on the private helper — the
   unsupported-extension assertion still holds, as `read_table_columns` raises
   `readers.dispatch.UnknownFormat` for `.xlsx`.

### Phase 5 — follow-ups (not in this change)

- ~~**Sage charge-collapsed exports.**~~ **Done 2026-07-30.** Verified against Sage's `DOCS.md`
  (`combine_charge_states` default **`true`**) and `crates/sage-cli/src/runner.rs`
  (`charge.unwrap_or(-1)`): the default configuration is peptidoform-level, so the interim
  `^0\.15\.` pin was a proxy for a parameter whose default is the opposite of the assumption, and a
  default-configured 0.15.x submission would have *failed* rather than been excluded. Replaced by
  `Parameters.combine_charge_states` plus `levels.<level>.requires_search_parameters`, a real
  availability gate. Sage is now 2/2.
- **Alternate layer source spellings.** MaxQuant 1.5.2.8 silently loses `MS_MS_Count` because its
  column is `MS/MS Count`. A first-present-wins `Layer.sources` list would recover it; layer names
  must stay unique, so it cannot be expressed today.
- **Remaining DDA vendors.** i2MassChroQ (5), MSAngel, ProlineStudio, quantms — 8 submissions, each
  needing a new rule document and a parameter parser.

---

## Verification

1. `uv run pre-commit run --hook-stage pre-push --all-files` in `apb` and `apb_studio`.
2. Re-run the capability sweep over `dda_qexactive`: MaxQuant 9, Sage 1, AlphaPept 2 report
   `supported` with branches `('mudata', 'ion')`; i2MassChroQ, MSAngel, ProlineStudio and quantms
   keep their existing "no packaged parsing-rule document" diagnostic.
3. Convert one fixture per vendor and check `obs_names` against the module annotation `raw_file`
   values, so the annotate and proteobench stages can run.
4. Confirm no regression on the four modules already converting MaxQuant 2.6.7.0 / 2.6.3.0 /
   2.7.5.0, plus FragPipe, PEAKS, WOMBAT, DIA-NN and Spectronaut.
5. `CHANGES.md` entry in both `apb` and `apb_studio`.
