# Acquisition Method as a Rule Selector

Date: 2026-07-25

> Capture the acquisition method (DDA / DIA) in the parsed search parameters, and let a
> parsing-rule level declare search-parameter-conditional overrides so DIA-NN DDA reports
> get `Ms1_Normalised` in `X` while DIA reports get `Precursor_Normalised`.

## Requirements

**User story.** The acquisition method decides *conversion behaviour*. DIA-NN writes the
same report format whether it analysed DIA or DDA runs — identical headers, both
`Precursor.Normalised` and `Ms1.Normalised` present — so nothing in the file itself says
which quantity belongs in `X`. Today `diann/v2/rules.json` hard-codes `Ms1_Normalised`
as the ion `x_layer`, which is the DDA choice applied to every v2 conversion, DIA
included. APB must pick the right quantity at convert time.

Once converted, nothing downstream re-decides: `apb.proteobench` scoring reads `X`,
`obs`, `var` and nothing else. So acquisition is consumed during rule resolution only,
and the effective `ParseRule` handed to the converters stays flat — the converters never
learn the concept exists.

**Scope (in).**

- New required `acquisition_method` field on `Parameters`, values restricted to
  `DDA` / `DIA` / `unknown`.
- Real detection in the DIA-NN parameter parser only. Every other parser reports
  `unknown` for now.
- A generic, declarative override mechanism in `rules.json`, keyed on parsed
  search-parameter values, resolved at rule-materialization time.
- `diann/v2` ion level: DIA default restored to `Precursor_Normalised`, with
  `Ms1_Normalised` as the DDA override. `diann/v1` stays DIA-only (DIA-NN had no DDA
  mode then).
- An explicit *applied column mapping* in the descriptive summary, so a converted object
  documents which vendor column landed in `X` and in each layer.

**Scope (out).**

- No CLI override (`--acquisition …`). The parameter file is the only source for now.
- No `DIA-PASEF` / `DDA-PASEF` vocabulary. Acquisition sub-flavours do not change
  conversion behaviour today; add a value only when a document needs to key an override
  on it.
- No `search_parameter_overrides` in `base` — level-only until something needs it.
- No detection sweep across the other vendors' parsers. Note that DIA-NN is the case we
  *know* writes one format for both modes; other tools may write a different format per
  mode (FragPipe), which is a separate rules document selected by columns/version, not
  an override.

**Acceptance.**

- Converting the cached DIA-NN DDA submission
  (`test_data_download/json_dir/Results_quant_ion_DDA_Astral/300beac4bd267751972cf484bb1cdee2fda0b3a4/`,
  `input_file.parquet` + `param_0..txt`, DIA-NN 2.6.0) puts `Ms1_Normalised` in `X`.
- Converting a DIA-NN DIA submission puts `Precursor_Normalised` in `X`.
- No parameter file, or `acquisition_method = "unknown"` → DIA default, i.e. the level's
  own `x_layer`, no override applied.
- `apb summary` on either output states which vendor column landed in `X`.

## Design

### `Parameters.acquisition_method`

```python
AcquisitionMethod = Literal["DDA", "DIA", "unknown"]
```

Required field, no `None`, no `Optional` handling at call sites — `"unknown"` is the
honest, typed value for a tool whose parameters do not state the mode. This follows the
existing finding in [TODO_params_model_review.md](TODO_params_model_review.md), *"required
fields must not be represented as optional nullable fields"*.

DIA-NN detection markers, confirmed against the cached DDA log
(`param_0..txt`, DIA-NN 2.6.0 Enterprise):

- command line contains `--dda`
- log line `All runs will be analysed as DDA runs`
- (corroborating, not needed) `WARNING: QuantUMS cannot be used on DDA data, disabled`

Everything else the DIA-NN parser sees stays `DIA`, which is DIA-NN's default acquisition
method.

### `search_parameter_overrides` in the rule level

The condition lives in the JSON as data, inside the level that owns the affected field —
no `levels` re-nesting:

```json
"levels": {
  "ion": {
    "axis": { "x_layer": "Precursor_Normalised" },
    "search_parameter_overrides": [
      { "when_search_parameters": { "acquisition_method": "DDA" },
        "axis": { "x_layer": "Ms1_Normalised" } }
    ]
  }
}
```

Both names say where the value comes from: APB already calls these *search parameters* in
`uns` (`search_parameters_path`, `search_parameters_version_status`,
`write_search_parameters`), so the rule schema uses the same word rather than a bare
`when`.

- `when_search_parameters` is **equality only**, on validated `Parameters` field names —
  a mapping, not an expression language.
- Merge order: `base` → level → every matching override, in source order, later wins.
  Reuses `_merge_rule_dicts`; the materialized result is still a flat `ParseRule`.
- An override body is a `RuleFragment`, so it can carry more than `axis` if a future case
  needs a different `missing_values` or a `required` layer.
- Generic by construction: any `Parameters` field can key an override. `acquisition_method`
  is simply the first one.

Selection stays layered: the `software_version` regex picks the *document*, the level key
picks the *level*, and the overrides then apply on top. Acquisition is a rule selector
like `software_version` — but expressed inside the document rather than by duplicating a
~130-line document to change one string.

### Applied column mapping

Conversion already writes a stage-owned `uns['anndata_proteomics']['descriptive_summary']`
with components (`quantification`, `annotation`, `fasta`, …) that `apb summary` prints. Add
a `column_mapping` component, derived from the effective rule so it cannot drift:

```json
"column_mapping": {
  "X": {"layer": "Ms1_Normalised", "source": "Ms1.Normalised"},
  "layers": {"Precursor_Normalised": "Precursor.Normalised", "Q_Value": "Q.Value"},
  "obs": {"Run": "Run"},
  "var": {"Precursor_Charge": "Precursor.Charge", "ProForma_ion": "computed:proforma_ion"}
}
```

This is the essential provenance — *what landed where* — not "which override fired". No
`applied_search_parameter_overrides` key; the mapping already shows the outcome, and
`rule_json` keeps the full effective rule for anyone who wants the detail.

### Alternatives set aside

- **Separate documents per acquisition** (`diann/v2/dda/rules.json`) mirroring `v1`/`v2`:
  duplicates a whole document to change one string, invites drift, and needs
  `_discovery.py` changes for the extra nesting level.
- **Acquisition-keyed `x_layer`** (`"x_layer": {"DIA": …, "DDA": …}`): minimal, but
  special-cases one field; the first DDA difference that is not `x_layer` forces a schema
  change.

## Implementation plan

- [ ] Add `AcquisitionMethod` and the required `acquisition_method` field to
      `params/model.py`; extend `_SERIES_FIELDS` round-tripping and the params CSV
      expectations in `tests/params/*.csv`.
- [ ] Detect `--dda` / `All runs will be analysed as DDA runs` in
      `params/parsers/diann.py`; default `DIA`. Set `unknown` in every other parser in
      `params/parsers/`.
- [ ] Copy the cached DDA log + a report subset into `tests/params/` /
      `tests/data/` as the DDA fixture (source: the `300beac4…` submission above).
- [ ] Add `search_parameter_overrides` (list of `{when_search_parameters, …RuleFragment}`)
      to `RuleFragment` in `rules/schema.py`; validate `when_search_parameters` keys
      against `Parameters.model_fields` and values against those fields' types.
- [ ] Extend `ParseRuleDocument.effective_rule(level, search_parameters=None)` to apply
      matching overrides after the base→level merge, keeping the return type `ParseRule`.
- [ ] Thread the resolved `Parameters` (already parsed in
      `converters/pipeline.py:resolve_parameters`) into rule materialization —
      `_select_rule` / `select_rule` / `convertible_levels` / `build_mudata` — so the
      converters receive an already-overridden flat rule.
- [ ] Flip `diann/v2/rules.json` ion `x_layer` to `Precursor_Normalised` and add the DDA
      override.
- [ ] Make `validate_rule_source` materialize and validate every override variant, so a
      mistyped override `x_layer` fails at `apb validate`, not at conversion time.
- [ ] Add the `column_mapping` component in `readers/summary.py`, derived from the
      effective rule; surface it in `apb summary`.
- [ ] Update `docs/json_schema.md` (override block, merge order, validation) and
      `docs/parameter_parsers.md` (`acquisition_method`); log in `CHANGES.md`.

**Files touched.** `params/model.py`, `params/parsers/diann.py` + the other
`params/parsers/*.py`, `rules/schema.py`, `rules/loader.py`, `converters/pipeline.py`,
`readers/summary.py`, `parsing_rules/diann/v2/rules.json`,
`parsing_rules/_schema/*.schema.json` (regenerated via `apb export-schema`),
`docs/json_schema.md`, `docs/parameter_parsers.md`, `CHANGES.md`.

**Test strategy.**

- Unit: DIA-NN parser returns `DDA` for the new fixture and `DIA` for the existing
  `tests/params/DIANN_*.log.txt`; every other parser returns `unknown`.
- Unit: override merge — `acquisition_method="DDA"` yields `x_layer="Ms1_Normalised"`,
  `"DIA"` and `"unknown"` and *no parameters* all yield `Precursor_Normalised`.
- Unit: `apb validate` rejects an override naming an undeclared layer and an unknown
  `when_search_parameters` key.
- Integration: convert the cached DDA Astral submission and a DIA submission; assert `X`
  provenance via the `column_mapping` summary component.
- Gate: `uv run pre-commit run --hook-stage pre-push --all-files`.

## Open questions

- Which DIA-NN version first supports `--dda`? `v1` (1.8.x) is assumed DIA-only; if a
  1.9.x document exists later, it may need the same override.
- Do other vendors have a usable acquisition marker in their parameter files, or does
  `unknown` stay permanent for them? Revisit when a second vendor needs an override.
