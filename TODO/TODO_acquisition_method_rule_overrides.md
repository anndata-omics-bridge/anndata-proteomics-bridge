# Acquisition Method as a Rule Selector

Date: 2026-07-25
Status: Implemented

> Capture the acquisition method (DDA / DIA) in the parsed search parameters, and let a
> parsing-rule level declare search-parameter-conditional overrides so DIA-NN DDA reports
> get `Ms1_Normalised` in `X` while DIA reports get `Precursor_Normalised`.

## Requirements

**User story.** The acquisition method decides *conversion behaviour*. DIA-NN writes the
same report format whether it analysed DIA or DDA runs — identical headers, both
`Precursor.Normalised` and `Ms1.Normalised` present — so nothing in the file itself says
which quantity belongs in `X`. Before this change, `diann/v2/rules.json` hard-coded
`Ms1_Normalised` as the ion `x_layer`, which applied the DDA choice to every v2
conversion, DIA included. APB must pick the right quantity at convert time.

Once converted, nothing downstream re-decides the quantitative source:
`apb.proteobench` scores `X` (while also reading the stored effective rule for semantic
column roles). Acquisition is consumed during rule resolution only, and the effective
`ParseRule` handed to converters stays flat — the converters never learn the concept
exists.

**Scope (in).**

- New non-nullable `acquisition_method` field on `Parameters`, values restricted to
  `DDA` / `DIA` / `unknown`, with a backward-compatible `"unknown"` default.
- Real detection in the DIA-NN parameter parser only. Every other parser reports
  `unknown` for now.
- A declarative, level-only axis override mechanism in `rules.json`, keyed on parsed
  search-parameter values and resolved at rule-materialization time. The first schema
  deliberately overrides only axis fields; list-valued rule fragments need separate,
  identity-aware merge semantics before they can safely become overrideable.
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
- No layer/column/modification overrides. APB's document merger concatenates object arrays;
  pretending that it can update an existing named layer would create duplicate declarations.
- No detection sweep across the other vendors' parsers. Note that DIA-NN is the case we
  *know* writes one format for both modes; other tools may write a different format per
  mode (FragPipe), which is a separate rules document selected by columns/version, not
  an override.

**Acceptance.**

- Converting the cached DIA-NN DDA submission
  (`test_data_download/json_dir/Results_quant_ion_DDA_Astral/300beac4bd267751972cf484bb1cdee2fda0b3a4/`,
  `input_file.parquet` + `param_0..txt`, DIA-NN 2.6.0) puts `Ms1_Normalised` in `X`.
- Converting a DIA-NN DIA submission puts `Precursor_Normalised` in `X`.
- `acquisition_method = "unknown"` → DIA default, i.e. the level's own `x_layer`, no
  override applied.
- No parameter file → DIA default in the direct materialization API and explicit
  `--rule-config` flow. Ordinary packaged `apb convert` still requires `--params` to select
  the software-version document.
- `apb summary` on either output states which vendor column landed in `X`.

## Design

### `Parameters.acquisition_method`

```python
AcquisitionMethod = Literal["DDA", "DIA", "unknown"]
```

The field is non-nullable but defaults to `"unknown"`:

```python
acquisition_method: AcquisitionMethod = "unknown"
```

This keeps call sites free of `None` handling while preserving compatibility with stored
AnnData search-parameter payloads and CSV fixtures created before the field existed. Series
deserialization must handle this field before the generic missing-string normalization:
`"unknown"` is intentionally a real enum value here even though legacy parameter fields treat
that token as missing.

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
- Merge order: `base` → level → every matching axis override, in source order, later wins.
  Reuses `_merge_rule_dicts` only for dictionaries/scalars; the materialized result is still
  a flat `ParseRule`.
- Keep the schema non-recursive and enforce level-only placement with three types:
  `RuleFragment` (ordinary body), `SearchParameterOverride` (condition plus partial axis),
  and `LevelRuleFragment` (ordinary body plus the override list). `base` remains a plain
  `RuleFragment`.
- Generic by construction: any `Parameters` field can key an override. `acquisition_method`
  is simply the first one. Condition values are validated and normalized against the
  corresponding Pydantic field type before equality matching.

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
  "X": {"layer": "Ms1_Normalised", "source": "Ms1.Normalised", "source_kind": "column"},
  "layers": {
    "Ms1_Normalised": {"source": "Ms1.Normalised", "source_kind": "column"},
    "Precursor_Normalised": {"source": "Precursor.Normalised", "source_kind": "column"},
    "Q_Value": {"source": "Q.Value", "source_kind": "column"}
  },
  "obs": {"Run": "Run"},
  "var": {"Precursor_Charge": "Precursor.Charge", "ProForma_ion": "computed:proforma_ion"}
}
```

`layers` lists every materialized layer, including the layer also referenced by `X`.
Long-format sources have `source_kind="column"`; wide-format sources have
`source_kind="pattern"` because the declared source is a sample-capturing regex rather than
one literal vendor column. This is the essential provenance — *what landed where* — not
"which override fired". No `applied_search_parameter_overrides` key; the mapping already
shows the outcome, and `rule_json` keeps the full effective rule for anyone who wants the
detail.

### Alternatives set aside

- **Separate documents per acquisition** (`diann/v2/dda/rules.json`) mirroring `v1`/`v2`:
  duplicates a whole document to change one string, invites drift, and needs
  `_discovery.py` changes for the extra nesting level.
- **Acquisition-keyed `x_layer`** (`"x_layer": {"DIA": …, "DDA": …}`): minimal, but
  special-cases one field; the first DDA difference that is not `x_layer` forces a schema
  change.

## Implementation plan

- [x] Add `AcquisitionMethod` and
      `acquisition_method: AcquisitionMethod = "unknown"` to `params/model.py`; make
      `from_series` preserve the meaningful `"unknown"` value and test legacy payloads/CSVs
      that omit the field. `_SERIES_FIELDS` already derives itself from `model_fields`.
- [x] Detect `--dda` / `All runs will be analysed as DDA runs` in
      `params/parsers/diann.py`; default `DIA`. Set `unknown` in every other parser in
      `params/parsers/`.
- [x] Add a compact, source-faithful DDA log under `tests/params/`; exercise the full
      `300beac4…` report directly from the canonical cache in the integration test.
- [x] Add non-recursive `SearchParameterOverride` and `LevelRuleFragment` models in
      `rules/schema.py`. Override bodies contain a partial `axis`; validate
      `when_search_parameters` keys against `Parameters.model_fields` and normalize values
      against those fields' types.
- [x] Extend `ParseRuleDocument.effective_rule(level, search_parameters=None)` to apply
      matching overrides after the base→level merge, keeping the return type `ParseRule`.
- [x] Thread the resolved `Parameters` (already parsed in
      `converters/pipeline.py:resolve_parameters`) into rule materialization —
      `_select_rule` / `select_rule` / `convertible_levels` / `build_mudata` — so the
      converters receive an already-overridden flat rule. A library call that supplies
      `params_path` but no `ParameterResolution` parses it once before selection.
- [x] Apply the same materialization in both explicit `--rule-config` branches and pass
      parameters through the CLI's `convertible_levels` call.
- [x] Flip `diann/v2/rules.json` ion `x_layer` to `Precursor_Normalised` and add the DDA
      override.
- [x] Make `validate_rule_source` materialize the default and every reachable compatible
      override combination, so a mistyped override `x_layer` fails at `apb validate`, not
      at conversion time.
- [x] Add the `column_mapping` component in `readers/summary.py`, derived from the
      effective rule; surface it in `apb summary`.
- [x] Update `docs/json_schema.md` (override block, merge order, validation) and
      `docs/parameter_parsers.md` (`acquisition_method`); log in `CHANGES.md`.

**Files touched.** `params/model.py`, `params/parsers/diann.py` + the other
`params/parsers/*.py`, `rules/schema.py`, `rules/loader.py`, `converters/pipeline.py`,
`scripts/cli.py`, `readers/summary.py`, `parsing_rules/diann/v2/rules.json`,
`parsing_rules/_schema/*.schema.json` (regenerated via `apb export-schema`),
`docs/json_schema.md`, `docs/parameter_parsers.md`, `CHANGES.md`.

**Test strategy.**

- Unit: DIA-NN parser returns `DDA` for the new fixture and `DIA` for the existing
  `tests/params/DIANN_*.log.txt`; every other parser returns `unknown`.
- Unit: override merge — `acquisition_method="DDA"` yields `x_layer="Ms1_Normalised"`,
  `"DIA"` and `"unknown"` and *no parameters* all yield `Precursor_Normalised`.
- Unit: `apb validate` rejects an override naming an undeclared layer and an unknown
  `when_search_parameters` key.
- Unit: packaged and explicit-rule-config single-level/MuData paths all receive the
  materialized rule.
- Integration: convert the cached DDA Astral submission and the cached DIA Astral
  submission; assert both `column_mapping` provenance and actual equality of `X` to the
  selected named layer.
- Gate: `uv run pre-commit run --hook-stage pre-push --all-files`.

## Open questions

- Which DIA-NN version first supports `--dda`? `v1` (1.8.x) is assumed DIA-only; if a
  1.9.x document exists later, it may need the same override.
- Do other vendors have a usable acquisition marker in their parameter files, or does
  `unknown` stay permanent for them? Revisit when a second vendor needs an override.
