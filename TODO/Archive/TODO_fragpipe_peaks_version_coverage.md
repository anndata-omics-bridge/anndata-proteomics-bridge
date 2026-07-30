# FragPipe / PEAKS version coverage, PEAKS wide-schema correctness, and missing-version contract

> Make every cached FragPipe and PEAKS ion submission convert correctly, not merely pass rule
> recognition: broaden the known version families, read delimiter-mislabeled PEAKS exports, prevent
> PEAKS summary columns from becoming observations, and record genuine missing-version selection
> explicitly.

**Date:** 2026-07-02 · **Corrected after live review:** 2026-07-24 · **Source:** `apb_studio`
full-corpus `make run` plus direct APB reproductions over the current `test_data_download` cache.
**Owner:** apb (readers, wide conversion, parsing rules, parameter provenance).
**Status:** implemented and verified 2026-07-24.

---

## Verified ground truth

The corpus contains five cached FragPipe and seven cached PEAKS submissions.

### Version coverage

| Vendor / case | Cached submissions | Parsed version | Current result |
|---|---:|---|---|
| FragPipe | 2 | `22.1-build02` | covered |
| FragPipe | 1 | `22.0` | `no rule covers software version` |
| FragPipe | 2 | `23.0` | `no rule covers software version` |
| PEAKS | 4 | `13 20250520` | `no rule covers software version` |
| PEAKS | 1 | `13 20250515` | `no rule covers software version` |
| PEAKS | 2 | `None` | one reaches the existing unversioned rule; one fails column matching |

The versioned FragPipe inputs fit the existing FragPipe ion rule. The versioned PEAKS inputs pass
the current subset recognizer when read with the correct delimiter. The regexes are therefore too
narrow for the known versions:

- FragPipe: `^22\.1-build02$`
- PEAKS: `^13$`

### The PEAKS `None` cases are not the same failure

There are two cached PEAKS parameter files whose parser succeeds with
`Parameters.software_version is None`:

```
# DIA-PASEF — comma-delimited .csv
test_data_download/json_dir/Results_quant_ion_DIA_diaPASEF/5691d485356c1abcd8efd2320a1666752a19f50b/

# DIA/AIF — comma-delimited content mislabeled as .txt
test_data_download/json_dir/Results_quant_ion_DIA_AIF/b5fddd9b5d27918e8d31ec07bcf599cbd214027a/
```

`resolve_rule_locator("peaks", "ion", None)` already returns the sole PEAKS rule document. The
DIA/AIF failure occurs earlier in practice because `read_table()` treats every `.txt` as TSV. It
reads that comma-delimited file as `112140 × 1`; reading it as CSV produces 38 columns and the rule
passes recognition. This is a reader bug, not a missing-version resolver bug.

### Passing subset recognition is not sufficient for PEAKS

The previous plan incorrectly treated extra wide columns as harmless. `matches()` only decides
whether required inputs are present. During conversion, `convert_wide()` currently builds the
observation axis from the union of sample tokens matched by every layer. The PEAKS patterns are
broad enough to capture aggregate and auxiliary headers as samples:

- `Group N Normalized Area`
- `Condition A/B Normalized Area`
- `Best AScore`
- optional-layer names with `.raw` / `_raw` suffixes that do not match the x-layer token
- one malformed DIA-PASEF auxiliary header that concatenates two run names

Concrete direct conversions with the current rule:

| Input | Real run-level columns | Current observations | Wrong additions |
|---|---:|---:|---|
| DDA `9d136133…` | 6 | 9 | `Group 1`, `Group 2`, `Best` |
| DIA/AIF `b5fddd9…` (read as CSV) | 6 | 12 | `Group 1`–`Group 6` |
| DIA/Astral `aa1d53e…` (`13 20250515`) | 6 | 14 | two groups plus six `_raw` auxiliary-layer duplicates |
| DIA-PASEF `5691d485…` (version absent) | 6 | 13 | six groups plus one malformed auxiliary-layer token |

Broadening the PEAKS version regex without fixing this would turn clean version failures into
silently wrong AnnData objects. The PEAKS sample-axis fixes must land in the same change.

### Parameter-state ambiguity

`param_version()` currently returns `None` for two distinct states:

1. parameter parsing succeeded and the file genuinely omitted the version;
2. parameter parsing raised because the file was malformed or belonged to another tool.

Only state 1 is the new missing-version contract. State 2 must retain the existing best-effort
parse-error policy and `search_parameters_error`; it must not gain a new fallback among
non-equivalent version variants accidentally.

---

## Requirements

### In scope

- FragPipe ion versions `22.x` and `23.x`.
- PEAKS ion major version 13, including `13 <build-date>` and `13.x` forms.
- Correct comma/tab detection for generic `.txt` tabular input.
- Correct PEAKS run-axis extraction across all seven cached submissions.
- Column-based rule selection only when parsing succeeded and the software version is genuinely
  absent.
- Distinct parameter states for version present, version missing, and parameter parse error.
- Rule-selection and missing-version provenance in `uns`.
- Search-parameter storage that always writes every model field and represents missing values as
  JSON `null`.

### Out of scope

- DIA-NN parser-error degradation beyond preserving its current behavior and
  `search_parameters_error`.
- WOMBAT level scaffolding, already fixed in `apb_studio`.
- General vendor delimiter detection beyond comma/tab text files.
- New public conversion APIs. Any selection/parameter-resolution carrier should remain internal.

---

## Design

### 1. Text delimiter handling

Keep extension-directed reading for `.csv`, `.tsv`, and `.parquet`. For `.txt`, detect comma versus
tab from the file content and then pass the detected delimiter to pandas.

- Use one canonical internal detector, based on the standard-library CSV dialect machinery or an
  equivalent existing facility.
- Reuse it for both full table reads and header-only cache inspection.
- Preserve UTF-8 BOM handling.
- Do not add a PEAKS-only reader or rename/copy the input file as a workaround.

This fixes the root cause in `readers/`, and removes the duplicated assumption in
`tests/conftest.py::_read_headers` that every non-Parquet input is tab-delimited.

### 2. Wide observation-axis contract

For a wide rule, `axis.x_layer` defines the observation axis.

- Derive `sample_order` only from x-layer matches, not from the union of all layers.
- Optional-layer tokens that are absent from the x-layer axis must never create observations.
- Ignore such auxiliary-only tokens with a WARNING that identifies the layer and token.
- Reindex every other layer onto the x-layer observation axis; missing optional measurements remain
  missing.

This is the general upstream fix for auxiliary layers expanding the observation axis.

### 3. PEAKS layer patterns

Tighten the PEAKS layer regexes so only run-level columns match.

- Require the run-level `LFQ_` prefix.
- Exclude `Group`, `Condition`, and `Best` summary columns structurally through the positive
  run-prefix match, not a growing list of negative special cases.
- Capture optional `.raw` / `_raw` suffixes outside the `sample` group so a run has the same sample
  token in the x, m/z, RT, and AScore layers.
- Keep AScore optional because the DIA exports do not provide it.

The expected PEAKS observations are the six run-level samples for each cached module. Exact names
must agree with the independent module annotations/run list.

### 4. Version regexes

- FragPipe ion: `^22\.1-build02$` → `^2[23]\.`.
- PEAKS ion: `^13$` → `^13(?:\.|\s|$)`.

This deliberately covers known major families while keeping new majors strict. FragPipe `24.x` and
PEAKS `14.x` remain uncovered even when their columns happen to fit.

### 5. Parameter resolution and rule selection

Parse the parameter file once in the orchestration path and preserve three internal states:

- **version present:** select only a rule whose version regex matches; an uncovered version is a
  hard failure with no column fallback.
- **version genuinely missing:** evaluate the `(slug, level)` variants with the complete
  header-compatibility predicate, including fragment label-column requirements; select exactly one
  matching variant and fail on zero or multiple matches.
- **parameter parse error:** preserve the pre-existing unversioned/equivalent-rule behavior and
  `search_parameters_error`, but do not use the new missing-version-only fallback to choose among
  non-equivalent variants.

The complete compatibility predicate should live in one converter-layer helper and be reused by
normal selection, missing-version selection, and `matching_rules()`. Do not put DataFrame/header
semantics into the rule-document loader.

### 6. Provenance

Every orchestration-selected output records:

```
uns["anndata_proteomics"]["rule_selection_method"]
```

with one of:

- `"software_version"`
- `"columns"`
- `"version_unavailable"` for preserved parse-error behavior
- `"rule_config"` for explicit `--rule-config`

Parameter version state is recorded separately:

```
uns["anndata_proteomics"]["search_parameters_version_status"]
```

with `"present"`, `"missing"`, or `"parse_error"`. A missing version logs a WARNING once per
conversion. A parse error continues to carry `search_parameters_error`.

For a single-level conversion, write these keys on the AnnData object. For MuData conversion, write
them on the MuData root and each generated modality so either storage unit retains its provenance.

### 7. Search-parameter storage

Remove `exclude_none=True` from `write_search_parameters()` so the JSON object always contains
every `Parameters` field.

- A genuinely absent `software_version` remains `None` in the model and JSON `null` on disk.
- Do not inject an `"ERROR: ..."` string into `software_version`: that would turn an error/status
  into version data and surface it as a real version in descriptive summaries.
- The structured `search_parameters_version_status` key and WARNING provide the missingness signal.
- Selection must use the parsed model value before storage, never a reporting sentinel.

`software_version` therefore remains optional-but-observed in
`TODO/TODO_params_model_review.md`; it must not become required while real supported parameter
files omit it.

---

## Acceptance criteria

- All five cached FragPipe submissions convert at ion level.
- All seven cached PEAKS submissions convert at ion level.
- Every cached PEAKS output has exactly the independently expected run observations:
  - no `Group *`, `Condition *`, or `Best` observations;
  - no duplicated `.raw` / `_raw` observations;
  - optional layers align to the x-layer axis.
- The DIA/AIF comma-delimited `.txt` input is read as 38 columns and converts without renaming it.
- FragPipe `22.0`, `22.1-build02`, and `23.0` resolve; fabricated `24.0` fails.
- PEAKS `13`, `13 20250515`, `13 20250520`, and a representative `13.x` resolve; fabricated `14.0`
  fails.
- Both real PEAKS version-missing inputs:
  - select the rule by columns;
  - log a WARNING;
  - store `software_version` as JSON `null`;
  - store every `Parameters` field;
  - record `rule_selection_method = "columns"` and
    `search_parameters_version_status = "missing"`.
- A malformed/wrong-tool parameter file is classified as `parse_error`, not `missing`, and does not
  gain the new non-equivalent-variant column fallback.
- Existing `search_parameters_error` behavior remains intact.
- `pytest tests/` and the full pre-commit gate (Ruff, Pyright, coverage) pass.

---

## Implementation plan

- [x] **Add failing regression tests first.**
  - comma-delimited `.txt` and tab-delimited `.txt`;
  - wide x-layer-only observation axis;
  - auxiliary-only layer tokens do not create observations and emit a warning;
  - PEAKS summary columns and raw suffix normalization;
  - exact observation names for the four concrete failures above.
- [x] **Fix text reading upstream.**
  - add/reuse one delimiter detector in `readers/tabular.py`;
  - route `.txt` through it in `readers/dispatch.py`;
  - reuse the same logic for header-only test-data inspection.
- [x] **Fix the general wide-axis contract.**
  - derive observations from `axis.x_layer`;
  - reindex auxiliary layers;
  - warn and ignore auxiliary-only sample tokens.
- [x] **Tighten PEAKS layer regexes.**
  - require `LFQ_`;
  - normalize `.raw` / `_raw` in the capture;
  - verify all seven cached PEAKS schemas.
- [x] **Broaden the known version regexes.**
  - FragPipe `^2[23]\.`;
  - PEAKS `^13(?:\.|\s|$)`.
- [x] **Separate parameter resolution states.**
  - parse once;
  - distinguish present, genuinely missing, and parse error;
  - add missing-version-only unique column selection without broadening parse-error fallback.
- [x] **Centralize complete header compatibility.**
  - include the fragment label-column requirement;
  - reuse the helper in all selection paths.
- [x] **Write provenance.**
  - set `rule_selection_method` and `search_parameters_version_status`;
  - cover AnnData, MuData root, and MuData modalities.
- [x] **Store all parameter keys.**
  - serialize `None` as JSON `null`;
  - add a direct wire-format assertion against `Parameters.model_fields`;
  - verify read/write and h5ad/h5mu round-trips.
- [x] **Cross-reference the model review.**
  - update `TODO/TODO_params_model_review.md` so `software_version` remains optional-but-observed.
- [x] **Run gates.**
  - focused reader, wide-converter, rule-resolution, parameter-storage, and cached-corpus tests;
  - full `pytest tests/`;
  - full pre-commit suite.

**Expected implementation files:**

- `src/anndata_proteomics/readers/tabular.py`
- `src/anndata_proteomics/readers/dispatch.py`
- `src/anndata_proteomics/converters/wide.py`
- `src/anndata_proteomics/converters/pipeline.py`
- `src/anndata_proteomics/converters/assemble.py`
- `src/anndata_proteomics/rules/loader.py` only if an internal candidate-enumeration helper is needed
- `src/anndata_proteomics/params/anndata_io.py`
- `src/anndata_proteomics/parsing_rules/fragpipe/rules.json`
- `src/anndata_proteomics/parsing_rules/peaks/rules.json`
- `tests/`
- `TODO/TODO_params_model_review.md`

No parser should synthesize a PEAKS version, no input should be renamed to force a reader, and no
test should accept “non-empty output” as proof of correct PEAKS conversion.

## Verification

- All five cached FragPipe and all seven cached PEAKS submissions convert to exactly six run
  observations.
- Both PEAKS version-missing submissions select by columns and store `software_version: null`.
- Full test/coverage gate: 602 passed, 4 skipped, 100% branch coverage.
- Changed-line coverage: 100%.
- Ruff, strict Pyright, dependency declarations, package smoke, and strict MkDocs all pass.
