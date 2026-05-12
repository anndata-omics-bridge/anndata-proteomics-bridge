# Plan: Parameter Parsing And Modification Handling

## Goal

Move reusable proteomics parameter-file parsing and modification handling from
ProteoBench into APB (`anndata_proteomics_bridge`).

These topics are interconnected:

1. Parameter files define searched fixed/variable modifications and other search
   settings.
2. Quantification result files contain modified sequence / peptidoform columns
   that need consistent parsing and normalization.
3. The same normalized modification identities should support both SDRF metadata
   export and ProForma sequence rendering.

This is a planning TODO only. Do not start implementation until the migration
plan is approved.

## Target Ownership

- APB owns reusable vendor parameter parsers.
- APB owns reusable proteomics modification parsing and normalization.
- ProteoBench owns submission workflows, UI, benchmarking behavior, and
  ProteoBench-specific presentation.
- ProteoBench is **not** changed in this stage. Wiring ProteoBench to consume
  APB (and eventually removing ProteoBench's duplicated parsers) is a separately
  approved future stage. This plan delivers a standalone, fully-tested
  implementation in APB only.

## Standards Context

### SDRF-Proteomics

SDRF-Proteomics uses repeated `comment[modification parameters]` columns for
searched modifications. The value format is an SDRF key=value object, not
ProForma.

Example:

```text
NT=Oxidation;AC=UNIMOD:35;TA=M;MT=variable;PP=Anywhere
```

Fields:

- `NT`: name term, e.g. `Oxidation`
- `AC`: accession, e.g. `UNIMOD:35`
- `TA`: target amino acid or site, e.g. `M`
- `MT`: modification type, e.g. `fixed` or `variable`
- `PP`: position qualifier, e.g. `Anywhere`

Parse order-insensitively, but export in canonical order:

```text
NT=<name>;AC=<accession>;MT=<fixed|variable>;TA=<target>;PP=<position>
```

### ProForma

ProForma is appropriate for modified peptide/proteoform sequences. It is the
right output format for localized peptidoforms from quantification result
columns, but it should not be the only internal representation because searched
parameter-file modifications are often rules, not concrete sequence
occurrences.

### Vocabulary Target

Use **Unimod** as the preferred accession vocabulary for search-engine
modifications because it is also the practical SDRF-Proteomics target for these
examples. Preserve PSI-MOD accessions when the source data already provides them
or when no Unimod equivalent is known.

Do **not** map directly to UniProt as the primary target. UniProt features may
be useful downstream, but they are not the right interchange target for search
settings or ProForma peptidoforms.

## Current ProteoBench State

### Parameter Parsers

ProteoBench currently owns reusable parser implementations under:

```text
ProteoBench/proteobench/io/params/
```

Relevant quant parser modules include:

- `alphadia.py`
- `alphapept.py`
- `diann.py`
- `fragger.py`
- `i2masschroq.py`
- `maxquant.py`
- `msangel.py`
- `peaks.py`
- `proline.py`
- `quantms.py`
- `sage.py`
- `spectronaut.py`
- `wombat.py`

De novo parsers exist in ProteoBench but are **out of scope** for this plan.

The parser contract is currently centered on `ProteoBenchParameters`, with
expected outputs tested by CSV fixtures in:

```text
ProteoBench/test/test_parse_params_*.py
ProteoBench/test/params/
```

ProteoBench also has a parameter overview matrix:

```text
ProteoBench/docs/parsing_overview.tsv
```

### Modification Rules In ProteoBench TOMLs

ProteoBench TOMLs with `[modifications_parser]` live under:

```text
ProteoBench/proteobench/io/parsing/io_parse_settings/
```

There are 60 TOML files with modification parser configuration. Grouped by tool:

| Tool | Files | Unique dictionaries | Parse column(s) | Pattern(s) |
|---|---:|---:|---|---|
| `alphapept` | 2 | 1 | `Modified sequence` | `([a-z]+)` |
| `diann` | 8 | 1 | `Sequence` | `\(([^()]*)\)` |
| `fragpipe` | 2 | 1 | `Modified Sequence` | `(?<=\[).+?(?=\])` |
| `fragpipe_DIA` | 5 | 2 | `Modified Sequence` | `(?<=\[).+?(?=\])`, `\[([^]]+)\]` |
| `i2massChroQ` | 2 | 1 | `proforma` | `(?<=\[).+?(?=\])` |
| `maxdia` | 5 | 2 | `Modified sequence` | `\([^()]*\)|\([^()]*\([^()]*\)[^()]*\)` |
| `maxquant` | 2 | 1 | `Modified sequence` | `\([^()]*\)|\([^()]*\([^()]*\)[^()]*\)` |
| `metamorpheus` | 2 | 1 | `Modified sequence` | `\[(.*?)\]` |
| `msaid` | 5 | 1 | `Sequence` | `\[(.*?)\]` |
| `peaks` | 9 | 1 | `Sequence` | `(?<=\().+?(?=\))` |
| `sage` | 3 | 1 | `Sequence` | `(?<=\[).+?(?=\])`, `\[([^]]+)\]` |
| `spectronaut` | 6 | 1 | `Sequence` | `\[(.*?)\]` |
| `wombat` | 3 | 1 | `Sequence` | `(?<=\[).+?(?=\])` |

Initial APB migration should cover tools with packaged APB rules:

- DIA-NN: ProteoBench `parse_settings_diann.toml`
- FragPipe: ProteoBench `parse_settings_fragpipe.toml` and relevant
  `parse_settings_fragpipe_DIA.toml`
- MaxQuant / MaxDIA: ProteoBench `parse_settings_maxquant.toml` and
  `parse_settings_maxdia.toml`
- PEAKS: ProteoBench `parse_settings_peaks.toml`
- Spectronaut: ProteoBench `parse_settings_spectronaut.toml`
- WOMBAT: ProteoBench `parse_settings_wombat.toml`

APB-unsupported quant tools should be inventoried but deferred unless explicitly
added to scope. De novo tools are explicitly out of scope.

## Current APB Gaps

APB currently has no parameter parser package, but `src/anndata_proteomics/params/`
already exists as an empty placeholder directory — this is the intended target for the
migration.

APB parsing-rule schema currently has no `modifications_parser` or
`modifications` section. Current APB TOMLs preserve vendor modified sequence
columns as `var` columns, for example:

- `parsing_rules/diann/parse_diann_ion_1.toml`
- `parsing_rules/fragpipe/parse_fragpipe_ion_1.toml`
- `parsing_rules/maxquant/parse_maxquant_ion_1.toml`
- `parsing_rules/peaks/parse_peaks_ion_1.toml`
- `parsing_rules/spectronaut/parse_spectronaut_ion_1.toml`
- `parsing_rules/wombat/parse_wombat_peptidoform_1.toml`

These files need first-class rule sections for modification parsing and
normalization.

### Integration With Existing Rule Schema

The existing pydantic rule schema lives in `src/anndata_proteomics/rules/schema.py`
with loader/validator/registry siblings. The proposed `[modifications]` section must be
added as an **optional** field on the top-level rule model so existing TOMLs continue to
validate unchanged. Conditional validation (parser-mode-specific required fields) should
follow the same pattern already used for `long` vs `wide` `layers` validation.

The existing DIA-NN rule (`parsing_rules/diann/parse_diann_ion_1.toml`) already keeps
both `Stripped_Sequence` and `Modified_Sequence` as var columns. After modification
normalization is added, `axis.var_keys` should be reconsidered: the canonical
identity key likely becomes `proforma_sequence` + `Precursor_Charge`, with the raw
vendor `Modified.Sequence` retained as provenance, not as the identity key.

## Proposed APB Package Shape

Suggested modules:

```text
src/anndata_proteomics/params/
  __init__.py
  model.py
  registry.py
  alphadia.py
  alphapept.py
  diann.py
  fragger.py
  i2masschroq.py
  maxquant.py
  msangel.py
  peaks.py
  proline.py
  quantms.py
  sage.py
  spectronaut.py
  wombat.py

src/anndata_proteomics/modifications/
  __init__.py
  model.py
  sdrf.py
  proforma.py
  vendor_tokens.py
  apply_rules.py
```

De novo parsers are out of scope; do not create `params/denovo/`.

## Parameter Model

Define an APB-owned typed parameter model before moving parser logic. The first
implementation target is:

```text
src/anndata_proteomics/params/model.py
```

This model must be strict Pydantic, not a dynamic dataclass. There should be no
`Any` in public parameter or modification models. If primitive types are not
expressive enough, introduce named domain models rather than weakening the
schema.

Candidate fields are the current ProteoBench parameter fields, including:

- software name and version
- search engine and version
- enzyme
- missed cleavages
- fixed modifications
- variable modifications
- precursor and fragment mass tolerances
- peptide length bounds
- precursor charge bounds
- FDR settings
- MBR / reanalysis setting
- quantification method
- protein inference
- scan window

Hard typing requirements:

- Use `pydantic.BaseModel` with `ConfigDict(extra="forbid")`.
- Do not use `Any`.
- FDR fields are numeric probabilities, not strings:
  - `ident_fdr_psm: Probability | None`
  - `ident_fdr_peptide: Probability | None`
  - `ident_fdr_protein: Probability | None`
- Charges are integers:
  - `min_precursor_charge: PositiveInt | None`
  - `max_precursor_charge: PositiveInt | None`
- Peptide lengths and missed cleavages are non-negative integers.
- m/z fields are non-negative floats:
  - `min_precursor_mz: NonNegativeFloat | None`
  - `max_precursor_mz: NonNegativeFloat | None`
  - `min_fragment_mz: NonNegativeFloat | None`
  - `max_fragment_mz: NonNegativeFloat | None`
- Tolerances should not be unstructured strings when parsed into the core model.
  Use a typed model such as `MassTolerance`.
- Predictor library, quantification method, protein inference, software names,
  and versions are strings or constrained enums where the vocabulary is stable.
- Unknown or unsupported parsed fields should be represented as typed warnings
  or typed `UnparsedParameter` entries, not free-form `Any`.

Suggested domain models:

```python
class Probability(BaseModel):
    value: float  # 0 <= value <= 1

class MassTolerance(BaseModel):
    lower: float | None = None
    upper: float | None = None
    value: float | None = None
    unit: Literal["ppm", "Da", "Th"]
    mode: Literal["absolute", "range", "automatic"]

class ChargeRange(BaseModel):
    minimum: PositiveInt | None = None
    maximum: PositiveInt | None = None

class MzRange(BaseModel):
    minimum: NonNegativeFloat | None = None
    maximum: NonNegativeFloat | None = None
```

Validation requirements:

- FDR values must be `0 <= value <= 1`.
- Charges must be positive integers.
- m/z values must be non-negative.
- Range models must enforce `minimum <= maximum` when both bounds are present.
- Mass tolerances must be non-negative unless the model represents a signed
  lower/upper range where `lower <= upper`.
- Parser adapters may accept vendor strings, but they must convert to typed
  models before returning APB parameter objects.
- Preserve current ProteoBench CSV-compatible serialization during migration via
  explicit serializers, not by weakening the model.
- Avoid adding public API beyond the parser entry points needed by ProteoBench
  and APB tests.

## Modification Models

Use separate models for searched modifications and sequence occurrences.
These are Pydantic models and must not use `Any`.

### SearchedModification

For parameter files and SDRF export:

- `name`
- `accession`
- `mod_type`: `fixed`, `variable`, or `unknown`
- `target`
- `position`
- `mass_delta`
- `source`

This model exports to SDRF key=value strings.
Validation rules:

- `accession` must be a controlled vocabulary accession such as `UNIMOD:35` or
  `MOD:00425` when present.
- `mod_type` is a constrained literal, not a string bag.
- `mass_delta` is a float when present.
- `target` and `position` should be constrained literals or small typed value
  objects once the supported vocabulary is known.

### ModificationOccurrence

For modified peptide/peptidoform sequences:

- `name`
- `accession`
- `target_residue`
- `sequence_index`
- `position`
- `mass_delta`
- `source_token`

This model exports to ProForma when localization is known.
Validation rules:

- `sequence_index` must be non-negative when present.
- `accession` follows the same controlled-vocabulary rule as
  `SearchedModification`.
- Unknown tokens are carried explicitly in `unknown_tokens`; they are not hidden
  in loosely typed fields.

### ModifiedSequence

For quantification results:

- `stripped_sequence`
- `proforma_sequence`
- `occurrences`
- `source_sequence`
- `unknown_tokens`

## APB TOML Schema Addition

Add an optional top-level `[modifications]` section to APB parsing-rule TOMLs.

Suggested schema:

```toml
[modifications]
source_column = "Modified.Sequence"
sequence_column = "Stripped.Sequence"
output_column = "proforma_sequence"
parser = "token_regex"
token_pattern = "\\(([^()]*)\\)"
token_position = "after_residue"
case_sensitive = false
unknown_policy = "preserve"

[[modifications.map]]
token = "(unimod:35)"
name = "Oxidation"
accession = "UNIMOD:35"
target = "M"
position = "Anywhere"

[[modifications.map]]
token = "(unimod:4)"
name = "Carbamidomethyl"
accession = "UNIMOD:4"
target = "C"
position = "Anywhere"
```

Schema fields:

- `source_column`: vendor column to parse.
- `sequence_column`: optional stripped sequence column used for validation.
- `output_column`: normalized sequence column to write into the DataFrame before
  AnnData assembly, default `proforma_sequence`.
- `parser`: parsing mode. Initial values:
  - `token_regex`: extract vendor tokens with regex and map tokens.
  - `already_proforma`: validate/canonicalize existing ProForma.
  - `separate_mod_column`: use base sequence plus separate modification column.
- `token_pattern`: regex used by `token_regex`.
- `token_position`: how tokens attach to residues:
  - `before_residue`
  - `after_residue`
  - `n_term`
  - `c_term`
  - `embedded`
  - `unknown`
- `case_sensitive`: whether token matching is case-sensitive.
- `unknown_policy`: `preserve`, `drop`, or `error`.
- `modifications.map`: token-to-normalized-modification mapping.

Do not migrate ProteoBench entries as plain `token -> name` dictionaries.
Migrate them as structured identity records:

- `token`: exact vendor token or mass token found in the result file
- `name`: canonical display name
- `accession`: controlled vocabulary accession
- `target`: residue or terminus target when known
- `position`: localization constraint such as `Anywhere`, `Protein N-term`,
  `Peptide N-term`, `C-term`, or `unknown`
- `mass_delta`: optional numeric delta for mass-token mappings

Keep old ProteoBench names (`before_aa`, `isalpha`, `isupper`) out of the APB
public schema unless they describe a real domain concept. They are parser
implementation details and should be replaced with clearer fields such as
`token_position` and `case_sensitive`.

### FragPipe Mapping Example

ProteoBench source:

```toml
[modifications_parser]
"parse_column" = "Modified Sequence"
"before_aa" = false
"isalpha" = true
"isupper" = true
"pattern" = "(?<=\\[).+?(?=\\])"
"modification_dict" = {
  "57.0215" = "Carbamidomethyl",
  "57.0216" = "Carbamidomethyl",
  "15.9949" = "Oxidation",
  "-17.026548" = "Gln->pyro-Glu",
  "-18.010565" = "Glu->pyro-Glu",
  "42.0106" = "Acetyl"
}
```

APB target:

```toml
[modifications]
source_column = "Modified Sequence"
parser = "token_regex"
token_pattern = "(?<=\\[).+?(?=\\])"
token_position = "after_residue"
case_sensitive = false
unknown_policy = "preserve"
output_column = "proforma_sequence"

[[modifications.map]]
token = "57.0215"
name = "Carbamidomethyl"
accession = "UNIMOD:4"
target = "C"
position = "Anywhere"
mass_delta = 57.0215

[[modifications.map]]
token = "57.0216"
name = "Carbamidomethyl"
accession = "UNIMOD:4"
target = "C"
position = "Anywhere"
mass_delta = 57.0216

[[modifications.map]]
token = "15.9949"
name = "Oxidation"
accession = "UNIMOD:35"
target = "M"
position = "Anywhere"
mass_delta = 15.9949

[[modifications.map]]
token = "-17.026548"
name = "Gln->pyro-Glu"
accession = "UNIMOD:28"
target = "Q"
position = "N-term"
mass_delta = -17.026548

[[modifications.map]]
token = "-18.010565"
name = "Glu->pyro-Glu"
accession = "UNIMOD:27"
target = "E"
position = "N-term"
mass_delta = -18.010565

[[modifications.map]]
token = "42.0106"
name = "Acetyl"
accession = "UNIMOD:1"
target = "N-term"
position = "N-term"
mass_delta = 42.0106
```

This structured mapping supports multiple exports:

- ProForma for localized peptidoform strings, e.g. `M[UNIMOD:35]PEPTIDE`.
- SDRF searched-modification metadata, e.g.
  `NT=Oxidation;AC=UNIMOD:35;MT=variable;TA=M;PP=Anywhere`.
- Internal AnnData provenance in `uns` without losing the original vendor token.

One caveat: quantification result columns usually identify observed/localized
modification occurrences, but do not reliably say whether a modification was
searched as fixed or variable. `MT=fixed|variable` should therefore come from
parameter-file parsing when available. Result-column normalization can keep
`mod_type = "unknown"` unless the tool rule has authoritative knowledge.

## Implementation Workflow

Work proceeds in **segments**. For each segment, the loop is:

1. Implement the code for the segment.
2. Self-review the diff (naming, scope, no premature abstractions, no unrelated changes).
3. Write tests against ProteoBench fixtures in
   `/Users/wolski/projects/anndata_bridge/ProteoBench/test/params/` where applicable.
4. Run the tests:
   `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/<relevant>`.
5. Re-review the resulting code, fix any issues found, rerun tests.
6. For segments that touch the conversion pipeline or packaged TOMLs, also run
   `uv run python tools/generate_report.py` and confirm it still completes.
7. Commit the segment.
8. Move to the next segment.

### Segment List

- **S1 — Parameter model.** `params/model.py` with `Parameters` (Pydantic).
  Tests for construction + serialization.
- **S2 — Common helpers + first parser (Sage).** `params/_common.py` (tolerance,
  bool, numeric, enzyme helpers). Port `proteobench/io/params/sage.py` →
  `params/sage.py`. Test against `ProteoBench/test/params/sage_parameterfile.json`
  with the `sage_parameterfile.csv` expected output.
- **S3 — Modification models + normalization.** `modifications/model.py`
  (`SearchedModification`, `ModificationOccurrence`, `ModifiedSequence`),
  `modifications/apply_rules.py` (token→model), `modifications/sdrf.py`,
  `modifications/proforma.py`. Unit tests for round-trip token → model → ProForma
  and model → SDRF key=value.
- **S4 — TOML schema extension.** Add optional `[modifications]` block to
  `rules/schema.py` with conditional validation. Schema tests.
- **S5 — Pipeline integration.** Converter accepts optional `params_path`.
  Apply modification normalization between `read_table` and assemble. Populate
  `uns['<software>']['search_parameters']`. Update DIA-NN packaged rule to use
  `proforma_sequence` in `axis.var_keys`. End-to-end test + `generate_report.py`.
- **S6 — Port remaining parameter parsers.** One commit per parser, each
  validated against its ProteoBench fixture pair. Order: structured formats
  first (AlphaPept YAML, AlphaDIA, quantms, MetaMorpheus, WOMBAT, MaxQuant XML),
  then text/log/report formats (DIA-NN log, FragPipe `.workflow`, PEAKS,
  Spectronaut, MSAID, MSAngel, ProlineStudio, i2MassChroQ).
- **S7 — Migrate packaged rule TOMLs.** Add `[modifications]` to DIA-NN,
  FragPipe, MaxQuant, PEAKS, Spectronaut, WOMBAT rules. Run
  `tools/generate_report.py` end-to-end on every packaged rule.

If a segment turns out larger than expected, split it; do not merge segments.
Each segment ends in a commit that leaves the tree in a green state.

## Migration Steps

1. **Inventory ProteoBench sources**
   - Inventory every parameter parser, fixture, expected CSV, and documented
     parameter in `docs/parsing_overview.tsv`.
   - Extract all `[modifications_parser]` TOML sections.
   - Deduplicate modification mapping dictionaries by tool.

2. **Define APB models**
   - Create `src/anndata_proteomics/params/model.py`.
   - Define the APB parameter model there first, before moving any parser logic.
   - Define `SearchedModification`, `ModificationOccurrence`, and
     `ModifiedSequence`.
   - Define SDRF and ProForma export helpers as separate adapters.

3. **Extend APB TOML schema**
   - Add `ModificationRule`, `ModificationMapEntry`, and `Modifications`
     Pydantic models.
   - Add validation:
     - `token_regex` requires `source_column`, `token_pattern`, and mappings.
     - `already_proforma` requires `source_column` but no token map.
     - `separate_mod_column` requires source columns defined by the rule.
     - `unknown_policy` must be explicit.

4. **Port one parser first**
   - Port one low-risk parameter parser first, preferably `sage.py`, with its
     fixtures and tests.
   - Keep expected output byte-for-byte equivalent where practical.
   - Ensure fixed/variable modifications become `SearchedModification` objects.

5. **Add modification normalization code**
   - Implement vendor token extraction and mapping.
   - Add ProForma rendering where localization is known.
   - Add SDRF rendering for searched modifications.
   - Reuse established libraries if available and practical; do not implement a
     full ProForma parser from scratch unless no dependency fits the needed
     scope.

6. **Apply normalization before conversion**
   - In the conversion pipeline, after `read_table` and before
     `recognize/convert`, apply modification normalization for the matched rule.
   - Add normalized columns, do not silently replace source vendor columns.
   - Update `axis.var_keys` for APB rules to use normalized ProForma columns
     where possible.

7. **Migrate APB packaged TOMLs**
   - DIA-NN: map `(unimod:35)`, `(unimod:1)`, `(unimod:4)`.
   - MaxQuant / MaxDIA: map `(ox)`, `(ac)`, `(oxidation (m))`,
     `(acetyl (protein n-term))`, and related MaxDIA variants.
   - FragPipe / FragPipe DIA: map numeric mass tokens such as `57.0215`,
     `15.9949`, `-17.026548`, `-18.010565`, `42.0106`.
   - PEAKS: map `+57.02`, `+15.99`, `-17.026548`, `-18.010565`, `+42.01`.
   - Spectronaut: map `[oxidation (m)]`, `[acetyl (protein n-term)]`,
     `[carbamidomethyl (c)]`.
   - WOMBAT: map numeric tokens such as `160`, `147`, `80`, `111`, `43`.

8. **Port remaining quant parameter parsers**
   - JSON/YAML-like formats first: Sage, AlphaDIA, AlphaPept, Wombat, quantms.
   - XML formats next: MaxQuant / MaxDIA.
   - Text/log/table formats next: DIA-NN, FragPipe, PEAKS, Spectronaut,
     MSAngel, ProlineStudio, i2MassChroQ.

9. **Add registry and CLI only after library API stabilizes**
   - Add a parser registry mapping software names and file signatures to parser
     functions.
   - Add CLI support only after the API is stable, e.g.:

     ```bash
     anndata-proteomics parse-params <parameter-file> --software Sage
     ```

10. **Equivalence verification (APB-side only)**
    - Run APB parsers against the same fixtures used by ProteoBench parser
      tests and compare serialized output for equivalence.
    - ProteoBench itself is **not** modified at this stage: no shim, no
      dependency pin, no deletion of ProteoBench parsers.
    - Wiring ProteoBench to consume APB is deferred to a separately approved
      future stage.

## Test Strategy

- Port existing ProteoBench expected CSV fixtures into APB tests or consume them
  from ProteoBench during the initial migration.
- Add unit tests for the normalized modification model:
  - vendor token to internal model
  - internal model to SDRF key=value encoding
  - internal model to ProForma where unambiguous
  - unsupported/ambiguous cases preserve the original source string
- Add unit tests for vendor token parsing for each APB-supported tool.
- Add one registry test per parser to verify software detection and parser
  dispatch.
- Add converter tests proving `var` identity uses normalized modified sequence
  columns where configured.
- Add compatibility tests in ProteoBench proving existing `extract_params`
  calls still return the same serialized values.

## Out Of Scope For First Pass

- Rewriting ProteoBench UI.
- Changing submitted JSON schemas.
- Changing benchmark scoring.
- Rich parameter validation beyond current parsed fields.
- Converting every modification to ProForma if localization or target
  information is ambiguous.
- De novo parser migration (explicitly out of scope).

## Parameter ↔ Quant Result Integration

Searched modifications from parameter files (`SearchedModification`) are **per-experiment**,
not per-row. They do not fit `obs` or `var`; they should land in
`uns['<software>']['search_parameters']` alongside other parsed search settings. The
`column_roles` ADR in `anndata_omics_bridge/docs/adr_tool_specific_views.md` is the
relevant precedent for `uns` shape.

The conversion pipeline currently does not consume parameter files at all. Two options
for how it gets there:

- **Coupled.** A converter accepts an optional `params_path` argument; if provided,
  the parser is run and result is stored in `uns`. Same call site as `read_table`.
- **Decoupled.** Parameter parsing is a standalone API (`parse_params(path, software=...)`)
  with results merged into an existing AnnData via a separate helper. Useful when the
  parameter file and quant file ship separately, or when users only want SDRF metadata.

The plan should pick one as the default and note the other as future work.

## Dependency Direction (This Stage: APB-Only)

At this stage, the work is **additive on the APB side only**. ProteoBench is not
modified, not refactored, and nothing is deleted from it. ProteoBench continues to
run exactly as it does today against its own in-tree parsers.

Concrete rules for this stage:

- APB stays free of any `proteobench` imports.
- APB's `params/` package is implemented as a full, standalone parser layer with
  its own tests and fixtures (fixtures may be copied from ProteoBench, but APB
  does not import ProteoBench code).
- ProteoBench's `pyproject.toml` is **not** modified to depend on APB yet.
- ProteoBench's `proteobench/io/params/*.py` is **not** modified, shimmed, or
  deleted.
- Equivalence is proven by running APB parsers against the same fixtures that
  ProteoBench tests use and comparing serialized output — without touching
  ProteoBench's runtime.

A later, separately approved stage will handle wiring ProteoBench to consume APB
(shim layer, dependency pin, eventual deletion of duplicated code). That stage is
out of scope here.

## Resolved Decisions (2026-05-11)

- **Model type.** Parameter model uses Pydantic, matching the existing
  `rules/schema.py` validation layer.
- **Pipeline integration.** Coupled: the converter takes an optional `params_path`
  argument. When provided, `parse_params` runs and the result lands in
  `uns['<software>']['search_parameters']`. No separate merge API in the first pass.
- **Canonical var identity.** ProForma is the canonical `var` modified-sequence
  column for every tool where localization is available. `axis.var_keys` for APB
  rules updates to `proforma_sequence` + charge (or equivalent). The original
  vendor `Modified.Sequence` column is retained as provenance in `columns.var`,
  not as an identity key.
- **Unknown tokens.** Default `unknown_policy = "preserve"`: log a warning and
  keep the original vendor token verbatim inside the ProForma string. `error` and
  `drop` remain available per rule for stricter contexts.
- **Vocabulary.** Canonicalize to Unimod accessions when known; preserve PSI-MOD
  only when no Unimod equivalent exists. (Aligned with the Standards Context
  section.)
- **SDRF export.** Layered exporter on top of parsed parameters
  (`anndata_proteomics.modifications.sdrf`), not part of the parser API itself.

## Token Lookup Rule

Each `modifications.map` entry is identified by the tuple
`(mass_delta, target, position)`, not by `mass_delta` alone. At parse time the
matcher reads the token's mass, the adjacent residue from the sequence, and the
position (`Anywhere`, `N-term`, `C-term`) and matches against this tuple. This
makes mass-only tokens unambiguous in practice:

- `42.0106` on K matches Acetyl-K; `42.0106` at N-term matches Acetyl-Nterm.
- `79.9663` on S/T matches Phospho; `79.9663` on Y matches Phospho-Y (or Sulfo-Y
  if the map entry targets Y with `Sulfo`).
- `-17.026548` on Q at N-term matches Gln→pyro-Glu.

If no map entry matches the tuple, `unknown_policy` applies.

## Design Decisions Still Open

- Should APB support auto-detection of parameter-file software, or require
  explicit `software=` selection first? (Default assumption: explicit first;
  auto-detection registry added later if needed.)

## Verification

Targeted APB checks:

```bash
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider \
  tests/test_modifications_*.py \
  tests/test_params_*.py \
  tests/test_rule_validate.py \
  tests/test_converters_e2e.py
```

Compatibility checks after ProteoBench delegation:

```bash
pytest ProteoBench/test/test_parse_params_*.py
pytest ProteoBench/test/test_*parsing*.py
```

## References

- SDRF-Proteomics specification:
  https://sdrf.quantms.org/specification.html
- SDRF terms reference:
  https://sdrf.quantms.org/sdrf-terms.html
- HUPO-PSI ProForma:
  https://www.psidev.info/proforma
- ProForma 2.0 paper/preprint:
  https://arxiv.org/abs/2109.11352
