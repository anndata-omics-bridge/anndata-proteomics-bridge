# TODO: Replace Nullable Layer Sources With Required `source`

> **DONE (2026-06-23).** Implemented: `Layer.source` (single required field; `categories`
> now defaults to an empty dict), `Modifications` is a `parser`-discriminated union
> (`TokenRegexModifications` / `AlreadyProformaModifications` / `SeparateModColumnModifications`),
> and `Fragments` is a `label_strategy`-discriminated union (`PositionalFragments` /
> `ColumnLabeledFragments`). Packaged TOMLs, converters, recognizer, JSON schema, docs, and
> tests were migrated. Optional top-level sections (`modifications`/`fragments`/
> `sample_name_cleanup`) remain `| None` by design. Kept for design history.

## Goal

Simplify the parsing-rule TOML layer schema by replacing the nullable pair:

```toml
source_column = "..."
column_pattern = "..."
```

with one required string:

```toml
source = "..."
```

The interpretation of `source` is determined by the existing top-level
`input_shape`:

- `input_shape = "long"`: `source` is an exact vendor column name.
- `input_shape = "wide"`: `source` is a regex pattern matching sample columns.

This removes the `str | None = None` shape from `Layer`. More generally, parsing-rule
schema fields should not be modeled as nullable strings when the real contract is
"different TOML shapes for different modes".

## Current Code

This is the current code in `src/anndata_proteomics/rules/schema.py` that this TODO is about.

```python
class Layer(_Strict):
    name: str
    encoding_mode: EncodingMode = "numeric"
    categories: dict[str, int] | None = None
    source_column: str | None = None
    column_pattern: str | None = None

    @model_validator(mode="after")
    def _factor_requires_categories(self) -> Layer:
        if self.encoding_mode == "factor" and not self.categories:
            raise ValueError(
                f"Layer {self.name!r}: encoding_mode='factor' requires non-empty 'categories'."
            )
        return self
```

Current TOML mapping for `Layer`:

```toml
[[layers]]
name = "FG_Quantity"
source_column = "FG.Quantity"
```

```toml
[[layers]]
name = "Intensity"
column_pattern = "^(?P<sample>.+) Intensity$"
```

```toml
[[layers]]
name = "Match_Type"
column_pattern = "^(?P<sample>.+) Match Type$"
encoding_mode = "factor"
categories = { "unmatched" = 0, "MS/MS" = 1, "MBR" = 2 }
```

Current modification-section code:

```python
class Modifications(_Strict):
    source_column: str
    parser: ModificationParser = "token_regex"
    token_pattern: str | None = None
    token_position: TokenPosition = "after_residue"
    case_sensitive: bool = False
    unknown_policy: UnknownPolicy = "preserve"
    sequence_column: str | None = None
    output_column: str = "proforma_sequence"
    map: list[ModificationMapEntry] = Field(default_factory=list)
```

Current TOML mapping for `Modifications`:

```toml
[modifications]
source_column = "Modified.Sequence"
parser = "token_regex"
token_pattern = "\\((?P<token>[^)]+)\\)"
```

```toml
[modifications]
source_column = "EG.ModifiedSequence"
parser = "already_proforma"
output_column = "proforma_sequence"
```

```toml
[modifications]
source_column = "Modified peptide"
sequence_column = "Peptide"
parser = "separate_mod_column"
```

Current fragment-section code:

```python
class Fragments(_Strict):
    value_columns: list[str] = Field(min_length=1)
    label_column: str | None = None
    delimiter: str = ";"
    label_output: str = "fragment_label"
```

Current TOML mapping for `Fragments`:

```toml
[fragments]
value_columns = ["Fragment.Quant.Raw", "Fragment.Quant.Corrected"]
label_column = "Fragment.Info"
delimiter = ";"
label_output = "fragment_label"
```

Current top-level parse-rule code:

```python
class ParseRule(_Strict):
    schema_version: str
    file_version: str
    software_name: str
    software_version: str
    input_shape: InputShape
    quantification_level: QuantificationLevel
    axis: Axis
    columns: Columns
    layers: list[Layer] = Field(min_length=1)
    sample_name_cleanup: SampleNameCleanup | None = None
    modifications: Modifications | None = None
    fragments: Fragments | None = None
```

Current TOML mapping for optional top-level sections:

```toml
[sample_name_cleanup]
pattern = "_Intensity$"
```

```toml
[modifications]
source_column = "EG.ModifiedSequence"
parser = "already_proforma"
```

```toml
[fragments]
value_columns = ["Fragment.Quant.Raw"]
```

## Rationale

The current `Layer` model has:

```python
source_column: str | None = None
column_pattern: str | None = None
```

and `ParseRule` validators enforce the XOR:

- long rules require `source_column` and forbid `column_pattern`
- wide rules require `column_pattern` and forbid `source_column`

That is correct behavior, but the model shape is looser than the TOML contract.
A layer always has a source; what changes is whether the source string is an
exact column name or a regex over wide matrix columns.

Using a neutral `source` name keeps the TOML compact while avoiding the misleading
overload of `source_column` for regex patterns. Internally, the pydantic model should
still make the mode distinction explicit instead of storing mode-specific fields as
`str | None`.

## Current Nullable TOML Model Inventory

This TODO is scoped to `src/anndata_proteomics/rules/schema.py`, i.e. the parsing-rule
TOML model. It does not cover the separate search-parameter pydantic model.

### Exact `type | None` Occurrences

#### `Layer.categories`

Code:

```python
class Layer(_Strict):
    name: str
    encoding_mode: EncodingMode = "numeric"
    categories: dict[str, int] | None = None
    source_column: str | None = None
    column_pattern: str | None = None
```

TOML section:

```toml
[[layers]]
name = "Match_Type"
encoding_mode = "factor"
categories = { "unmatched" = 0, "MS/MS" = 1, "MBR" = 2 }
```

Meaning: `categories` is only meaningful when a layer has `encoding_mode = "factor"`.
For numeric layers it is absent.

Desired model: avoid `dict | None`; use an empty dict default and keep the existing validator
that factor layers must have non-empty categories.

#### `Layer.source_column`

Code:

```python
class Layer(_Strict):
    name: str
    encoding_mode: EncodingMode = "numeric"
    categories: dict[str, int] | None = None
    source_column: str | None = None
    column_pattern: str | None = None
```

TOML section:

```toml
[[layers]]
name = "FG_Quantity"
source_column = "FG.Quantity"
```

Meaning: exact vendor value column for long-format rules.

Desired model: remove. Replace with required `source: str`.

#### `Layer.column_pattern`

Code:

```python
class Layer(_Strict):
    name: str
    encoding_mode: EncodingMode = "numeric"
    categories: dict[str, int] | None = None
    source_column: str | None = None
    column_pattern: str | None = None
```

TOML section:

```toml
[[layers]]
name = "Intensity"
column_pattern = "^(?P<sample>.+) Intensity$"
```

Meaning: regex matching wide-format matrix columns and extracting `sample`.

Desired model: remove. Replace with required `source: str`, interpreted as a regex when
`input_shape = "wide"`.

#### `Modifications.token_pattern`

Code:

```python
class Modifications(_Strict):
    source_column: str
    parser: ModificationParser = "token_regex"
    token_pattern: str | None = None
    token_position: TokenPosition = "after_residue"
    case_sensitive: bool = False
    unknown_policy: UnknownPolicy = "preserve"
    sequence_column: str | None = None
    output_column: str = "proforma_sequence"
    map: list[ModificationMapEntry] = Field(default_factory=list)
```

TOML section:

```toml
[modifications]
source_column = "Modified.Sequence"
parser = "token_regex"
token_pattern = "\\((?P<token>[^)]+)\\)"
```

Meaning: required for `parser = "token_regex"`, forbidden for
`parser = "already_proforma"`.

Desired model: remove nullable field by splitting `Modifications` into parser-specific
entities. `TokenRegexModifications` has required `token_pattern: str`; other parser models do
not have this field.

#### `Modifications.sequence_column`

Code:

```python
class Modifications(_Strict):
    source_column: str
    parser: ModificationParser = "token_regex"
    token_pattern: str | None = None
    token_position: TokenPosition = "after_residue"
    case_sensitive: bool = False
    unknown_policy: UnknownPolicy = "preserve"
    sequence_column: str | None = None
    output_column: str = "proforma_sequence"
    map: list[ModificationMapEntry] = Field(default_factory=list)
```

TOML section:

```toml
[modifications]
source_column = "Modified peptide"
sequence_column = "Peptide"
parser = "separate_mod_column"
```

Meaning: only meaningful for parser shapes where the unmodified sequence is separate from the
modification-token column.

Desired model: remove nullable field by making `sequence_column: str` required on the
parser-specific model that needs it; other parser models should not expose the field.

#### `Fragments.label_column`

Code:

```python
class Fragments(_Strict):
    value_columns: list[str] = Field(min_length=1)
    label_column: str | None = None
    delimiter: str = ";"
    label_output: str = "fragment_label"
```

TOML section:

```toml
[fragments]
value_columns = ["Fragment.Quant.Raw", "Fragment.Quant.Corrected"]
label_column = "Fragment.Info"
delimiter = ";"
label_output = "fragment_label"
```

Meaning: if present, fragment labels come from a packed vendor column. If absent, old DIA-NN
fragment labels are positional (`frag_0`, `frag_1`, ...).

Desired model: remove nullable field by splitting fragment label strategies. Use an explicit
`label_strategy = "column"` variant with required `label_column: str`, and a
`label_strategy = "positional"` variant with no `label_column`.

#### `ParseRule.sample_name_cleanup`

Code:

```python
class ParseRule(_Strict):
    ...
    layers: list[Layer] = Field(min_length=1)
    sample_name_cleanup: SampleNameCleanup | None = None
    modifications: Modifications | None = None
    fragments: Fragments | None = None
```

TOML section:

```toml
[sample_name_cleanup]
pattern = "_Intensity$"
```

Meaning: optional top-level section, only meaningful for wide rules.

Desired model: optional section. This is object-level section absence, not a nullable scalar.
May remain `SampleNameCleanup | None` unless the whole parse-rule model is split into shape
variants with an explicit wide-only field.

#### `ParseRule.modifications`

Code:

```python
class ParseRule(_Strict):
    ...
    layers: list[Layer] = Field(min_length=1)
    sample_name_cleanup: SampleNameCleanup | None = None
    modifications: Modifications | None = None
    fragments: Fragments | None = None
```

TOML section:

```toml
[modifications]
source_column = "EG.ModifiedSequence"
parser = "already_proforma"
output_column = "proforma_sequence"
```

Meaning: optional top-level section. It is required only when computed columns need
modification parsing, such as `how = "proforma_sequence"` or `how = "stripped_sequence"`.

Desired model: optional section can remain nullable. Internally, the `Modifications` section
itself should be parser-specific and avoid nullable scalar fields.

#### `ParseRule.fragments`

Code:

```python
class ParseRule(_Strict):
    ...
    layers: list[Layer] = Field(min_length=1)
    sample_name_cleanup: SampleNameCleanup | None = None
    modifications: Modifications | None = None
    fragments: Fragments | None = None
```

TOML section:

```toml
[fragments]
value_columns = ["Fragment.Quant.Raw", "Fragment.Quant.Corrected"]
label_column = "Fragment.Info"
delimiter = ";"
label_output = "fragment_label"
```

Meaning: optional top-level section for packed fragment expansion. Only valid when
`quantification_level = "fragment"`.

Desired model: optional section can remain nullable. Internally, the `Fragments` section
should be strategy-specific and avoid `label_column: str | None`.

### `str | None = None` Patterns

| Model | Field | Current meaning | Proposed action |
| --- | --- | --- | --- |
| `Layer` | `source_column` | Exact vendor column for long-format layers. Required when `input_shape = "long"`, forbidden when `input_shape = "wide"`. | Remove and replace with required `source: str`. |
| `Layer` | `column_pattern` | Regex for wide-format sample columns. Required when `input_shape = "wide"`, forbidden when `input_shape = "long"`. | Remove and replace with required `source: str`. |
| `Modifications` | `token_pattern` | Regex used only by `parser = "token_regex"`. Must be absent for `parser = "already_proforma"`. | Remove nullable field by splitting `Modifications` into parser-specific models. |
| `Modifications` | `sequence_column` | Helper/source sequence column for parser shapes where sequence and modification tokens are separate. | Remove nullable field by making it required only on the parser-specific model that needs it. |
| `Fragments` | `label_column` | Packed fragment identity column for column-labeled fragment exports; absent for positional older DIA-NN labels. | Remove nullable field by splitting fragment label strategies. |

### Other `| None = None` Patterns

| Model | Field | Current meaning | Proposed action |
| --- | --- | --- | --- |
| `Layer` | `categories: dict[str, int] | None` | Required only for factor-encoded layers. | Change to `dict[str, int] = Field(default_factory=dict)` and keep the non-empty validator for `encoding_mode = "factor"`. |
| `ParseRule` | `sample_name_cleanup: SampleNameCleanup | None` | Optional TOML section, only valid for wide rules. | Keep nullable. |
| `ParseRule` | `modifications: Modifications | None` | Optional TOML section, required only when computed ProForma/stripped-sequence columns need it. | Keep nullable. |
| `ParseRule` | `fragments: Fragments | None` | Optional TOML section for packed fragment expansion. Only valid for `quantification_level = "fragment"`. | Keep nullable. |

### Summary Of What Should Change

Remove nullable model shape where it hides a required TOML concept:

1. Replace `Layer.source_column` and `Layer.column_pattern` with required `source: str`.
   Model long and wide rule shapes explicitly, so `source` is always present.
2. Replace `Layer.categories: dict[str, int] | None = None` with an empty-dict default.
3. Replace parser-dependent nullable modification fields with parser-specific models.
4. Replace fragment-label nullable shape with explicit fragment-label variants.

Keep nullable fields where absence represents an optional TOML section or a real format variant:

1. `ParseRule.sample_name_cleanup`
2. `ParseRule.modifications`
3. `ParseRule.fragments`

Even these optional sections can be revisited later, but they are object-level section presence,
not nullable scalar fields.

## Target Pydantic Shape

### Layer Sources

Use one required TOML key, `source`, but model long and wide rules as different entities.

```python
class LayerBase(_Strict):
    name: str
    source: str
    encoding_mode: EncodingMode = "numeric"
    categories: dict[str, int] = Field(default_factory=dict)


class LongParseRule(ParseRuleBase):
    input_shape: Literal["long"]
    layers: list[LayerBase] = Field(min_length=1)


class WideParseRule(ParseRuleBase):
    input_shape: Literal["wide"]
    layers: list[LayerBase] = Field(min_length=1)
```

The distinction is then owned by the rule shape:

- `LongParseRule`: `layer.source` is an exact vendor column.
- `WideParseRule`: `layer.source` is a regex and must contain a `sample` named group.

Do not reintroduce `source_column: str | None` or `column_pattern: str | None`.

### Modification Parser Variants

The current `Modifications` model mixes multiple parser contracts:

```python
token_pattern: str | None = None
sequence_column: str | None = None
```

Replace that with parser-specific entities:

```python
class TokenRegexModifications(_Strict):
    parser: Literal["token_regex"]
    source_column: str
    token_pattern: str
    token_position: TokenPosition = "after_residue"
    case_sensitive: bool = False
    unknown_policy: UnknownPolicy = "preserve"
    output_column: str = "proforma_sequence"
    map: list[ModificationMapEntry] = Field(min_length=1)


class AlreadyProformaModifications(_Strict):
    parser: Literal["already_proforma"]
    source_column: str
    output_column: str = "proforma_sequence"


class SeparateModColumnModifications(_Strict):
    parser: Literal["separate_mod_column"]
    source_column: str
    sequence_column: str
    token_position: TokenPosition = "after_residue"
    case_sensitive: bool = False
    unknown_policy: UnknownPolicy = "preserve"
    output_column: str = "proforma_sequence"
    map: list[ModificationMapEntry] = Field(default_factory=list)


Modifications = Annotated[
    TokenRegexModifications
    | AlreadyProformaModifications
    | SeparateModColumnModifications,
    Field(discriminator="parser"),
]
```

This means if `sequence_column` is only meaningful for one parser, that parser gets a required
`sequence_column: str`; other parsers do not have the field at all.

### Fragment Label Variants

The current fragment model has:

```python
label_column: str | None = None
```

Replace this with explicit variants:

```python
class PositionalFragments(_Strict):
    label_strategy: Literal["positional"]
    value_columns: list[str] = Field(min_length=1)
    delimiter: str = ";"
    label_output: str = "fragment_label"


class ColumnLabeledFragments(_Strict):
    label_strategy: Literal["column"]
    value_columns: list[str] = Field(min_length=1)
    label_column: str
    delimiter: str = ";"
    label_output: str = "fragment_label"


Fragments = Annotated[
    PositionalFragments | ColumnLabeledFragments,
    Field(discriminator="label_strategy"),
]
```

This makes positional DIA-NN fragments explicit instead of encoding them as
`label_column = None`.

## Target TOML

Long-format layer:

```toml
input_shape = "long"

[[layers]]
name = "FG_Quantity"
source = "FG.Quantity"
```

Wide-format layer:

```toml
input_shape = "wide"

[[layers]]
name = "Intensity"
source = "^(?P<sample>.+) Intensity$"
```

## Implementation Plan

1. Update `Layer` in `src/anndata_proteomics/rules/schema.py`:

   ```python
   class Layer(_Strict):
       name: str
       source: str
       encoding_mode: EncodingMode = "numeric"
       categories: dict[str, int] = Field(default_factory=dict)
   ```

2. Remove `source_column` and `column_pattern` from the pydantic model.

3. Split parser-dependent models:

   - modification parser variants instead of `token_pattern: str | None`
   - fragment label variants instead of `label_column: str | None`

4. Update validation:

   - long rules: `source` is interpreted as an exact vendor column
   - wide rules: `source` is compiled as a regex and must contain `(?P<sample>...)`
   - keep `axis.x_layer` and factor-category validators unchanged

5. Update converters and recognizers:

   - long conversion pivots by `layer.source`
   - wide conversion matches headers with `layer.source`
   - recognition checks exact columns for long rules and regex matches for wide rules

6. Migrate all packaged TOMLs:

   - replace each `source_column = "..."` with `source = "..."`
   - replace each `column_pattern = "..."` with `source = "..."`
   - keep `[modifications].source_column` unchanged; this is not a layer source
   - add `label_strategy = "positional"` or `label_strategy = "column"` to `[fragments]`
   - ensure `[modifications]` sections match the parser-specific shape

7. Regenerate JSON schema:

   ```bash
   uv run anndata-proteomics export-schema
   ```

8. Update docs:

   - `docs/toml_schema.md`
   - `docs/parsing_architecture.md`
   - any TODO/HOWTO files that mention `layers.source_column` or
     `layers.column_pattern`

9. Update tests:

   - model tests for missing `source`
   - long rule validation uses `source`
   - wide rule validation uses `source`
   - packaged TOML validation
   - converter/recognizer tests

## Non-Goals

- Do not rename `[modifications].source_column`; that field really is an exact
  vendor column carrying modification tokens.
- Do not add per-layer `long = true` or `source_kind`; `input_shape` is already
  the rule-level discriminator.
- Do not split into `LongLayer` and `WideLayer` unless the single `source`
  approach proves insufficient.

## Validation Commands

```bash
uv run pytest tests/test_rule_models.py tests/test_rule_loader.py tests/test_packaged_rules.py tests/test_json_schema_validation.py
uv run pytest tests/test_converters_long.py tests/test_converters_wide.py tests/test_recognize.py
uv run ruff check src/anndata_proteomics/rules/schema.py src/anndata_proteomics/converters tests
```
