# Review TOML Structure For Explicit Columns

## Problem

The current TOML structure makes `proforma_sequence` look like a normal input
column, but it is not present in vendor quantification tables. It is produced by
APB's modification normalization step.

This makes rules hard to reason about because `axis.var_keys` and `columns.var`
mix three different concepts:

- columns selected directly from the input table;
- normalized columns produced by APB preprocessing, e.g. `proforma_sequence`;
- derived identity columns, e.g. ProForma ion `<modified-sequence>/<charge>`.

The goal is to make the TOML explicit and crystal clear.

## Proposed Direction

Separate selected columns from computed columns.

Sketch:

```toml
[axis]
obs_keys = ["Run"]
var_keys = ["ProForma"]
x_layer = "Precursor_Normalised"

[axis.duplicates]
mode = "error"

[columns.obs.select]
Run = "Run"

[columns.var.select]
Modified_Sequence = "Modified.Sequence"
Charge = "Precursor.Charge"
Protein_Ids = "Protein.Ids"

[[columns.var.compute]]
name = "Peptidoform"
from = ["Modified_Sequence"]
how = "proforma_sequence"

[[columns.var.compute]]
name = "ProForma"
from = ["Peptidoform", "Charge"]
how = "proforma_ion"
```

Meaning:

- `columns.var.select` names original input-table columns that are copied into
  `adata.var`.
- Values in `select` must be vendor-table columns. APB-normalized columns such
  as `proforma_sequence` and `stripped_sequence` are not allowed there.
- `columns.var.compute` names columns that APB derives.
- `Peptidoform` is computed from the modification-normalized ProForma sequence.
- `ProForma` for ion rules is computed as `<Peptidoform>/<charge>`.
- `axis.var_keys = ["ProForma"]` then makes `adata.var_names` the ProForma ion
  identifier.

## Decisions

- `select` values may reference original vendor-table columns only, plus the
  wide-file `"<sample>"` placeholder for obs.
- Computed columns reference already declared output names from `select` or
  earlier `compute` entries, e.g. `from = ["Peptidoform", "Charge"]`.
- `how` is a string enum: `proforma_sequence`, `stripped_sequence`, and
  `proforma_ion`.
- Both obs and var use the symmetric `[columns.<axis>.select]` structure.
- No backward compatibility is kept for old `[columns.var]`, `[columns.obs]`, or
  top-level `[duplicates]`; new rules must use the explicit structure.

## Proposed Schema Shape

Python model sketch:

```python
class ColumnCompute(BaseModel):
    name: str
    from_: list[str] = Field(alias="from", min_length=1)
    how: Literal["proforma_sequence", "stripped_sequence", "proforma_ion"]

class Axis(BaseModel):
    obs_keys: list[str]
    var_keys: list[str]
    x_layer: str
    duplicates: Duplicates = Field(default_factory=Duplicates)

class ColumnGroup(BaseModel):
    select: dict[str, str] = Field(default_factory=dict)
    compute: list[ColumnCompute] = Field(default_factory=list)

class Columns(BaseModel):
    obs: ColumnGroup
    var: ColumnGroup
```

Possible restriction: `columns.obs.compute` may be forbidden initially if there
is no use case.

## Migration Example

Required structure:

```toml
[axis]
var_keys = ["ProForma"]

[axis.duplicates]
mode = "error"

[columns.var.select]
Precursor_Charge = "Precursor.Charge"
Protein_Ids = "Protein.Ids"
Modified_Sequence = "Modified.Sequence"

[[columns.var.compute]]
name = "Peptidoform"
from = ["Modified_Sequence"]
how = "proforma_sequence"

[[columns.var.compute]]
name = "ProForma"
from = ["Peptidoform", "Precursor_Charge"]
how = "proforma_ion"
```

## Implementation Plan After Approval

1. Update all packaged TOML rules first.
2. Update `rules/schema.py` to model `select` and `compute` explicitly.
3. Update converters to:
   - build selected columns from `columns.<axis>.select`;
   - compute derived columns before axis-index construction;
   - use `axis.var_keys` against the final explicit output columns.
4. Move duplicate handling under `axis.duplicates`.
5. Update JSON schema export.
6. Update tests for:
   - TOML parsing;
   - invalid compute definitions;
   - ProForma ion formatting;
   - packaged conversion outputs.

## Acceptance Criteria

- No TOML rule makes a computed value look like an original vendor column.
- Ion `adata.var_names` are ProForma ion identifiers: `<modpeptideseq>/<charge>`.
- Peptidoform rules remain sequence-level ProForma without `/charge`.
- Duplicate handling is expressed as `[axis.duplicates]`.
- The JSON schema reflects the new explicit structure.
