# JSON Schema Reference for Parsing Rules

APB stores one self-contained `rules.json` document per software-version group under
`src/anndata_proteomics/parsing_rules/<vendor>/`. Every document has the same source
shape: metadata, one shared `base`, and a `levels` map. The map keys define the
quantification levels; fragments do not repeat `quantification_level`.

Two Pydantic models in `rules/schema.py` are authoritative:

- `ParseRuleDocument` validates the stored source document.
- `ParseRule` validates each effective rule after a level is merged over `base`.

Their generated schemas are `parsing_rules/_schema/parse_rule_document.schema.json`
and `parse_rule.schema.json`. Run `apb validate path/to/rules.json` for one document,
or `apb validate` for all packaged documents.

## Source-document shape

This minimal long-format document contains one ion level:

```json
{
  "schema_version": "0.1",
  "file_version": "1",
  "software_name": "MyTool",
  "software_version": "^1\\..*",
  "base": {
    "input_shape": "long",
    "axis": {
      "obs_keys": ["Run"],
      "duplicates": {"mode": "error"}
    },
    "columns": {
      "obs": {"select": {"Run": "R.FileName"}}
    }
  },
  "levels": {
    "ion": {
      "axis": {
        "var_keys": ["Precursor_Id"],
        "x_layer": "Intensity"
      },
      "columns": {
        "var": {"select": {"Precursor_Id": "PEP.Id"}}
      },
      "layers": [
        {"name": "Intensity", "source": "PEP.Quantity"}
      ]
    }
  }
}
```

Document fields:

| Field | Meaning |
|---|---|
| `schema_version` | Version of the source/effective configuration schema. |
| `file_version` | Revision of this software-version document. |
| `software_name` | Human-readable software name. |
| `software_version` | Regex matched against the version parsed from parameters. |
| `base` | Partial rule body shared by every level in the document. |
| `levels` | Non-empty map from level name to its partial rule body. |

Allowed level keys are `ion`, `fragment`, `peptidoform`, `peptide`, and `protein`.
Unknown keys are rejected throughout the document.

Every vendor uses this shape, including single-level vendors. Version coverage remains
explicit: DIA-NN, for example, has `v1/rules.json` and `v2/rules.json`; FragPipe has one
vendor-root `rules.json`. A document contains no paths to other configuration files.

## Base-to-level merge

APB merges the selected `levels.<level>` fragment over `base`, injects the document
metadata and level name, then validates the result as `ParseRule`.

- Objects deep-merge; level keys win.
- Scalars replace base values.
- Scalar arrays such as `obs_keys` and `var_keys` replace base arrays.
- Arrays of objects such as `layers`, `compute`, and modification `map` append in
  base-to-level order.

Put a field in `base` only when it applies to every level in that document. Typical
base content is `input_shape`, observation identity, common selected columns, sample
cleanup, and modification parsing. Feature identity, quantitative layers, and fragment
handling usually belong to a level.

Editing `base` requires every effective level to remain valid. Editing a level validates
that merged level. `apb validate` always checks the source document and every level.

## Search-parameter axis overrides

A level may conditionally override its partial `axis` using validated search parameters:

```json
{
  "levels": {
    "ion": {
      "axis": {
        "var_keys": ["ProForma_ion"],
        "x_layer": "Precursor_Normalised"
      },
      "search_parameter_overrides": [
        {
          "when_search_parameters": {
            "acquisition_method": "DDA"
          },
          "axis": {
            "x_layer": "Ms1_Normalised"
          }
        }
      ]
    }
  }
}
```

Conditions use equality over fields declared by the typed `Parameters` model. Unknown
field names and values that do not validate as the corresponding parameter type are
rejected. APB merges `base`, then the selected level, then every matching axis override
in source order. The result is validated and handed to converters as an ordinary flat
`ParseRule`; the conditional source structure is not stored in the effective rule.

Overrides are level-only: `base` is a plain `RuleFragment`, while level fragments own
`search_parameter_overrides`. Override bodies intentionally contain only a partial
`axis`. Object arrays such as layers and computed columns have append semantics in the
base-to-level merger and therefore cannot safely update an existing named declaration.

`apb validate` checks the default rule and every compatible override combination.
Ordinary packaged conversion parses search parameters before materialization. Explicit
`--rule-config` conversion also applies overrides when `--params` is supplied and uses
the level default when it is not.

## Effective rule shape

After merging, every effective `ParseRule` has these required fields:

| Field | Meaning |
|---|---|
| `schema_version` | Copied from the document. |
| `file_version` | Copied from the document. |
| `software_name` | Copied from the document. |
| `software_version` | Copied from the document. |
| `input_shape` | `long` or `wide`. |
| `quantification_level` | Injected from the selected `levels` key. |
| `axis` | Observation/feature identity and primary layer. |
| `columns` | Columns selected or computed into `obs` and `var`. |
| `column_roles` | Optional semantic locations needed by canonical-data consumers. |
| `layers` | At least one vendor-reported measurement. |

## Axis

`axis.obs_keys` names selected observation columns that identify one run.
`axis.var_keys` names selected or computed columns that identify one feature.
`axis.x_layer` must match a declared layer name.

`axis.duplicates.mode` is one of `error` (default), `aggregate`, `keep_first`, or
`keep_all_as_raw_table`. `error` rejects repeated observation-feature cells,
`aggregate` sums their layer values, and `keep_first` retains the first non-null
layer value in input order. `keep_all_as_raw_table` is reserved but not implemented.
Duplicate aggregation resolves repeated observations of an already valid level; it
does not derive another quantitative level.

## Columns

`columns.obs.select` and `columns.var.select` map APB output names to exact vendor
input columns:

```json
{
  "obs": {"select": {"R_FileName": "R.FileName"}},
  "var": {
    "select": {
      "FG_Charge": "FG.Charge",
      "EG_ModifiedSequence": "EG.ModifiedSequence"
    }
  }
}
```

Left-hand names preserve the vendor's words, case, and namespace while replacing
separators with underscores. Do not replace them with cross-vendor semantic aliases.

`columns.var.compute` is ordered. Each object has `name`, `from`, and `how`:

| `how` | Required `name` | Purpose |
|---|---|---|
| `coalesce` | Any declared output | First non-null source value in declared order. |
| `join_nonempty` | Any declared output | Non-null/non-empty source values joined with the required `separator`. |
| `stripped_sequence` | `ProForma_peptide` | Sequence without modifications. |
| `proforma_sequence` | `ProForma_peptidoform` | Normalized modified sequence. |
| `proforma_ion` | `ProForma_ion` | Peptidoform plus charge. |
| `proforma_fragment` | `ProForma_fragment` | Ion plus fragment label. |

All computed columns may depend on earlier computed columns. `coalesce` and
`join_nonempty` require at least two sources; `separator` is required only for
`join_nonempty`. A generic compute may intentionally replace a selected output with
the same name. ProForma computed columns retain their reserved APB identifier names,
not input-column aliases.

`column_roles` identifies the APB output column that carries a downstream
semantic role without repeating the vendor input name. Its currently supported
field is `protein_accessions`, which must name a declared `var` column:

```json
{
  "column_roles": {
    "protein_accessions": "Protein_Ids"
  }
}
```

Declare the role in the level that owns that column. Consumers resolve it from
the stored effective rule and continue to use `X` as the quantitative matrix;
they do not parse the vendor table again.

## Layers

A layer is one reported measurement varying over `obs × var`:

```json
{
  "name": "FG_Quantity",
  "source": "FG.Quantity",
  "encoding_mode": "numeric",
  "missing_values": [0],
  "required": false
}
```

- In a long rule, `source` is an exact vendor column.
- In a wide rule, `source` is a regex over measurement headers and must contain a
  named `(?P<sample>...)` group.
- `encoding_mode` defaults to `numeric`. String-valued matrix data uses `factor`
  with a non-empty `categories` map.
- `missing_values` lists numeric vendor sentinels that APB replaces with `NaN`
  before matrix assembly. It is layer-specific and must not include valid zeros or
  factor category codes.
- `required` defaults to false, but `axis.x_layer` is always required.

Store a value as a layer only when it varies by sample for the same feature.
Feature-invariant values belong in `var`; sample-invariant values belong in `obs`.

## Wide rules

Wide documents normally put the synthetic sample column in `base`:

```json
{
  "input_shape": "wide",
  "axis": {"obs_keys": ["sample"]},
  "columns": {"obs": {"select": {"sample": "<sample>"}}},
  "sample_name_cleanup": {"pattern": ""}
}
```

Each wide level then declares feature columns and layer regexes such as
`"^(?P<sample>.+) Intensity$"`. `sample_name_cleanup` is optional and forbidden on
long rules.

## Modification normalization

The optional `modifications` object converts vendor tokens to accessions:

```json
{
  "source_column": "Modified.Sequence",
  "parser": "token_regex",
  "token_pattern": "\\(([^()]*)\\)",
  "token_position": "after_residue",
  "case_sensitive": false,
  "unknown_policy": "preserve",
  "output_column": "proforma_sequence",
  "map": [
    {"token": "UniMod:35", "accession": "UNIMOD:35"}
  ]
}
```

`parser` must be `token_regex`; `map` must contain at least one entry. Canonical
metadata comes from `modifications/unimod_registry.toml`. That internal reference
dataset remains TOML and is not a parsing-rule configuration.

## Fragment rules

Packed fragment exports declare a `fragments` object with one label strategy:

- `column` requires `label_column` and reads labels from packed vendor tokens.
- `positional` synthesizes `frag_0`, `frag_1`, and so on.

Both require `value_columns` and may set `delimiter` and `label_output`. Fragment
levels compute `ProForma_fragment` from `ProForma_ion` plus `label_output` and use it
in `axis.var_keys`.

## Locations and version selection

```text
parsing_rules/
  diann/
    v1/rules.json
    v2/rules.json
  fragpipe/rules.json
  spectronaut/rules.json
```

The registry discovers vendor-root `rules.json` and `v*/rules.json`. The required
`software_version` regex selects the matching document. The folder is organization,
not an additional fallback rule; a version mismatch is an error.

`apb convert data.tsv LEVEL --rule-config rules.json` selects one level from an
external document and writes AnnData. Without `LEVEL`, APB converts every level whose
required columns match the data and always writes MuData, including when only one
level matches. Supplying `--params` materializes any matching search-parameter axis
overrides before column matching and conversion.

## Adding a vendor, version, or level

1. Copy the closest `rules.json` document.
2. Keep the existing version grouping unless the supported version coverage truly changes.
3. Update document metadata and place common behavior in `base`.
4. Add each real quantitative level under `levels`.
5. Map exact vendor columns into `obs`, `var`, and layers.
6. Configure ProForma computations and modification mappings where sequences exist.
7. Run `apb validate`, then the packaged-rule, recognition, and conversion tests.
