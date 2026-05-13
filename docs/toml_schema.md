# TOML Schema Reference for Parsing Rules

This document is the contract reference for the parsing-rule TOML files
shipped under `src/anndata_proteomics/parsing_rules/<vendor>/`. Each
file maps one vendor quantification export into AnnData. The pydantic
validator lives at `src/anndata_proteomics/rules/schema.py` and is
authoritative — when this doc and the validator disagree, the validator
wins.

## Background concepts

### AnnData orientation

- `obs` = samples / runs
- `var` = quantified features
- `X` = primary quantitative matrix
- `layers` = additional per-sample-per-feature matrices

(`uns` is written by the parser at conversion time — provenance, the
serialized rule, parameter records, factor mappings — and is not
declared in the TOML.)

### Axis keys

- `axis.obs_keys` defines the columns that become the AnnData
  observation axis.
- `axis.var_keys` defines the columns that become the AnnData variable
  axis.
- Pragmatic rule for `var_keys`: use the smallest set of declared output
  column names that avoids duplicates in `var`.

### Obs-axis conventions

- **Wide rules**: `obs_keys = ["sample"]` and `sample = "<sample>"`
  under `[columns.obs.select]`. The `sample` token is captured by the
  `(?P<sample>...)` group of each layer's `column_pattern`.
- **Long rules**: `obs_keys` is the single column that uniquely
  identifies a run within the vendor's export. Preserve the vendor's
  natural name on the LHS — DIA-NN uses `Run`, MaxQuant `Raw_File`,
  Spectronaut `R_FileName`. These are vendor-specific identifiers, not
  synonyms for "sample"; no canonical run-id name is imposed.
- Additional obs-side annotations (`Experiment`, `Fraction`,
  `R_Condition`, …) may appear in `[columns.obs.select]` even when not
  in `axis.obs_keys`. They enrich `adata.obs` but do not participate in
  uniqueness.

### Var-axis naming convention

Per the [HUPO-PSI ProForma 2.0 spec](https://github.com/HUPO-PSI/ProForma),
ProForma covers three quantification levels:
- bare peptide `PEPTIDE` — a degenerate ProForma with no modifications;
- peptidoform `M[UNIMOD:35]PEPTIDE` — sequence + mods;
- ion `M[UNIMOD:35]PEPTIDE/2` — peptidoform + charge (spec §7.1,
  optional extension).

Three reserved compute names mirror those three levels:

| Compute name | `how` | Meaning |
|---|---|---|
| `ProForma_peptide` | `stripped_sequence` | bare sequence, no mods |
| `Peptidoform` | `proforma_sequence` | sequence + mods |
| `ProForma` | `proforma_ion` | peptidoform + `/charge` |

How rules use them:
- **Peptidoform-level rules** (`quantification_level = "peptidoform"`):
  `var_keys = ["ProForma"]` produced by `how = "proforma_sequence"`.
  Optionally also expose `ProForma_peptide`.
- **Ion-level rules** (`quantification_level = "ion"`):
  `var_keys = ["ProForma"]` produced by `how = "proforma_ion"`. The
  ion compute chains through a `Peptidoform` intermediate (local to
  the TOML); `ProForma_peptide` may be exposed alongside.

Why `ProForma_peptide` is computed from the modified-sequence column
rather than the vendor's "peptide" column: stripping the modification
tokens with one controlled algorithm gives a consistent result across
vendors. Vendor "peptide" columns disagree (case, flanking residues,
presence-or-absence), so deriving from the modification-bearing column
keeps a single source of truth.

Schema invariants worth knowing:
- Any `how = "proforma_ion"` compute must appear in `axis.var_keys`
  (`schema.py:198-206`).
- `how = "proforma_ion"` requires exactly two source columns
  (peptidoform intermediate + charge).
- `how = "proforma_sequence"` and `how = "stripped_sequence"` each
  require exactly one source column and a `[modifications]` block.

### Column naming

- Right-hand-side names must preserve the exact vendor column names.
- Left-hand-side names may be cleaner internal names.
- `[columns.*.select]` is strictly for values that exist in the input
  table, plus the wide-file placeholder `"<sample>"`. APB-derived
  values (the column whose name matches `modifications.output_column`,
  and the reserved literal `stripped_sequence`) must never appear under
  `select`; declare them via `[[columns.var.compute]]`.

### Numeric vs string layers

`layers` are numeric storage. String-valued matrix fields (e.g.
FragPipe `Match Type`) are encoded as integer factors via
`encoding_mode = "factor"` + a `categories` mapping. The factor
mapping is stored both in the TOML and in `uns` at conversion time.

### Version fields

- `schema_version` — version of the TOML schema itself.
- `file_version` — version of this specific parsing-rule TOML.
- `software_version` — vendor software version metadata, when known.

Filename convention:

- `parse_<software>_<quantification_level>_<file_version>.toml`
- The folder and filename use the lowercase vendor short-name
  (`diann/parse_diann_ion_1.toml`); the `software_name` value inside
  the TOML preserves the canonical spelling (`"DIA-NN"`).
- `quantification_level` in the filename must match the in-TOML
  `quantification_level` value. `tests/test_packaged_rules.py`
  enforces this.

### Duplicate handling

Each TOML defines `[axis.duplicates] mode` for duplicate `(obs, var)`
entries. Allowed values: `error`, `aggregate`, `keep_first`,
`keep_all_as_raw_table`. The schema default is `error`; override only
when a software-specific aggregation rule is intentionally needed.

### Long vs wide

Long and wide rules share the same top-level concepts. The only
difference is how a layer finds its source data:

- long: `layers.source_column`
- wide: `layers.column_pattern`

Wide rules still must define `obs`. In the minimal case, `obs` is
created from the `<sample>` token extracted by the layer regex.
Vendor-derived obs names are kept as they are at this stage; richer
obs annotation belongs to a later stage.

## TOML schema

This section enumerates every section / key the validator recognises.

### Common entries

These entries are valid for both long and wide rules.

- `schema_version` — string. Example: `"0.1"`.
- `file_version` — string. Example: `"1"`.
- `software_name` — string, human-readable software identifier.
  Example: `"FragPipe"`.
- `software_version` — string, optional.
- `input_shape` — `"long"` or `"wide"`.
- `quantification_level` — `"ion"`, `"peptidoform"`, `"peptide"`, or
  `"protein"`. Must match the filename token.

- `[axis]`
  - `obs_keys` — array of strings. Must reference declared output
    column names from `[columns.obs.select]`. (Computed obs columns
    are not supported.)
  - `var_keys` — array of strings. Must reference declared output
    column names from `[columns.var.select]` or
    `[[columns.var.compute]]`. The schema requires that any
    `how = "proforma_ion"` compute appears in `var_keys`.
  - `x_layer` — string. Must match one `layers.name`.
- `[axis.duplicates]`
  - `mode` — `"error"` (default), `"aggregate"`, `"keep_first"`,
    `"keep_all_as_raw_table"`.

- `[columns.obs.select]` — key-value mapping
  `Internal_Name = "Vendor column"`, or
  `Internal_Name = "<sample>"` for wide rules.

- `[columns.var.select]` — key-value mapping
  `Internal_Name = "Vendor column"`. Values must be original
  input-table columns; reserved derived names
  (`modifications.output_column` and `stripped_sequence`) must not
  appear here.

- `[[columns.var.compute]]` — APB-derived var columns
  - `name` — output column name.
  - `from` — list of declared var column names (from `select` or
    earlier `compute` entries).
  - `how` — one of `"proforma_sequence"`, `"stripped_sequence"`,
    `"proforma_ion"`. The first two require a `[modifications]` block
    and exactly one source column. `proforma_ion` is ion-level only,
    requires exactly two source columns, and its `name` must appear in
    `axis.var_keys`.

- `[[layers]]`
  - `name` — internal layer name (required).
  - `encoding_mode` — `"numeric"` (default) or `"factor"`.
  - `categories` — `{ "value" = code, … }`. Required and non-empty
    when `encoding_mode = "factor"` (enforced by `schema.py`).
  - `source_column` — required for long, forbidden for wide.
  - `column_pattern` — required for wide, forbidden for long.

### Long-only entries

- `layers.source_column` — vendor column to pivot into this layer.

### Wide-only entries

- `layers.column_pattern` — regex that finds the source columns for
  this layer. Must expose a named capture group `(?P<sample>...)`
  when sample names are encoded in the column names. The `sample`
  capture describes the *vendor column shape*, not the *user's sample
  naming* — use `.+` for the sample token. If you need to strip or
  rewrite the captured sample names, use `[sample_name_cleanup]`.

- `[sample_name_cleanup]` — optional, wide-only (forbidden on long).
  - `pattern` — regex used to extract or normalize a basename /
    raw-file-like obs name from the captured `<sample>` tokens.

### Modifications

Every TOML carries a `[modifications]` block. It turns embedded
modification tokens in a vendor sequence column (e.g.
`Modified.Sequence`) into a ProForma-normalised output column.

- `[modifications]`
  - `source_column` — vendor column containing modification tokens.
  - `parser` — `"token_regex"` (default), `"already_proforma"`, or
    `"separate_mod_column"`. Each parser has consistency constraints
    enforced by the validator:
    - `token_regex` requires `token_pattern` AND at least one
      `[[modifications.map]]` entry.
    - `already_proforma` forbids both `token_pattern` and `map`.
    - `separate_mod_column` requires `source_column` and accepts an
      optional `map`.
  - `token_pattern` — regex whose first capture group is the vendor
    token. Required for `token_regex`. Example: `"\\(([^()]*)\\)"`.
  - `token_position` — where the token sits relative to its residue.
    One of `"before_residue"`, `"after_residue"` (default), `"n_term"`,
    `"c_term"`, `"embedded"`, `"unknown"`.
  - `case_sensitive` — bool, default `false`.
  - `unknown_policy` — what to do with tokens not in `map`. One of
    `"preserve"` (default), `"drop"`, `"error"`.
  - `sequence_column` — used only by `separate_mod_column`; the
    column carrying the stripped (unmodified) sequence.
  - `output_column` — name of the derived column that
    `how = "proforma_sequence"` exposes. Default
    `"proforma_sequence"`.
- `[[modifications.map]]` — one entry per vendor token.
  - `token` — the vendor's token string as it appears between the
    `token_pattern` delimiters.
  - `accession` — UNIMOD/MOD accession (`UNIMOD:35`, `MOD:00425`).
    The validator requires the `UNIMOD:N` or `MOD:N` shape.

The canonical fields a downstream consumer expects (`name`, `target`,
`position`, `mass_delta`) are NOT carried per-tool. They are filled at
rule-load time from
`src/anndata_proteomics/modifications/unimod_registry.toml` so all
tools agree on what e.g. `UNIMOD:35` means. Adding a new accession
requires adding it to that registry first.

Vendor tokens may be numeric mass deltas (`"15.9949"`, `"+57.02"`),
named labels (`"Acetyl (Protein N-term)"`), or UniMod-style strings
(`"UniMod:35"`) — see the worked examples below.

## Long example (DIA-NN)

```toml
schema_version = "0.1"
file_version = "1"
software_name = "DIA-NN"
software_version = "2.3.0"
input_shape = "long"
quantification_level = "ion"

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
Precursor_Charge = "Precursor.Charge"
Protein_Ids = "Protein.Ids"
Genes = "Genes"

[[columns.var.compute]]
name = "Peptidoform"
from = ["Modified_Sequence"]
how = "proforma_sequence"

[[columns.var.compute]]
name = "ProForma"
from = ["Peptidoform", "Precursor_Charge"]
how = "proforma_ion"

[[layers]]
name = "Precursor_Normalised"
source_column = "Precursor.Normalised"

[[layers]]
name = "Q_Value"
source_column = "Q.Value"

[modifications]
source_column = "Modified.Sequence"
parser = "token_regex"
token_pattern = "\\(([^()]*)\\)"
token_position = "after_residue"
case_sensitive = false
unknown_policy = "preserve"
output_column = "proforma_sequence"

[[modifications.map]]
token = "UniMod:1"
accession = "UNIMOD:1"

[[modifications.map]]
token = "UniMod:35"
accession = "UNIMOD:35"
```

## Wide example (FragPipe)

```toml
schema_version = "0.1"
file_version = "1"
software_name = "FragPipe"
software_version = "23.0"
input_shape = "wide"
quantification_level = "ion"

[axis]
obs_keys = ["sample"]
var_keys = ["ProForma"]
x_layer = "Intensity"

[axis.duplicates]
mode = "error"

[columns.obs.select]
sample = "<sample>"

[columns.var.select]
Peptide_Sequence = "Peptide Sequence"
Modified_Sequence = "Modified Sequence"
Charge = "Charge"
Protein_ID = "Protein ID"
Gene = "Gene"

[[columns.var.compute]]
name = "Peptidoform"
from = ["Modified_Sequence"]
how = "proforma_sequence"

[[columns.var.compute]]
name = "ProForma"
from = ["Peptidoform", "Charge"]
how = "proforma_ion"

[[layers]]
name = "Intensity"
column_pattern = "^(?P<sample>.+) Intensity$"

[[layers]]
name = "Spectral_Count"
column_pattern = "^(?P<sample>.+) Spectral Count$"

[[layers]]
name = "Match_Type"
column_pattern = "^(?P<sample>.+) Match Type$"
encoding_mode = "factor"
categories = { "unmatched" = 0, "MS/MS" = 1, "MBR" = 2 }

[modifications]
source_column = "Modified Sequence"
parser = "token_regex"
token_pattern = "\\[([^\\]]+)\\]"
token_position = "after_residue"
case_sensitive = false
unknown_policy = "preserve"
output_column = "proforma_sequence"

[[modifications.map]]
token = "57.0215"
accession = "UNIMOD:4"

[[modifications.map]]
token = "15.9949"
accession = "UNIMOD:35"
```

## Software families already shipped

- Long: DIA-NN, Spectronaut, MaxQuant (evidence-like).
- Wide: FragPipe, PEAKS, WOMBAT.

This is why the rule schema supports both `source_column` and
`column_pattern`.

## Adding a new vendor

1. Copy the closest-matching shipped TOML into a new
   `parsing_rules/<vendor>/parse_<vendor>_<level>_1.toml`.
2. Update `software_name`, `software_version`, `input_shape`,
   `quantification_level`.
3. Replace `[columns.var.select]` RHS values with the new vendor's
   actual column names; keep the LHS clean.
4. Update layer `source_column` (long) or `column_pattern` (wide).
   Keep the sample-token in wide regexes as `.+`.
5. Adjust `[modifications]`: pick the right `token_pattern` for the
   vendor's token syntax, then enumerate the vendor's tokens under
   `[[modifications.map]]`. Add any new UNIMOD accessions to
   `modifications/unimod_registry.toml`.
6. Run the packaged-rules test:
   `pytest tests/test_packaged_rules.py -k <vendor>`.
