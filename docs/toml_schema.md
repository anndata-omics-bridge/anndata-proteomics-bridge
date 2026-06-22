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

Four reserved compute names mirror the quantification levels:

| Compute name | `how` | Meaning |
|---|---|---|
| `ProForma_peptide` | `stripped_sequence` | bare sequence, no mods |
| `ProForma_peptidoform` | `proforma_sequence` | sequence + mods |
| `ProForma_ion` | `proforma_ion` | peptidoform + `/charge` |
| `ProForma_fragment` | `proforma_fragment` | ion + `/fragment_label` |

How rules use them:
- **Peptidoform-level rules** (`quantification_level = "peptidoform"`):
  `var_keys = ["ProForma_peptidoform"]` produced by
  `how = "proforma_sequence"`. Optionally also expose
  `ProForma_peptide`.
- **Ion-level rules** (`quantification_level = "ion"`):
  `var_keys = ["ProForma_ion"]` produced by `how = "proforma_ion"`,
  chained from a `ProForma_peptidoform` intermediate. Optionally also
  expose `ProForma_peptide`.
- **Fragment-level rules** (`quantification_level = "fragment"`):
  `var_keys = ["ProForma_fragment"]` produced by `how = "proforma_fragment"`,
  chained from a `ProForma_ion` intermediate (so here `proforma_ion` is an
  *intermediate*, not the var key). Requires a `[fragments]` block (see below).
  Grammar: `ProForma_fragment = "{peptidoform}/{charge}/{fragment_label}"`, e.g.
  `M[UNIMOD:35]PEPTIDE/2/b4-unknown^1` — note the `/` separator carries charge
  after the peptidoform and the fragment label after the ion.

**Protein-level rules** (`quantification_level = "protein"`) use a plain vendor
column as `var_keys` (e.g. `Protein_Group`) with no ProForma compute.

Why `ProForma_peptide` is computed from the modified-sequence column
rather than the vendor's "peptide" column: stripping the modification
tokens with one controlled algorithm gives a consistent result across
vendors. Vendor "peptide" columns disagree (case, flanking residues,
presence-or-absence), so deriving from the modification-bearing column
keeps a single source of truth.

Schema invariants worth knowing:
- The compute `name` is pinned by `how`: `stripped_sequence` →
  `ProForma_peptide`, `proforma_sequence` → `ProForma_peptidoform`,
  `proforma_ion` → `ProForma_ion`. The validator rejects any other
  name.
- At **ion** level a `how = "proforma_ion"` compute must appear in
  `axis.var_keys`; at **fragment** level it is an intermediate and must
  *not* be a var key (the var key is `ProForma_fragment`).
- `how = "proforma_ion"` requires exactly two source columns
  (peptidoform intermediate + charge).
- `how = "proforma_fragment"` is fragment-level only, requires exactly two
  source columns (a `ProForma_ion` intermediate + the `[fragments].label_output`
  column), and must appear in `axis.var_keys`.
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

- `parse_<software>_<quantification_level>.toml` inside a vendor (and, for version-specific
  DIA-NN levels, a `vN`/`vN_M` subfolder — see "Version folders" below). Single-version vendors may
  keep a legacy flat `parse_<software>_<quantification_level>_<n>.toml`.
- The folder and filename use the lowercase vendor short-name (`diann/parse_diann_ion.toml`,
  `diann/v1/parse_diann_protein.toml`); the `software_name` value inside the TOML preserves the
  canonical spelling (`"DIA-NN"`).
- `quantification_level` in the filename must match the in-TOML `quantification_level` value.
  `tests/test_packaged_rules.py` enforces this.

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
- `quantification_level` — `"ion"`, `"peptidoform"`, `"peptide"`,
  `"protein"`, or `"fragment"`. Must match the filename token.
  `"fragment"` requires a `[fragments]` block.

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
    `"proforma_ion"`, `"proforma_fragment"`. The first two require a
    `[modifications]` block and exactly one source column. `proforma_ion`
    requires exactly two source columns and is valid at ion level (where
    its `name` must appear in `axis.var_keys`) and fragment level (as an
    intermediate). `proforma_fragment` is fragment-level only, requires
    exactly two source columns (a `ProForma_ion` intermediate +
    `[fragments].label_output`), and its `name` must appear in
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

### Fragments

Only for `quantification_level = "fragment"`. Some vendors (DIA-NN) do not emit one
row per fragment; instead they pack per-fragment values as parallel, delimiter-joined
lists inside each precursor row (`Fragment.Info`, `Fragment.Quant.Raw`, …, aligned by
index and often terminated by a trailing delimiter). The `[fragments]` block tells
`convert()` to **explode** those lists into one row per fragment *before* column
materialization and the long-format pivot, so the rest of the pipeline is reused
unchanged.

- `[fragments]`
  - `value_columns` — array of the parallel packed value columns to split
    (e.g. `["Fragment.Quant.Raw", "Fragment.Correlations"]`). All same length per row; mismatch raises.
  - `label_column` — *optional* packed column whose tokens identify each fragment
    (e.g. `"Fragment.Info"`, tokens like `b4-unknown^1/327.16`). When **omitted** (older DIA-NN
    with no `Fragment.Info`), labels are **positional**: `frag_0`, `frag_1`, … by index within the
    precursor.
  - `delimiter` — token separator, default `";"`.
  - `label_output` — name of the column the explode produces (the token before `/` of
    `label_column`, or `frag_<i>` when positional), default `"fragment_label"`. Use it as a source
    for a `how = "proforma_fragment"` compute.

Each packed value column appears **twice** on purpose: in `[fragments].value_columns`
(so explode knows to split it) and in `[[layers]]` with a matching `source_column` (so
the now-scalar column is pivoted into a layer like any other long column).

**Caveats.** Explode multiplies the row count by the fragments-per-precursor (~12 for
DIA-NN), so converting a full report into a dense fragment matrix is memory-heavy: a
single 6-run AIF report builds ~827k features and peaks around ~6.5 GB (the matrix is
~90% dense, so this is genuinely large data, not overhead — the converter scatters
directly into the dense matrix rather than via `pivot_table`, and trims unused columns
before exploding). For large reports, prefer converting per run / filtering precursors
first; a chunked-streaming builder would be the next step for full-scale fragment work.
DIA-NN fragment columns also vary by version/config — some exports drop `Fragment.Info`
or carry a reduced set of `Fragment.Quant.*` columns, so a fragment rule may not fit
every DIA-NN file.

### Multiple levels per vendor

A vendor's single export can back several levels (DIA-NN's `report.tsv` backs ion, peptidoform,
peptide, protein, and fragment). Ship one TOML per level. Because every level reads the same
columns, header-based `recognize()` cannot pick a *level* and returns `None` for such a file; the
level is selected explicitly.

### Version folders (DIA-NN columns change across versions)

DIA-NN's `report.tsv` columns differ by version: `Fragment.Quant.Corrected` is 1.7–1.8 only;
`Fragment.Quant.Raw` / `Fragment.Correlations` and `PG.Normalised` / `PG.Quantity` exist ≤1.9.2 and
are gone in 2.x; `Fragment.Info` is absent from all pure DIA-NN exports. So a single rule per level
cannot fit every version. Version-dependent levels live in **version subfolders** selected by the
software version:

```text
parsing_rules/diann/
  parse_diann_ion.toml          # version-agnostic levels live at the vendor root
  parse_diann_peptidoform.toml
  parse_diann_peptide.toml
  v1/ parse_diann_protein.toml  # 1.x: PG.Normalised/PG.Quantity present
      parse_diann_fragment.toml #      positional fragment (no Fragment.Info)
  v2/ parse_diann_protein.toml  # 2.x: PG.MaxLFQ + Genes.MaxLFQ only
```

- **The folder name is the selector.** `vN` covers major version `N`; finer `vN_M` (e.g. `v1_9`)
  covers `N.M`. `registry.resolve_rule_path(software, level, version)` picks the most-specific
  folder whose version prefixes the file's version and contains `parse_<sw>_<level>.toml`, else the
  vendor-root file (so version-agnostic levels and single-version vendors keep working). `None` when
  a version doesn't provide a level (e.g. fragment on 2.x).
- **The version comes from the param file**, which is **mandatory**: the GUI / `convert_one` parse
  the co-located param (`params.registry.parse_params`) for `software_version`, resolve the rule by
  version, then **validate** the data columns against it (`converters.recognize.matches`). A column
  mismatch is an **error** (verify the version / param file) — no silent fallback.

**Adding a variant** for a new DIA-NN format: drop `parse_diann_<level>.toml` into the right
`vN`/`vN_M` folder (create it if needed) with that version's columns. No code change — the resolver
finds it. `load_packaged_rule(software, level, version)` / `find_rule(...)` address a specific
version; pass `None` for version-agnostic root rules.

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
var_keys = ["ProForma_ion"]
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
name = "ProForma_peptidoform"
from = ["Modified_Sequence"]
how = "proforma_sequence"

[[columns.var.compute]]
name = "ProForma_ion"
from = ["ProForma_peptidoform", "Precursor_Charge"]
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
var_keys = ["ProForma_ion"]
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
name = "ProForma_peptidoform"
from = ["Modified_Sequence"]
how = "proforma_sequence"

[[columns.var.compute]]
name = "ProForma_ion"
from = ["ProForma_peptidoform", "Charge"]
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

DIA-NN ships all five levels from one `report.tsv` —
`parse_diann_{ion,peptidoform,peptide,protein,fragment}_1.toml` — demonstrating the
one-TOML-per-level pattern.

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
