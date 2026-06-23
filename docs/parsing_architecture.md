# Parsing architecture (UML)

The parsing subsystem lives under `src/anndata_proteomics/` and is split into five
subpackages. This doc gives a top-level dependency map, then **one diagram + a short
description per module**, then the two end-to-end flows that tie them together.

Two distinct "parsing" concerns to keep separate:

- **`params/`** parses a vendor **search-parameter file** (whole-experiment settings) into a
  typed `Parameters` record.
- **`modifications/`** parses **peptide modification strings** (a single sequence's mods) and
  models searched modifications for SDRF.

---

## Module dependency overview

Arrows mean **"imports / depends on"**. `rules` and `readers` are leaves; `converters` is the
orchestrator that pulls everything together.

```mermaid
flowchart TD
    converters --> params
    converters --> modifications
    converters --> rules
    converters --> readers
    params --> modifications
    modifications --> rules
    modifications --> unimod[(unimod_registry.toml)]

    classDef leaf fill:#eef2ff,stroke:#9aa7d8;
    class rules,readers leaf;
```

| Module | One-line role |
|---|---|
| `params/` | Vendor parameter-file → typed `Parameters`. |
| `modifications/` | Vendor modified-sequence → ProForma + modification models. |
| `rules/` | The TOML parsing-rule schema (`ParseRule`) + its loader/registry. |
| `readers/` | Read a tabular quant file → `DataFrame`. |
| `converters/` | `DataFrame` + `ParseRule` → `AnnData`. |

---

## `params/` — vendor parameter-file parsing

**Parses a vendor search-parameter file into one typed `Parameters` record** — the model for
parameter-file parsing.

- **Inputs:** DIA-NN log/cfg, MaxQuant `mqpar.xml`, Sage JSON, AlphaPept / WOMBAT YAML,
  FragPipe `.workflow`, PEAKS / Spectronaut text, MSAID, MetaMorpheus.
- **Entry point:** each vendor module exposes `extract_params(source) -> Parameters`.
- **Dispatch:** `registry.py` looks up the parser by software name.
- **Reading:** `_common.py` centralizes file I/O (`read_text` / `read_lines`).
- **AnnData I/O:** `anndata_io.py` reads/writes a `Parameters` into `adata.uns`.

See [parameter_parsers.md](parameter_parsers.md) for the per-vendor breakdown (input formats,
parse techniques, and the three modification-mapping families).

```mermaid
classDiagram
    direction LR
    class _Strict {
        <<pydantic, extra=forbid>>
    }
    class Parameters {
        +str|None software_name
        +str|None software_version
        +str|None enzyme
        +bool|None enable_match_between_runs
        +int|None allowed_miscleavages
        +to_series() Series
        +from_series(series) Parameters$
    }
    class Probability {
        +float value
    }
    class MassTolerance {
        +float|None value
        +str|None unit
        +str mode
        +parse(value) MassTolerance$
    }
    class UnparsedParameter {
        +str name
        +ScalarValue value
        +str|None source
    }
    class SearchedModification {
        <<from modifications/>>
    }
    _Strict <|-- Parameters
    _Strict <|-- Probability
    _Strict <|-- MassTolerance
    _Strict <|-- UnparsedParameter
    Parameters *-- "0..3" Probability : ident_fdr
    Parameters *-- "0..2" MassTolerance : precursor / fragment
    Parameters *-- "*" UnparsedParameter : unparsed_parameters
    Parameters ..> SearchedModification : fixed_mods / variable_mods
    note for Parameters "~20 more scalar fields: charges, m/z, peptide lengths, FDR, quant method. Built by each vendor extract_params(); validators canonicalize enzyme, FDR, tolerance, mods."
```

---

## `modifications/` — modification-string parsing & models

**Normalizes peptide modifications** — modification identities/strings, distinct from `params/`
(whole-file settings). Two jobs:

- **ProForma normalization:** `apply_rules.apply_rule(seq, rule)` turns a vendor
  modified-sequence string (e.g. `PEPM(ox)TIDE`, `_[ac]PEP…`) into a canonical ProForma string
  plus localized `ModificationOccurrence`s. `pipeline.py` builds the `ModificationRule` from a
  rules `Modifications` section, resolving each `MapEntry` against the bundled `unimod_registry`.
- **Searched modifications:** `SearchedModification` models fixed/variable mods from parameter
  files, for SDRF export (`sdrf.py`).
- **Rendering:** `proforma.py` renders the ProForma string.

```mermaid
classDiagram
    direction LR
    class ModType {
        <<enum>>
        fixed
        variable
        unknown
    }
    class SearchedModification {
        +str name
        +str|None accession
        +ModType mod_type
        +str|None target
        +float|None mass_delta
    }
    class ModificationOccurrence {
        +str name
        +str|None accession
        +int|None sequence_index
        +str|None position
        +float|None mass_delta
        +str|None source_token
    }
    class ModifiedSequence {
        +str stripped_sequence
        +str proforma_sequence
        +list~str~ unknown_tokens
    }
    class MapEntry {
        <<frozen>>
        +str token
        +str name
        +str|None accession
        +float|None mass_delta
    }
    class ModificationRule {
        <<frozen>>
        +str token_pattern
        +str token_position
        +str unknown_policy
        +str output_column
    }
    class UnimodEntry {
        +str accession
        +str name
        +str target
        +float mass_delta
    }
    SearchedModification --> ModType
    ModifiedSequence *-- "*" ModificationOccurrence : occurrences
    ModificationRule *-- "*" MapEntry : entries
    note for ModificationRule "apply_rule(seq, rule) returns a ModifiedSequence. pipeline._to_runtime_rule builds the rule from a rules Modifications section, filling each MapEntry from the UnimodEntry registry."
```

---

## `rules/` — the TOML parsing-rule schema

**Defines and loads the TOML parsing-rule schema** that tells the converters how to turn a
vendor table into AnnData. A leaf subpackage (imports no other subpackage).

- `schema.py` — the pydantic `ParseRule` (axis keys, column select/compute, layers, optional
  modifications section) with cross-field validation.
- `loader.py` — parse + validate a TOML file.
- `registry.py` — find packaged rules by `(software, level, version)`.
- `validate.py` — validate rule files.

```mermaid
classDiagram
    direction LR
    class ParseRule {
        +str schema_version
        +str software_name
        +str software_version
        +str input_shape
        +str quantification_level
    }
    class Axis {
        +list~str~ obs_keys
        +list~str~ var_keys
        +str x_layer
    }
    class Duplicates
    class Columns
    class ColumnGroup {
        +dict select
    }
    class ColumnCompute {
        +str name
        +list~str~ from_
        +str how
    }
    class Layer {
        +str name
        +str source
        +str encoding_mode
        +dict categories
    }
    class Modifications {
        <<union by parser>>
        +str source_column
        +str parser
        +str output_column
    }
    class TokenRegexModifications {
        +str token_pattern
        +str token_position
        +bool case_sensitive
        +str unknown_policy
    }
    class AlreadyProformaModifications
    class SeparateModColumnModifications {
        +str sequence_column
        +str token_position
        +bool case_sensitive
        +str unknown_policy
    }
    class Fragments {
        <<union by label_strategy>>
        +str label_strategy
        +list~str~ value_columns
        +str delimiter
        +str label_output
    }
    class PositionalFragments
    class ColumnLabeledFragments {
        +str label_column
    }
    class ModificationMapEntry {
        +str token
        +str accession
    }
    class SampleNameCleanup
    ParseRule *-- Axis
    ParseRule *-- Columns
    ParseRule *-- "1..*" Layer
    ParseRule *-- "0..1" Modifications
    ParseRule *-- "0..1" Fragments
    ParseRule *-- "0..1" SampleNameCleanup
    Axis *-- Duplicates
    Columns *-- "2" ColumnGroup : obs / var
    ColumnGroup *-- "*" ColumnCompute : compute
    Modifications <|-- TokenRegexModifications
    Modifications <|-- AlreadyProformaModifications
    Modifications <|-- SeparateModColumnModifications
    Fragments <|-- PositionalFragments
    Fragments <|-- ColumnLabeledFragments
    TokenRegexModifications *-- "*" ModificationMapEntry : map
    SeparateModColumnModifications *-- "*" ModificationMapEntry : map
```

---

## `readers/` — tabular file reading

**Reads a quant table into a pandas `DataFrame`.** A leaf subpackage.

- `tabular.py` — per-format readers (csv / tsv / parquet).
- `dispatch.read_table` — picks a reader by file extension.

```mermaid
flowchart LR
    path[path] --> rt[dispatch.read_table]
    rt -->|.csv| csv[tabular.read_csv]
    rt -->|.tsv / .txt| tsv[tabular.read_tsv]
    rt -->|.parquet| pq[tabular.read_parquet]
    csv --> df[(DataFrame)]
    tsv --> df
    pq --> df
```

---

## `converters/` — DataFrame + ParseRule → AnnData

**Turns a `DataFrame` + a `ParseRule` into an `AnnData`.**

- `assemble.convert` — orchestrates: optional modification normalization → column
  materialization → long/wide strategy → assemble.
- `long.py` / `wide.py` — the two shape strategies, each returning `ConversionPieces`.
- `assemble.to_anndata` — builds the `AnnData` (rule stored in `uns`).
- `recognize.py` — auto-picks a rule from table headers.
- When a `params_path` is given, parsed `Parameters` are attached to `uns`.

```mermaid
flowchart TD
    df[(DataFrame)] --> conv[assemble.convert]
    rule[ParseRule] --> conv
    headers[table headers] -.-> recog[recognize.recognize]
    recog -. selects .-> rule
    conv --> mods{rule.modifications?}
    mods -->|yes| apply[modifications.apply_modifications]
    mods -->|no| mat[_materialize_columns: select + compute]
    apply --> mat
    mat --> shape{input_shape}
    shape -->|long| long[long.convert_long]
    shape -->|wide| wide[wide.convert_wide]
    long --> pieces[ConversionPieces]
    wide --> pieces
    pieces --> asm[assemble.to_anndata]
    asm --> ad[(AnnData)]
    paramsp[params_path?] -. optional .-> attach[_attach_search_parameters]
    attach -. writes uns .-> ad
```

---

## End-to-end flows (cross-module)

### A. Vendor parameter file → `Parameters`

```mermaid
sequenceDiagram
    autonumber
    participant Caller
    participant Reg as params.registry
    participant V as params vendor
    participant Common as params._common
    participant Model as params.model.Parameters

    Caller->>Reg: parse_params(path, software)
    Reg->>V: get_parser(software) then extract_params(source)
    V->>Common: read_text / read_lines(source)
    Common-->>V: text / lines
    V->>V: vendor parse (regex, yaml, json, toml, xml)
    V->>Model: build Parameters from fields
    Model->>Model: run validators (enzyme map, FDR, tolerance, mods, ranges)
    Model-->>V: Parameters
    V-->>Caller: Parameters
```

### B. Rule TOML + table → `AnnData` (+ optional params)

```mermaid
sequenceDiagram
    autonumber
    participant Caller
    participant Loader as rules.loader
    participant Reader as readers.dispatch
    participant Conv as converters.assemble
    participant Mods as modifications.pipeline
    participant LW as converters.long/wide
    participant P as params

    Caller->>Loader: load_rule(toml_path)
    Loader-->>Caller: ParseRule
    Caller->>Reader: read_table(path)
    Reader-->>Caller: DataFrame
    Caller->>Conv: convert(df, rule, params_path)
    opt rule.modifications set
        Conv->>Mods: apply_modifications(df, rule.modifications)
        Mods-->>Conv: df plus proforma_sequence, stripped_sequence
    end
    Conv->>Conv: _materialize_columns (select plus compute)
    alt input_shape long
        Conv->>LW: convert_long(df, rule)
    else wide
        Conv->>LW: convert_wide(df, rule)
    end
    LW-->>Conv: ConversionPieces
    Conv->>Conv: to_anndata(pieces, rule)
    opt params_path given
        Conv->>P: _attach_search_parameters(adata, params_path, software)
        P-->>Conv: writes uns search_parameters
    end
    Conv-->>Caller: AnnData
```

> Storage keys: the rule is saved in `uns['anndata_proteomics']['rule_json']`; parsed
> parameters in `uns['anndata_proteomics']['search_parameters']` (with the source path in
> `…['search_parameters_path']`).

---

This document is hand-maintained; when adding a vendor, model field, or rule section, update
the relevant module diagram above. Sources are under `src/anndata_proteomics/`.
