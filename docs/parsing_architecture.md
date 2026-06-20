# Parsing architecture (UML)

Diagrams of APB's parsing subsystem, covering two related but distinct flows:

1. **Parameter parsing** — a vendor parameter file → a typed `Parameters` record.
2. **Table conversion** — a parsing-rule TOML + a vendor quant table → an `AnnData`
   (optionally attaching parsed search parameters into `uns`).

Diagrams are [Mermaid](https://mermaid.js.org/) and render in GitHub and most IDEs. Sources
referenced are under `src/anndata_proteomics/`.

---

## 1. Class diagram — data models

Three model families: **params** (`params/model.py`), **modifications**
(`modifications/model.py`, `apply_rules.py`, `unimod_registry.py`), and **rules**
(`rules/schema.py`), plus the `ConversionPieces` container.

```mermaid
classDiagram
    direction LR

    %% ---- params/model.py ----
    class _Strict {
        <<pydantic BaseModel, extra=forbid>>
    }
    class Probability {
        +float value
    }
    class MassTolerance {
        +float|None value
        +str|None unit
        +str mode
        +str|None label
        +parse(value) MassTolerance$
    }
    class UnparsedParameter {
        +str name
        +ScalarValue value
        +str|None source
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
    note for Parameters "plus ~20 more scalar fields: charges, m/z, peptide lengths, FDR, quant method"
    _Strict <|-- Probability
    _Strict <|-- MassTolerance
    _Strict <|-- UnparsedParameter
    _Strict <|-- Parameters
    Parameters *-- "0..3" Probability : ident_fdr
    Parameters *-- "0..2" MassTolerance : precursor/fragment tol
    Parameters *-- "*" SearchedModification : fixed_mods / variable_mods
    Parameters *-- "*" UnparsedParameter : unparsed_parameters

    %% ---- modifications/model.py ----
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
        +str|None position
        +float|None mass_delta
        +str|None source
    }
    class ModificationOccurrence {
        +str name
        +str|None accession
        +str|None target_residue
        +int|None sequence_index
        +str|None position
        +float|None mass_delta
        +str|None source_token
    }
    class ModifiedSequence {
        +str stripped_sequence
        +str proforma_sequence
        +str|None source_sequence
        +list~str~ unknown_tokens
    }
    SearchedModification --> ModType
    ModifiedSequence *-- "*" ModificationOccurrence : occurrences

    %% ---- modifications/apply_rules.py + unimod_registry.py ----
    class MapEntry {
        <<frozen>>
        +str token
        +str name
        +str|None accession
        +str|None target
        +str|None position
        +float|None mass_delta
    }
    class ModificationRule {
        <<frozen>>
        +str source_column
        +str token_pattern
        +str token_position
        +bool case_sensitive
        +str unknown_policy
        +str output_column
    }
    class UnimodEntry {
        +str accession
        +str name
        +str target
        +str position
        +float mass_delta
    }
    ModificationRule *-- "*" MapEntry : entries

    %% ---- converters/_pieces.py ----
    class ConversionPieces {
        +ndarray X
        +DataFrame obs
        +DataFrame var
        +dict layers
        +dict uns
    }
```

The rule-schema models (`rules/schema.py`) compose as follows:

```mermaid
classDiagram
    direction LR
    class ParseRule {
        +str schema_version
        +str software_name
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
        +str encoding_mode
        +str|None source_column
        +str|None column_pattern
    }
    class Modifications {
        +str source_column
        +str parser
        +str|None token_pattern
        +str unknown_policy
        +str output_column
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
    ParseRule *-- "0..1" SampleNameCleanup
    Axis *-- Duplicates
    Columns *-- "2" ColumnGroup : obs / var
    ColumnGroup *-- "*" ColumnCompute : compute
    Modifications *-- "*" ModificationMapEntry : map
```

> `Parameters.fixed_mods/variable_mods` hold `SearchedModification` (defined in
> `modifications/model.py` but re-exported through `params` for SDRF use). The runtime
> `ModificationRule` (in `apply_rules.py`) is the resolved form of the TOML `Modifications`
> section — `modifications/pipeline._to_runtime_rule` fills each `MapEntry` from the
> `UnimodEntry` registry.

---

## 2. Flow — vendor parameter file → `Parameters`

```mermaid
sequenceDiagram
    autonumber
    participant Caller
    participant Reg as params.registry
    participant V as params.<vendor>
    participant Common as params._common
    participant Model as params.model.Parameters

    Caller->>Reg: parse_params(path, software)
    Reg->>Reg: get_parser(software) -> extract_params
    Reg->>V: extract_params(source)
    V->>Common: read_text / read_lines(source)
    Common-->>V: text / lines
    V->>V: vendor parse (regex / yaml / json / toml / xml)
    V->>Model: Parameters(**fields)
    Model->>Model: field validators (enzyme map, FDR>=1,<br/>MassTolerance.parse, mod ProForma, ranges)
    Model-->>V: Parameters
    V-->>Caller: Parameters
```

## 2b. Flow — rule TOML + table → `AnnData` (+ optional params)

```mermaid
sequenceDiagram
    autonumber
    participant Caller
    participant Loader as rules.loader
    participant Reader as readers.dispatch
    participant Conv as converters.assemble.convert
    participant Mods as modifications.pipeline
    participant LW as converters.long/wide
    participant Asm as converters.assemble.to_anndata
    participant Params as params.registry + anndata_io

    Caller->>Loader: load_rule(toml) -> ParseRule
    Caller->>Reader: read_table(path) -> DataFrame
    Caller->>Conv: convert(df, rule, params_path=?)
    opt rule.modifications is not None
        Conv->>Mods: apply_modifications(df, rule.modifications)
        Mods->>Mods: _to_runtime_rule (unimod resolve) + apply_rule per row
        Mods-->>Conv: df + proforma_sequence / stripped_sequence
    end
    Conv->>Conv: _materialize_columns (select + compute)
    alt input_shape == long
        Conv->>LW: convert_long(df, rule) -> ConversionPieces
    else wide
        Conv->>LW: convert_wide(df, rule) -> ConversionPieces
    end
    Conv->>Asm: to_anndata(pieces, rule)
    Asm-->>Conv: AnnData (uns['anndata_proteomics']['rule_json'])
    opt params_path provided
        Conv->>Params: _attach_search_parameters(adata, params_path, software)
        Params->>Params: parse_params -> Parameters; write_search_parameters
        Note over Params: uns['anndata_proteomics']['search_parameters'(_path)]
    end
    Conv-->>Caller: AnnData
```

> Rule auto-detection: `converters.recognize.recognize(headers) -> ParseRule | None`
> picks the unique packaged rule whose `matches(headers, rule)` holds, when the caller
> doesn't supply a rule explicitly.

---

## 3. Package / component overview

```mermaid
flowchart TD
    subgraph readers
        dispatch[dispatch.read_table]
        tabular[tabular.read_csv/tsv/parquet]
        dispatch --> tabular
    end

    subgraph rules
        schema[schema.ParseRule + Axis/Columns/Layer/Modifications]
        loader[loader.load_rule]
        rreg[registry.find_rule]
        loader --> schema
        loader --> rreg
    end

    subgraph params
        preg[registry.get_parser / parse_params]
        vendors[vendor extract_params x10]
        pmodel[model.Parameters / MassTolerance / Probability]
        pcommon[_common.read_text / read_lines]
        pio[anndata_io.write/read_search_parameters]
        preg --> vendors
        vendors --> pmodel
        vendors --> pcommon
        pio --> pmodel
    end

    subgraph modifications
        mpipe[pipeline.apply_modifications]
        mapply[apply_rules.apply_rule]
        mpro[proforma.render_proforma]
        munic[unimod_registry.resolve]
        mmodel[model.*]
        msdrf[sdrf.to/from_sdrf_value]
        mpipe --> mapply
        mpipe --> munic
        mapply --> mpro
        mapply --> mmodel
    end

    subgraph converters
        conv[assemble.convert]
        recog[recognize.recognize]
        lw[long.convert_long / wide.convert_wide]
        asm[assemble.to_anndata]
        pieces[_pieces.ConversionPieces]
        conv --> lw
        conv --> asm
        lw --> pieces
    end

    conv --> rules
    conv --> readers
    conv --> mpipe
    conv --> preg
    conv --> pio
    recog --> rules
    msdrf -.-> mmodel
```

---

## Notes

- `params/` is standalone (no `proteobench` imports); `extract_params(source) -> Parameters`
  is the uniform vendor entry point, dispatched via `params/registry.py`.
- Source acquisition is centralized in `params/_common.py` (`read_text` / `read_lines`);
  vendor modules only own the format-specific parse step.
- Parsed parameters live in `adata.uns['anndata_proteomics']['search_parameters']` (JSON), with
  the originating path in `…['search_parameters_path']`; the parsing rule is stored in
  `…['rule_json']`.
- This document is hand-maintained; when adding a vendor, model field, or rule section, update
  the relevant diagram above.
