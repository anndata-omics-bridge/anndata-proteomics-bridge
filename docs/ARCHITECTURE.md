# Architecture

APB is a Python library plus the `apb` CLI. It converts proteomics vendor tables
to AnnData or MuData using packaged JSON parsing rules, and can attach search
parameters and external annotations.

## Data Flow

**Figure:** the main conversion path. Rule JSON documents are validated into `ParseRule`
objects; vendor data becomes a `DataFrame`; optional parameter files, sample
annotation, and FASTA validation enrich the final AnnData/MuData object.
ProteoBench scoring consumes the sample design produced by annotation.

```mermaid
flowchart TD
    config[/software-version rules.json/] --> resolve[registry.resolve_rule_locator / find_rule]
    resolve --> load[loader.load_rule_document + effective_rule]
    load --> source([schema.ParseRuleDocument validated])
    load --> rule([schema.ParseRule validated])
    source -. export-schema .-> sourceschema[[parse_rule_document.schema.json]]
    rule -. export-schema .-> schemajson[[parse_rule.schema.json]]
    rule -. apb validate .-> val[validate: PASS / FAIL]

    data[/vendor data file/] --> readtable[readers.read_table]
    readtable --> df[(pandas.DataFrame)]
    df -. recognize auto-pick .-> rule

    df --> convert[converters.assemble.convert]
    rule --> convert
    convert --> mods{rule.modifications?}
    mods -->|yes| applymods[modifications.apply_modifications]
    mods -->|no| shape{input_shape}
    applymods --> shape
    shape -->|long| long[long.convert_long]
    shape -->|wide| wide[wide.convert_wide]
    long --> pieces[ConversionPieces]
    wide --> pieces
    pieces --> assemble[assemble.to_anndata + factors.encode_factor]
    assemble --> adata[(AnnData / MuData)]

    paramfile[/vendor parameter file/] --> parseparams[params.registry.parse_params]
    parseparams --> params([params.model.Parameters])
    params -. write_search_parameters .-> adata

    adata --> annobs[annotation.annotate_obs: obs axis]
    fastafile[/FASTA files/] --> annvar[var_fasta.annotate_var_from_fasta: protein varm]
    adata --> annvar
    fastafile --> valfasta[validate_fasta: peptide-derived varm]
    adata --> valfasta
    valfasta -. MuData .-> mulink[varp feature_mapping: peptide feature to protein feature]

    moduletoml[/ProteoBench module TOML/] --> pbscore[proteobench.score_quantification]
    annobs --> pbscore
    pbscore --> pbvarm[varm proteobench: feature statistics]
    pbscore --> pbuns[uns anndata_proteomics/proteobench: roles, mapping provenance, scores]

    classDef io fill:#eef2ff,stroke:#9aa7d8;
    class config,data,paramfile,fastafile,moduletoml io;
```

More detailed diagrams are in [parsing_architecture.md](parsing_architecture.md).
Search-parameter parser details are in [parameter_parsers.md](parameter_parsers.md).

## Package Map

| Area | Current role |
|---|---|
| `rules/` | Pydantic source/effective rule composition, document merge loader, registry, validator, and JSON Schema export. |
| `parsing_rules/` | Packaged JSON rules under `src/anndata_proteomics/parsing_rules/`. |
| `readers/` | File-extension dispatch to CSV, TSV/TXT, and Parquet readers. |
| `converters/` | Rule-driven long/wide conversion into AnnData; multi-level CLI conversion can assemble MuData. |
| `modifications/` | Parsing-rule modification schemas, vendor sequence normalization to ProForma, and searched-modification models. |
| `params/` | Typed search-parameter model, parser registry, AnnData storage helpers, and vendor parsers under `params/parsers/`. |
| `annotation/` | `obs` annotation, FASTA-derived protein `varm['fasta']`, peptide `varm['fasta_validation']`, and MuLink-compatible `varp['feature_mapping']`. |
| `fasta/` | FASTA parsing, typed decoy/contaminant configuration, protein metadata, and enzyme-aware theoretical peptide counts. |
| `proteobench/` | Typed module TOMLs, canonical role resolution, matrix-native HYE intermediates, compatible metrics, and storage orchestration. |
| `prozor` dependency | Backend-neutral Aho--Corasick matching and reusable protein-inference primitives; APB owns FASTA parsing and AnnData/MuData storage. |
| `scripts/` | The installed `apb` CLI. |

The core parsing packages follow one tested dependency order:
`rules → params → modifications`. `rules` also imports `modifications` directly when composing
`ParseRule`; `params` reuses canonical modification identities; `modifications` never imports
either higher layer. Empty package initializers keep these edges visible to static package
analysis.

## Packaged Rules

Packaged rule files live inside the Python package:

```text
src/anndata_proteomics/parsing_rules/
  _schema/parse_rule.schema.json
  _schema/parse_rule_document.schema.json
  alphapept/rules.json       # ion
  diann/v1/rules.json        # ion, fragment, protein
  diann/v2/rules.json        # ion, protein
  fragpipe/rules.json        # ion
  maxquant/rules.json        # ion
  peaks/rules.json           # ion
  sage/rules.json            # ion, peptidoform (parameter-gated)
  spectronaut/rules.json     # ion, fragment, protein
  wombat/rules.json          # ion, peptidoform
```

`rules.registry.resolve_rule_locator()` selects a software-version document by its
required regex and addresses a level inside it. Vendor-root and `v*/` documents use
the same source shape.

## CLI Surface

`pyproject.toml` installs one console script:

```bash
apb
```

The CLI subcommands are:

| Subcommand | Purpose |
|---|---|
| `apb validate [path ...]` | Validate rule JSON; with no paths, validate all packaged rules. |
| `apb list` | List packaged rules and their metadata. |
| `apb export-schema` | Regenerate the source-document and effective-rule schemas. |
| `apb convert <data> [level] --params <param-file>` | Convert vendor data to `.h5mu` or a selected `.h5ad` level. |
| `apb annotate <data> <annotations.toml/csv/tsv>` | Join external sample metadata onto `obs`. |
| `apb fasta <data> <proteome.fasta>` | Annotate proteins and, by default, validate every peptide-derived modality against FASTA; accession-dependent work uses only `column_roles.fasta_accessions` or an explicit override. |
| `apb proteobench <data> <module.toml>` | Score every annotated quantification level a container holds, each into its own `uns`/`varm`; requires `sample_name`, `condition`, and `column_roles.protein_assignment`, while FASTA remains optional. |

## Search Parameters

Generic code stays at `params/`:

```text
params/
  anndata_io.py
  model.py
  registry.py
  parsers/
```

Vendor-specific parser implementations live only under `params/parsers/`.
Their native `extract_params(...)` signatures may express vendor-specific
inputs or options. `params.registry` dispatches by software name through a
uniform callable that accepts either one source or an explicit source tuple.
MetaMorpheus, for example, receives its TOML and version-text files as a
two-source tuple; the registry never guesses a related filename.

## Current Limits

- Conversion coverage is limited to the packaged vendor/level rules above.
- MuLink edges can only target protein features already present in a MuData;
  exhaustive FASTA matches remain available in each peptide modality's validation table.
- `duplicates.mode = "keep_all_as_raw_table"` is reserved but not implemented.
- ProteoBench coverage currently implements the quantitative HYE metrics;
  PYE/plasma, de-novo, and entrapment variants are not implemented.
- Scoring is level-agnostic: the module TOML supplies the sample design and the
  per-species expected ratios, and the feature axis is `var_names` (the rule's
  joined `axis.var_keys`), so any converted level can be scored. Only ion and
  peptidoform have published ProteoBench modules, so fragment and protein scores
  are comparable across APB runs but not to the ProteoBench leaderboard.

## Adding Things

- New vendor table rule: add JSON under
  `src/anndata_proteomics/parsing_rules/<software>/`, then run `apb validate`.
- New vendor parameter parser: add `params/parsers/<vendor>.py`, register it in
  `params/registry.py`, and add parser fixtures/tests.
- New schema field: edit `rules/schema.py`, update tests and
  [json_schema.md](json_schema.md), then run `apb export-schema`.
