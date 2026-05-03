# Architecture (current state)

**Status as of 2026-05-02** (HEAD `09ee417`). This document describes what is *implemented today*. For the broader design and remaining steps, see [RESTART_PLAN.md](RESTART_PLAN.md).

## Data flow

```
TOML rule file
      │
      ▼                              schema.ParseRule  ◄── _export_schema ──► parse_rule.schema.json
loader.load_rule  ───────────────►  (validated rule)                          (IDE / CI consumers)
      ▲
      │ via
      │
registry.find_rule(software, level, version)
      ▲
      │ resolves against
      │
parsing_rules/<vendor>/parse_<software>_<level>_<version>.toml
                                                              ▲
                                                              │ enumerated by
                                                              │
                                              registry.iter_packaged_rules
                                                              │ used by
                                                              ▼
                                                  validate.validate_all_packaged

Vendor data file ──► readers.read_table ──► pandas.DataFrame
                          │ dispatches on extension to
                          ▼
                   readers.tabular.read_csv / read_tsv / read_parquet

DataFrame ──► converters.recognize ──► ParseRule (auto-pick from packaged set)
DataFrame + ParseRule ──► converters.convert ──► AnnData
                              │ dispatches on rule.input_shape to
                              ▼
                       converters.long.convert_long  /  converters.wide.convert_wide
                              │ then
                              ▼
                       converters.assemble.to_anndata  ◄── factors.encode_factor (per layer)
```

## Modules

### `rules/schema.py`

**Purpose** — pydantic models for the parsing-rule TOML format. Single source of truth for what a valid rule looks like; the JSON Schema is derived from these models.

**Public API**

- `ParseRule` — top-level container. Validates structure plus cross-field rules: long/wide layer consistency, `axis.x_layer` exists in `layers[*].name`, factor-encoded layers require non-empty `categories`, `sample_name_cleanup` only valid for wide rules, no unknown top-level keys (`extra='forbid'`).
- `Axis`, `Columns`, `Layer`, `Duplicates`, `SampleNameCleanup` — per-section sub-models, all strict.
- `InputShape`, `QuantificationLevel`, `EncodingMode`, `DuplicateMode` — `Literal` type aliases used as field types.

**Depends on** — `pydantic` only.

**Tests** — [tests/test_rule_models.py](../tests/test_rule_models.py)

### `rules/loader.py`

**Purpose** — read a TOML file and return a validated `ParseRule`.

**Public API**

- `load_rule(path) -> ParseRule` — `tomllib.loads` + `ParseRule.model_validate`. Raises `FileNotFoundError` on missing path; on pydantic `ValidationError`, attaches the file path as an exception note before re-raising.
- `load_packaged_rule(software, quantification_level, file_version="1") -> ParseRule` — sugar over `find_rule + load_rule`.

**Depends on** — `rules.schema`, `rules.registry`, stdlib `tomllib`.

**Tests** — [tests/test_rule_loader.py](../tests/test_rule_loader.py)

### `rules/registry.py`

**Purpose** — locate packaged TOMLs by `(software, level, version)` and enumerate them.

**Public API**

- `packaged_rules_root() -> Path` — `parsing_rules/` inside the installed package, via `importlib.resources` (works for editable installs and wheels).
- `iter_packaged_rules() -> Iterator[Path]` — sorted glob of all `parsing_rules/<vendor>/parse_*.toml`.
- `find_rule(software, quantification_level, file_version="1") -> Path` — resolve to a specific TOML.
- `RuleNotFound(LookupError)` — raised by `find_rule`; message lists what *is* available in the vendor folder.

**Depends on** — stdlib only (`importlib.resources`, `pathlib`).

**Tests** — [tests/test_rule_registry.py](../tests/test_rule_registry.py)

### `rules/validate.py`

**Purpose** — produce CLI-friendly batch validation results.

**Public API**

- `ValidationResult` — frozen dataclass `{path, ok, error, rule}`. `error` is a one-line summary string; `rule` is populated when `ok=True`.
- `validate_file(path) -> ValidationResult` — never raises; failures come back as `ok=False`.
- `validate_all_packaged() -> list[ValidationResult]` — walks every packaged rule.
- `main(argv=None) -> int` — `validate-rules` console entrypoint. Prints `PASS path` / `FAIL path: msg` per rule plus a summary line, returns 0 if all pass else 1.

**Depends on** — `rules.loader`, `rules.registry`, `rules.schema`.

**Tests** — [tests/test_rule_validate.py](../tests/test_rule_validate.py)

### `rules/_export_schema.py`

**Purpose** — emit `parse_rule.schema.json` from `ParseRule.model_json_schema()` as a side-output for IDE tooling (Even Better TOML / taplo) and CI sanity checks.

**Public API**

- `main()` — `export-rule-schema` console entrypoint; writes to `parsing_rules/_schema/parse_rule.schema.json`.

**Depends on** — `rules.schema`, stdlib `json`.

**Tests** — covered indirectly by `test_rule_models.test_json_schema_export_has_expected_top_level_properties`. [tests/test_json_schema_validation.py](../tests/test_json_schema_validation.py) additionally round-trips every packaged TOML through `jsonschema.validate(...)` against the generated schema (structural-parity smoke test). Pydantic remains the only source of truth for cross-field rules — JSON Schema is strictly weaker, by design.

### `parsing_rules/`

**Purpose** — packaged vendor TOML rules, one folder per vendor, plus the generated JSON Schema under `_schema/`.

**Current contents**

```
parsing_rules/
├── _schema/parse_rule.schema.json              generated, IDE-consumed
├── diann/parse_diann_ion_1.toml                long, ion
├── spectronaut/parse_spectronaut_ion_1.toml    long, ion
├── maxquant/parse_maxquant_ion_1.toml          long, ion
├── fragpipe/parse_fragpipe_ion_1.toml          wide, ion
├── peaks/parse_peaks_ion_1.toml                wide, ion
└── wombat/parse_wombat_peptidoform_1.toml      wide, peptidoform
```

**Filename convention** — `parse_<software>_<level>_<file_version>.toml`. The level token must match the TOML's `quantification_level` field; [tests/test_packaged_rules.py](../tests/test_packaged_rules.py) enforces this.

### `readers/tabular.py`

**Purpose** — generic per-extension readers; no vendor semantics, no rule application.

**Public API**

- `read_csv(path)` — comma-delimited, UTF-8 with BOM tolerance.
- `read_tsv(path)` — tab-delimited, UTF-8 with BOM tolerance.
- `read_parquet(path)` — via pyarrow.

Each is a thin wrapper around `pd.read_csv` / `pd.read_parquet` so the test suite has stable entry points and any future overrides land in one place.

**Depends on** — `pandas`, `pyarrow`.

**Tests** — [tests/test_readers_tabular.py](../tests/test_readers_tabular.py)

### `readers/dispatch.py`

**Purpose** — pick the right reader from a file's extension. Vendor differences across the 6 packaged TOMLs collapse to four extensions; this is the only piece of code that knows that.

**Public API**

- `read_table(path) -> pd.DataFrame` — dispatch by `path.suffix.lower()`.
- `EXTENSION_TO_READER` — registry mapping. `.txt` is treated as TSV (MaxQuant convention).
- `UnknownFormat(ValueError)` — raised when the extension is not registered.

**Depends on** — `readers.tabular`.

**Tests** — [tests/test_readers_dispatch.py](../tests/test_readers_dispatch.py); end-to-end coverage in [tests/test_readers_integration.py](../tests/test_readers_integration.py), which parametrizes over every packaged TOML and reads the matching test_data_download file (skips if the gitignored cache is absent).

## `scripts/cli.py` — `anndata-proteomics` umbrella CLI

**Purpose** — single user-facing CLI with subcommands. Built on `cyclopts`.

**Subcommands**

| Subcommand | Behavior | Exit code |
|---|---|---|
| `validate [path ...]` | Validate one or more TOML rules. With no path, walks all packaged rules. | 0 if all pass, 1 otherwise |
| `list` | List packaged rules: software, level, file_version, path. | 0 |
| `export-schema` | Regenerate `parse_rule.schema.json`. | 0 |
| `convert <data> <rule.toml>` | **STUB** — not yet implemented (RESTART_PLAN steps 5–10). | 2 (distinct from validation failure 1) |

**Implementation principle** — `validate` shares the PASS/FAIL formatter (`rules.validate._print_and_exit_code`) with the older `validate-rules` console script, so output is identical.

**Tests** — [tests/test_cli.py](../tests/test_cli.py)

## Console scripts

Wired in [pyproject.toml](../pyproject.toml) under `[project.scripts]`:

| Command | Module | Purpose |
|---|---|---|
| `anndata-proteomics` | `scripts.cli:main` | Umbrella CLI with `validate / list / export-schema / convert` subcommands |

### `converters/recognize.py`

**Purpose** — match a DataFrame's column headers to one of the packaged `ParseRule`s.

**Public API**

- `matches(headers, rule) -> bool` — does this rule plausibly fit?
- `recognize(headers) -> ParseRule | None` — unique match, or None on zero / multiple.

**Tests** — [tests/test_recognize.py](../tests/test_recognize.py)

### `converters/long.py`

**Purpose** — apply a long-format `ParseRule` to a DataFrame: pivot to per-layer (obs × var) matrices.

**Public API**

- `convert_long(df, rule) -> ConversionPieces` — full pipeline (build obs/var, pivot every layer, factor-encode where needed). Honors `duplicates.mode`. Coerces non-factor layers via `pd.to_numeric(errors='coerce')` so vendor sentinels like `"-"` become NaN rather than blowing up the pivot.

**Tests** — [tests/test_converters_long.py](../tests/test_converters_long.py)

### `converters/wide.py`

**Purpose** — apply a wide-format `ParseRule` to a DataFrame: extract sample tokens from column headers via each layer's `column_pattern`, build per-layer matrices.

**Public API**

- `convert_wide(df, rule) -> ConversionPieces` — extracts samples (union across layers, insertion-order), builds var from `[columns.var]`, gathers each layer's matching columns into an `(n_obs × n_var)` matrix. Applies `sample_name_cleanup.pattern` if present.

**Tests** — [tests/test_converters_wide.py](../tests/test_converters_wide.py)

### `converters/factors.py`

**Purpose** — encode string-valued layer data to integer codes per the TOML `categories` map.

**Public API**

- `encode_factor(series, categories, default=-1) -> Series[int64]` — unknowns and NaN → `default`.

**Tests** — [tests/test_converters_factors.py](../tests/test_converters_factors.py)

### `converters/assemble.py`

**Purpose** — assemble `ConversionPieces` into an `AnnData`, plus `convert(df, rule)` umbrella.

**Public API**

- `to_anndata(pieces, rule) -> ad.AnnData` — wraps the pieces; writes `uns['anndata_proteomics']` with `rule`, `schema_version`, `software_name`, `input_shape`, `quantification_level`.
- `convert(df, rule) -> ad.AnnData` — dispatches to `convert_long` or `convert_wide` based on `rule.input_shape`, then assembles.

**Tests** — [tests/test_converters_assemble.py](../tests/test_converters_assemble.py); end-to-end coverage for all 6 packaged vendors in [tests/test_converters_e2e.py](../tests/test_converters_e2e.py).

## Not yet implemented

- Hook `convert(df, rule)` into the `anndata-proteomics convert` CLI subcommand (still a stub).
- `converters/long.py`, `converters/wide.py` — apply a validated rule to a DataFrame.
- `converters/factors.py` — encode string-valued layers as integer factors per the TOML.
- `converters/assemble.py` — assemble `obs`, `var`, `X`, `layers`, `uns` into an `AnnData`.
- `cli.py` — user-facing CLI for validation / listing / conversion.

See [RESTART_PLAN.md](RESTART_PLAN.md) §"First Implementation Order" for the agreed sequence.

## Adding things

- **New vendor TOML** — drop a file at `parsing_rules/<software>/parse_<software>_<level>_<file_version>.toml`. `validate-rules` and the test suite pick it up automatically; no registry edits needed. The level in the filename must match the TOML's `quantification_level`.
- **New schema field** — edit [rules/schema.py](../src/anndata_proteomics/rules/schema.py), update the fixtures in [tests/test_rule_models.py](../tests/test_rule_models.py) and [docs/toml_schema.md](toml_schema.md), then run `export-rule-schema` to regenerate the JSON Schema.
- **New `quantification_level` value** (e.g. `protein`) — already in the `QuantificationLevel` literal; just author the TOMLs and matching tests.
