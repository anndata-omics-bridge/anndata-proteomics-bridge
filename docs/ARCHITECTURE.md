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
                                                              │
                                                              ▼
                                                  validate.main → CLI exit code
```

Future flow (not implemented): `ParseRule + vendor data file → readers → converters → AnnData`. See [RESTART_PLAN.md](RESTART_PLAN.md) steps 5–10.

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

**Tests** — covered indirectly by `test_rule_models.test_json_schema_export_has_expected_top_level_properties`.

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

## Console scripts

Wired in [pyproject.toml](../pyproject.toml) under `[project.scripts]`:

| Command | Module | Purpose |
|---|---|---|
| `validate-rules` | `rules.validate:main` | Sanity-check every packaged TOML; exit code suitable for CI |
| `export-rule-schema` | `rules._export_schema:main` | Regenerate `parse_rule.schema.json` from the pydantic models |

## Not yet implemented

- `readers/dispatch.py`, `readers/tabular.py` — vendor file (`.tsv` / `.csv` / `.parquet`) → `pandas.DataFrame`.
- `converters/recognize.py` — pick the right `ParseRule` from a vendor file's header.
- `converters/long.py`, `converters/wide.py` — apply a validated rule to a DataFrame.
- `converters/factors.py` — encode string-valued layers as integer factors per the TOML.
- `converters/assemble.py` — assemble `obs`, `var`, `X`, `layers`, `uns` into an `AnnData`.
- `cli.py` — user-facing CLI for validation / listing / conversion.

See [RESTART_PLAN.md](RESTART_PLAN.md) §"First Implementation Order" for the agreed sequence.

## Adding things

- **New vendor TOML** — drop a file at `parsing_rules/<software>/parse_<software>_<level>_<file_version>.toml`. `validate-rules` and the test suite pick it up automatically; no registry edits needed. The level in the filename must match the TOML's `quantification_level`.
- **New schema field** — edit [rules/schema.py](../src/anndata_proteomics/rules/schema.py), update the fixtures in [tests/test_rule_models.py](../tests/test_rule_models.py) and [docs/toml_schema.md](toml_schema.md), then run `export-rule-schema` to regenerate the JSON Schema.
- **New `quantification_level` value** (e.g. `protein`) — already in the `QuantificationLevel` literal; just author the TOMLs and matching tests.
