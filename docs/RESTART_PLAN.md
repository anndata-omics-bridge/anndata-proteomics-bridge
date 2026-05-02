# Restart Plan For `anndata_proteomics_bridge`



## Goal

Restart `anndata_proteomics_bridge` as a focused package that:

- ships parsing-rule TOMLs with the code
- validates those TOMLs with `pydantic`
- reads vendor quant files
- converts them into AnnData using the TOML rules

Out of scope for this restart:

- second-stage `obs` annotation
- adding conditions, factors, batches, or other later metadata enrichment
- project-specific downstream benchmarking logic

This project should stop at:

- vendor file + parsing TOML -> AnnData

## Main Decisions

- Keep the Python package name as `anndata_proteomics`
- Keep the repository name as `anndata_proteomics_bridge`
- Ship required TOMLs inside `src/anndata_proteomics/parsing_rules/`
- Group TOMLs by vendor, not by wide/long
- Keep wide vs long as a TOML property: `input_shape`
- Encode quantification level (`ion` / `peptidoform` / `peptide` / `protein`) both as a TOML
  field (`quantification_level`) **and** as a token in the filename
  (`parse_<software>_<level>_<file_version>.toml`); a test enforces the two stay in sync
- Do not implement second-stage `obs` annotation here

## Status (2026-05-01)

- ✅ Pre-restart `src/` cleared (commit `a4dec1c`)
- ✅ `rules/schema.py` with pydantic `ParseRule` + JSON Schema export to
  `parsing_rules/_schema/parse_rule.schema.json` (commit `acd64d8`)
- ✅ Six packaged TOMLs committed in `parsing_rules/` (commit `29c773c`):
  `diann`, `spectronaut`, `maxquant` (long, ion); `fragpipe`, `peaks` (wide, ion);
  `wombat` (wide, peptidoform). `proteome_discoverer` dropped — no submissions in
  the ProteoBench catalog (`raw_file_db_full.csv`).
- ✅ `rules/loader.py` + `rules/registry.py` + `rules/validate.py` (step 3 done).
- ✅ JSON Schema round-trip tests in `tests/test_json_schema_validation.py` (commit `4f4f252`).
- ✅ `scripts/cli.py` — `anndata-proteomics` umbrella CLI with `validate / list /
  export-schema / convert (stub)` subcommands (cli.py skeleton in place; `convert`
  remains a stub until readers/converters land).
- ⏭️ Next: `readers/dispatch.py` + `readers/tabular.py` (step 5 of "First Implementation Order")

## Proposed Package Layout

```text
anndata_proteomics_bridge/
├── pyproject.toml
├── README.md
├── docs/
│   └── RESTART_PLAN.md
├── test_data_download/
├── tests/
└── src/
    └── anndata_proteomics/
        ├── __init__.py
        ├── cli.py
        ├── parsing_rules/
        │   ├── _schema/             # generated parse_rule.schema.json (committed)
        │   ├── diann/
        │   ├── fragpipe/
        │   ├── maxquant/
        │   ├── peaks/
        │   ├── spectronaut/
        │   └── wombat/
        ├── rules/
        │   ├── __init__.py
        │   ├── loader.py
        │   ├── registry.py
        │   ├── schema.py
        │   └── validate.py
        ├── readers/
        │   ├── __init__.py
        │   ├── dispatch.py
        │   └── tabular.py
        └── converters/
            ├── __init__.py
            ├── assemble.py
            ├── factors.py
            ├── recognize.py
            ├── long.py
            └── wide.py
```

## What Goes In Each File

### `src/anndata_proteomics/cli.py`

Top-level command-line entrypoint. It should expose rule validation, rule listing, and file
conversion commands, and then dispatch into the package code.

## `src/anndata_proteomics/parsing_rules/`

Packaged TOML rules required at runtime. The rules should be grouped by vendor, with one folder per
vendor.

```text
parsing_rules/
├── _schema/             # generated parse_rule.schema.json (commit acd64d8)
├── diann/               # parse_diann_ion_1.toml
├── fragpipe/            # parse_fragpipe_ion_1.toml
├── maxquant/            # parse_maxquant_ion_1.toml
├── peaks/               # parse_peaks_ion_1.toml
├── spectronaut/         # parse_spectronaut_ion_1.toml
└── wombat/              # parse_wombat_peptidoform_1.toml
```

Filename convention: `parse_<software>_<quantification_level>_<file_version>.toml`. The level
must match the in-TOML `quantification_level` field — `tests/test_packaged_rules.py` enforces it.

## `src/anndata_proteomics/rules/`

Rule loading, validation, discovery, and schema definitions should live together. Keep this as one
coherent subsystem instead of splitting schema models into a separate top-level package.

### `src/anndata_proteomics/rules/schema.py`

Formal `pydantic` schema for the TOML files. It should define the common rule structure, the long
and wide variants, and the conditional validation rules. A side script (`_export_schema.py`) emits
`parsing_rules/_schema/parse_rule.schema.json` so editors that combine TOML with JSON Schema get
live validation in `.toml` files. Pydantic is the single source of truth; JSON Schema is derived.

### `src/anndata_proteomics/rules/loader.py`

Load TOML files from package resources or explicit paths. It should parse raw TOML into Python data
and return validated rule objects from `rules/schema.py`.

### `src/anndata_proteomics/rules/registry.py`

Discover and resolve the packaged rules. It should support listing available rules and locating the
right rule file by vendor or filename.


### `src/anndata_proteomics/rules/validate.py`

Validation entrypoints for one rule or all packaged rules. It should provide clear validation
errors for CLI use and tests.

## `src/anndata_proteomics/readers/`

Generic file-reading code should stay separate from rule logic and conversion logic. This part
should only read the vendor files into tables without interpreting vendor semantics.

### `src/anndata_proteomics/readers/dispatch.py`

Choose the right low-level reader based on file type or explicit options. It should dispatch by
extension and support delimiter overrides when needed.

### `src/anndata_proteomics/readers/tabular.py`

Generic tabular file reading into pandas DataFrames. It should handle csv, tsv, txt, and parquet
without adding vendor-specific parsing semantics.

## `src/anndata_proteomics/converters/`

Conversion code should take a loaded table plus a validated rule and produce the pieces needed to
build AnnData. This is also the right place for rule-based format recognition from vendor headers.

### `src/anndata_proteomics/converters/recognize.py`

Recognize the vendor file type from the input header and match it to the correct parsing rule when
possible. The goal is to avoid requiring the user to manually specify the vendor type if the file
can be recognized reliably.

### `src/anndata_proteomics/converters/long.py`

Apply validated long-format rules. It should build `obs`, `var`, and the layer-ready data from the
configured source columns and enforce the duplicate policy.

### `src/anndata_proteomics/converters/wide.py`

Apply validated wide-format rules. It should find layer columns by pattern, extract sample names,
build `obs` and `var`, and reshape the file into layer-ready matrices while keeping vendor-derived
obs names at this stage.

### `src/anndata_proteomics/converters/factors.py`

Factor-encode string-valued matrix data. It should convert configured string layers to integer
codes and keep the category mapping consistent with the TOML rule.

### `src/anndata_proteomics/converters/assemble.py`

Construct the final AnnData object from prepared pieces. It should assemble `obs`, `var`, `X`,
`layers`, and `uns` without adding vendor-specific logic.

## `test_data_download/`

Input examples for rule design and conversion tests. Start with a curated subset needed for tests,
not the full uncontrolled dataset.

## `tests/`

Tests for the rule schema and conversion behavior. They should cover rule validation, long and wide
conversion, factor encoding, and packaged-rule discovery.

Recommended test split:

```text
tests/
├── test_rule_models.py        # pydantic schema (in place; commit acd64d8)
├── test_packaged_rules.py     # filename ↔ field, long/wide coverage (commit 29c773c)
├── test_rule_loader.py        # rules/loader.py (in place)
├── test_rule_registry.py      # rules/registry.py (in place)
├── test_rule_validate.py      # rules/validate.py + main() exit code (in place)
├── test_long_conversion.py    # TODO: converters/long.py
├── test_wide_conversion.py    # TODO: converters/wide.py
├── test_factors.py            # TODO: converters/factors.py
└── test_cli.py                # TODO: cli.py
```

## TOML Organization

Use vendor-first grouping.

Do not start with separate top-level `long/` and `wide/` folders.

Reason:

- one vendor may have both shapes
- `input_shape` is already inside the TOML
- vendor grouping is the more useful first lookup key

## First Implementation Order

1. ✅ Clean the existing package structure without touching out-of-scope parts. (`a4dec1c`)
2. ✅ Add `rules/schema.py` with `pydantic` models. (`acd64d8`)
3. ✅ Add `rules/loader.py` and `rules/validate.py` (plus `rules/registry.py` and the
   `validate-rules` console script).
4. ✅ Move packaged TOMLs into `src/anndata_proteomics/parsing_rules/` — 6 vendors
   covering all software present in the test_data_download catalog. (`29c773c`)
5. ⏭️ Add generic tabular reading in `readers/`. **next**
6. Implement `converters/recognize.py`.
7. Implement `converters/long.py`.
8. Implement `converters/wide.py`.
9. Implement `converters/factors.py`.
10. Implement `converters/assemble.py`.
11. Add tests against a small curated `test_data_download` subset.

## Non-Goals For This Restart

- second-stage `obs` annotation
- experimental design joins
- condition assignment outside what is directly available in the vendor export
- downstream benchmarking metrics
- UI / notebook polish
