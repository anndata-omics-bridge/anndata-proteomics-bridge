# PLAN 2026-05-02 — JSON Schema round-trip tests + `anndata-proteomics` CLI

## Context

Two changes that fall out of the same conversation:

1. **JSON Schema round-trip validation.** Today we generate `parsing_rules/_schema/parse_rule.schema.json` from the pydantic models (`rules/_export_schema.py`) and only consume it in IDE tooling. Nothing in our test suite parses a TOML and runs it through `jsonschema.validate(...)`. That leaves two real (if narrow) gaps unguarded:
   - The exported JSON Schema could be malformed and we'd never notice — IDE autocomplete would silently break.
   - Pydantic and the JSON Schema could disagree on the *structural* checks (types, literals, required fields, `additionalProperties: false`), which would mean editors see a different set of valid TOMLs than what pydantic actually accepts.
   Pydantic stays strictly stronger because per-layer cross-field rules ("long → every layer has `source_column`", "factor encoding → non-empty `categories`", "x_layer matches a layer name") are not expressible in JSON Schema. So this plan does *not* try to enforce parity on those — only the structural rules.

2. **User-facing CLI** (`anndata-proteomics ...`). Today we have two narrow console scripts (`validate-rules`, `export-rule-schema`) that both walk packaged rules. There is no command to validate a single user-provided TOML, no command to list the packaged rules, and no place to land the future `convert` command. The user explicitly wants option B from the earlier discussion: a single umbrella CLI living under `src/anndata_proteomics/scripts/`, ready to grow as readers/converters land.

Both changes are small, additive, and don't touch the production validation path. After this plan: the schema export is sanity-checked, and `anndata-proteomics validate new.toml` works.

## Part 1 — JSON Schema round-trip validation

### Files

- **`tests/test_json_schema_validation.py`** (new) — three test groups described below.
- **`pyproject.toml`** — add `jsonschema>=4` to `[project.optional-dependencies] dev`.

### Tests

```python
# tests/test_json_schema_validation.py

SCHEMA_PATH = packaged_rules_root() / "_schema" / "parse_rule.schema.json"

def _load_schema() -> dict: ...

# 1. Well-formedness
def test_exported_schema_is_valid_draft_2020_12():
    """The generated parse_rule.schema.json must itself be a valid JSON Schema."""
    jsonschema.Draft202012Validator.check_schema(_load_schema())

# 2. Round-trip parity for packaged TOMLs (happy path)
@pytest.mark.parametrize("toml_path", list(iter_packaged_rules()), ids=...)
def test_packaged_rule_passes_json_schema(toml_path):
    """Every packaged TOML, parsed by tomllib, must validate against the JSON Schema."""
    data = tomllib.loads(toml_path.read_text())
    jsonschema.validate(instance=data, schema=_load_schema())

# 3. Structural negatives the JSON Schema *should* still catch
def test_json_schema_rejects_missing_required_field(tmp_path):
    """A valid TOML missing quantification_level must fail JSON Schema validation."""
    bad = tmp_path / "missing.toml"
    bad.write_text(LONG_EXAMPLE_MINUS_QUANT_LEVEL)
    data = tomllib.loads(bad.read_text())
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(instance=data, schema=_load_schema())

def test_json_schema_rejects_unknown_top_level_key(tmp_path):
    """extra='forbid' must round-trip to additionalProperties: false."""
    bad = tmp_path / "extra.toml"
    bad.write_text(LONG_EXAMPLE + '\nfoo = "bar"\n')
    data = tomllib.loads(bad.read_text())
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(instance=data, schema=_load_schema())

def test_json_schema_rejects_invalid_literal(tmp_path):
    """duplicates.mode = 'wrong' must fail because the JSON Schema enumerates allowed values."""
    bad = tmp_path / "wrong_mode.toml"
    bad.write_text(LONG_EXAMPLE.replace('mode = "error"', 'mode = "wrong"'))
    data = tomllib.loads(bad.read_text())
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(instance=data, schema=_load_schema())
```

Reuse `LONG_EXAMPLE` / `WIDE_EXAMPLE` constants from [tests/test_rule_models.py](../tests/test_rule_models.py) by importing them, so we don't duplicate fixtures.

### Documentation

- One short paragraph in [docs/ARCHITECTURE.md](../docs/ARCHITECTURE.md) under `rules/_export_schema.py` noting that the export is round-trip-checked by `tests/test_json_schema_validation.py` (structural parity only; pydantic stays the only source of truth for cross-field rules).

## Part 2 — `anndata-proteomics` CLI under `src/anndata_proteomics/scripts/`

### Files

- **`src/anndata_proteomics/scripts/__init__.py`** — empty per Coding Rules.
- **`src/anndata_proteomics/scripts/cli.py`** — single `cyclopts.App` named `anndata-proteomics`, with the subcommands below.
- **`tests/test_cli.py`** — invoke the cyclopts `App` programmatically, capture stdout, assert exit codes and output.
- **`pyproject.toml`** — add `anndata-proteomics = "anndata_proteomics.scripts.cli:main"` to `[project.scripts]`.
- **`docs/ARCHITECTURE.md`** — add a `scripts/` section listing the subcommands and their dispatchers.
- **`README.md`** — bump the Quickstart to use `anndata-proteomics validate <path>` as the primary entrypoint; keep `validate-rules` mentioned as the CI-friendly bulk-walk alias.

### Subcommands

```
anndata-proteomics validate [path ...]    Validate one or more TOMLs.
                                          With no path: walks all packaged rules
                                          (= current `validate-rules` behavior).
anndata-proteomics list                   List packaged rules: vendor, software,
                                          quantification_level, file_version, path.
anndata-proteomics export-schema          Regenerate parse_rule.schema.json
                                          (= current `export-rule-schema`).
anndata-proteomics convert <data> <toml>  STUB. Prints "not yet implemented;
                                          see RESTART_PLAN.md steps 5-10",
                                          returns exit 2.
```

`cyclopts` is already a dependency.

### Implementation sketch

```python
# src/anndata_proteomics/scripts/cli.py

from pathlib import Path
from cyclopts import App
from anndata_proteomics.rules import _export_schema
from anndata_proteomics.rules.registry import iter_packaged_rules
from anndata_proteomics.rules.validate import (
    validate_all_packaged, validate_file, ValidationResult
)

app = App(name="anndata-proteomics", help="anndata_proteomics CLI")


def _print_and_exit_code(results: list[ValidationResult]) -> int:
    """Shared PASS/FAIL printer used by validate(); same output format as
    `validate-rules`, so users see one consistent result format."""
    ...


@app.command
def validate(*paths: Path) -> int:
    """Validate one or more TOML rule files; defaults to all packaged."""
    if not paths:
        results = validate_all_packaged()
    else:
        results = [validate_file(p) for p in paths]
    return _print_and_exit_code(results)


@app.command(name="list")
def list_rules() -> int:
    """List packaged parsing rules."""
    for p in iter_packaged_rules():
        # parse filename: parse_<software>_<level>_<version>.toml
        parts = p.stem.split("_")
        software, level, version = parts[1], parts[-2], parts[-1]
        print(f"{software:14}  {level:12}  v{version:3}  {p}")
    return 0


@app.command(name="export-schema")
def export_schema_cmd() -> int:
    _export_schema.main()
    return 0


@app.command
def convert(data: Path, rule_toml: Path) -> int:
    """STUB until readers/ + converters/ land."""
    print(
        f"error: convert is not yet implemented; "
        f"see docs/RESTART_PLAN.md steps 5-10"
    )
    return 2


def main() -> int:
    return app() or 0
```

Notes:
- The shared `_print_and_exit_code` is the same logic currently inside `rules.validate.main`. We factor it out so both that older entry and the new cli.py use the same formatter — single source of formatting truth, avoids drift.
- `_export_schema.main()` already exists; cli.py just calls it.
- `convert` returns exit code 2 (≠ 0 and ≠ 1) so a CI script can distinguish "not implemented yet" from "validation failed" if it cares.

### What stays where

- **`validate-rules` and `export-rule-schema` console scripts stay** for now. Reasons:
  - They are tested and used. Removing them would break any CI line that already uses them.
  - The new `anndata-proteomics validate` and `anndata-proteomics export-schema` are a *superset*, not a replacement.
  - When `cli.py` matures (after `convert` lands for real), we revisit deprecating the narrow scripts.
- **`rules.validate.main()` stays** as the implementation behind `validate-rules`; cli.py calls into the shared formatter. No refactor that moves files between `rules/` and `scripts/`.

### Tests — `tests/test_cli.py`

```python
def test_cli_validate_no_args_walks_packaged(capsys):
    rc = app(["validate"], exit_on_error=False)
    out = capsys.readouterr().out
    assert rc == 0
    assert "PASS" in out
    assert "0 failed" in out

def test_cli_validate_single_path_happy(capsys):
    path = find_rule("diann", "ion")
    rc = app(["validate", str(path)], exit_on_error=False)
    out = capsys.readouterr().out
    assert rc == 0
    assert "PASS" in out

def test_cli_validate_single_path_bad(tmp_path, capsys):
    bad = tmp_path / "bad.toml"
    bad.write_text("not = valid [[")
    rc = app(["validate", str(bad)], exit_on_error=False)
    out = capsys.readouterr().out
    assert rc == 1
    assert "FAIL" in out

def test_cli_list_shows_six_rules(capsys):
    rc = app(["list"], exit_on_error=False)
    out = capsys.readouterr().out
    assert rc == 0
    assert "diann" in out and "wombat" in out
    assert out.count("\n") == 6

def test_cli_export_schema_writes_file(tmp_path, monkeypatch, capsys):
    # Easiest: just call and verify the file exists and round-trips.
    rc = app(["export-schema"], exit_on_error=False)
    assert rc == 0
    assert (packaged_rules_root() / "_schema" / "parse_rule.schema.json").exists()

def test_cli_convert_is_stub(capsys):
    rc = app(["convert", "data.tsv", "rule.toml"], exit_on_error=False)
    out = capsys.readouterr().out
    assert rc == 2
    assert "not yet implemented" in out.lower()
```

(Exact cyclopts invocation pattern depends on its API for in-process testing — I'll verify the right form when implementing; the principle is "call the App with argv list, capture stdout, check exit code".)

## Order of operations

1. Part 1 first (JSON Schema tests + `jsonschema` dev dep). Self-contained, no new public surface.
2. Part 2 second (CLI scaffold). Depends on Part 1 having clean tests as a baseline.
3. One commit per part:
   - `test: round-trip JSON Schema validation for packaged TOMLs`
   - `feat(scripts): anndata-proteomics CLI with validate / list / export-schema / convert (stub)`
4. Update `docs/RESTART_PLAN.md` only at the end of Part 2: note that the user-facing `cli.py` skeleton is in place (with `convert` still a stub).

## Verification

```bash
cd /Users/wolski/projects/anndata_bridge/anndata_proteomics_bridge
source .venv/bin/activate
uv pip install -e '.[dev]'

# Part 1 — JSON Schema tests
.venv/bin/python -m pytest tests/test_json_schema_validation.py -q

# Part 2 — full suite + CLI smoke
.venv/bin/python -m pytest tests/ -q
anndata-proteomics validate                                  # walks packaged, exit 0
anndata-proteomics validate src/anndata_proteomics/parsing_rules/diann/parse_diann_ion_1.toml
anndata-proteomics list
anndata-proteomics export-schema
anndata-proteomics convert dummy.tsv rule.toml; echo $?       # exit 2
```

## Out of scope

- No conversion logic — `convert` is a stub.
- No deprecation of `validate-rules` / `export-rule-schema`. Revisit once `convert` is real.
- No refactor of `rules/validate.py` beyond extracting the print-and-exit formatter into a small shared helper used by both `rules.validate.main()` and the new `scripts.cli.validate`.
- Cross-field-rule parity between pydantic and JSON Schema. Pydantic stays strictly stronger; that's the documented contract.
- No `scripts/` subfolder explosion — one `cli.py` for all subcommands until the surface justifies splitting (rough rule: ≥5 commands or any single command >100 LOC).
