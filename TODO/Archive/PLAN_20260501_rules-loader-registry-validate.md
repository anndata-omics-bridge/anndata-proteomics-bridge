# PLAN 2026-05-01 — `rules/loader.py` + `rules/registry.py` + `rules/validate.py` (RESTART_PLAN step 3)

## Context

Step 3 of [docs/RESTART_PLAN.md](../docs/RESTART_PLAN.md): factor the in-line `tomllib.loads(...) + ParseRule.model_validate(...)` walk that currently lives in [tests/test_packaged_rules.py](../tests/test_packaged_rules.py) into reusable library code, so:

- Future converters and the eventual `cli.py` have **one** documented entrypoint to load a TOML rule.
- A `python -m anndata_proteomics.rules.validate` command (and a `validate-rules` console script) can sanity-check every packaged rule with friendly errors — useful in CI and for a contributor adding a new vendor TOML.
- Discovery (`iter_packaged_rules`, `find_rule(software, level, version)`) is split from parsing/validation (`load_rule`) because RESTART_PLAN explicitly separates `loader.py` and `registry.py` (file → ParseRule vs (software, level, version) → path).

The current `ParseRule` schema in [src/anndata_proteomics/rules/schema.py](../src/anndata_proteomics/rules/schema.py) is the single source of truth. This plan adds *no* schema changes.

## Files to create

```
src/anndata_proteomics/rules/loader.py        load_rule(path) -> ParseRule, plus packaged convenience
src/anndata_proteomics/rules/registry.py      iter_packaged_rules / find_rule / RuleNotFound
src/anndata_proteomics/rules/validate.py      ValidationResult, validate_file, validate_all_packaged, main()
tests/test_rule_loader.py                     happy + error paths for loader
tests/test_rule_registry.py                   happy + error paths for registry
tests/test_rule_validate.py                   ValidationResult and the walk; main() exit code
```

## Files to modify

- [tests/test_packaged_rules.py](../tests/test_packaged_rules.py) — replace inline `tomllib + model_validate` with `validate_all_packaged()`; keep filename ↔ field test.
- [pyproject.toml](../pyproject.toml) — add `validate-rules = "anndata_proteomics.rules.validate:main"`.
- [docs/RESTART_PLAN.md](../docs/RESTART_PLAN.md) — flip step 3 ✅ once landed.

## API summary

`loader.load_rule(path) -> ParseRule` — TOML → validated rule, FileNotFoundError on missing path, ValidationError with file path attached on bad data.

`loader.load_packaged_rule(software, quantification_level, file_version="1") -> ParseRule` — sugar over `find_rule + load_rule`.

`registry.packaged_rules_root() -> Path` — uses `importlib.resources.files("anndata_proteomics") / "parsing_rules"`.

`registry.iter_packaged_rules() -> Iterator[Path]` — sorted glob of `*/parse_*.toml`.

`registry.find_rule(software, level, file_version="1") -> Path` — resolves to `parsing_rules/<software>/parse_<software>_<level>_<file_version>.toml`. Raises `RuleNotFound` (subclass of `LookupError`) listing what's available in the vendor folder.

`validate.ValidationResult` — frozen dataclass `{path, ok, error, rule}`.

`validate.validate_file(path) -> ValidationResult` — never raises.

`validate.validate_all_packaged() -> list[ValidationResult]`.

`validate.main(argv=None) -> int` — prints `PASS path` / `FAIL path: msg` per rule, then summary `N packaged rule(s) checked, K failed.`, returns 0 / 1.

## Verification

```bash
.venv/bin/python -c "from anndata_proteomics.rules.loader import load_packaged_rule; r = load_packaged_rule('wombat', 'peptidoform'); print(r.software_name, r.quantification_level, r.input_shape, len(r.layers))"
# WOMBAT peptidoform wide 2

validate-rules    # 6 PASS, exit 0
.venv/bin/python -m pytest tests/ -q
```

## Out of scope

- User-facing CLI surface (`cli.py`) — later.
- `readers/` (step 5).
- `proteome_discoverer/` packaged TOML (step 4 finishing, not 3).
- Recursive vendor subdirectories — convention is one level.
