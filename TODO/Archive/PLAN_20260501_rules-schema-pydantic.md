# PLAN 2026-05-01 — `rules/schema.py`: pydantic models for parsing-rule TOMLs (+ JSON Schema side-output)

## Context

[docs/RESTART_PLAN.md](../docs/RESTART_PLAN.md) §"First Implementation Order" puts pydantic schema first, before any vendor-specific work. [docs/toml_schema.md](../docs/toml_schema.md) is the authoritative spec — it defines the TOML rule shape, the long-vs-wide split, factor encoding for string-valued layers, and the duplicate-handling policy.

After looking at JSON+JSON-Schema vs TOML+pydantic in this session: **stay on TOML for authoring**. Regex-heavy values, inline tables (`categories = { "MS/MS" = 1 }`), and comments make TOML strictly nicer to hand-write than JSON. The IDE-tooling argument for JSON Schema is real, but pydantic gives it back for free via `model_json_schema()` — we generate a `parse_rule.schema.json` from the same models and ship it in the package, and editors that combine TOML with JSON Schema (Even Better TOML / taplo) get live autocomplete and validation in `.toml` files.

So the model and the schema are the same artefact: pydantic is the single source of truth, JSON Schema is a derived view.

This plan covers only `rules/schema.py` plus its tests and the schema-export script. Loader, registry, validate-CLI, readers, converters, and packaged TOMLs are later plans (matching RESTART_PLAN steps 3–11).

## Scope

**In scope:**
- `src/anndata_proteomics/rules/schema.py` — pydantic models matching the TOML spec.
- Cross-field validators for the conditional rules listed in [docs/toml_schema.md](../docs/toml_schema.md) §"Formal validation".
- A small generator that emits `parse_rule.schema.json` from the pydantic models.
- Unit tests covering happy paths and each conditional rule's failure mode.

**Out of scope (later plans, do NOT do here):**
- Loading TOML from disk → that's `rules/loader.py`.
- Discovering packaged rules → that's `rules/registry.py`.
- A `validate` CLI command → that's `rules/validate.py` + `cli.py`.
- Any reader or converter logic.
- Writing actual vendor TOMLs (e.g. for DIA-NN, FragPipe). The spec examples in `toml_schema.md` are sufficient as test fixtures here.

## Files to create

```
src/anndata_proteomics/rules/__init__.py            # empty per Coding Rules
src/anndata_proteomics/rules/schema.py              # pydantic models — main deliverable
src/anndata_proteomics/rules/_export_schema.py      # generates parse_rule.schema.json from the models
tests/test_rule_models.py                           # unit tests
parsing_rules/_schema/parse_rule.schema.json        # generated artefact, committed for IDE consumption
```

The `parsing_rules/` top-level directory is created here (currently doesn't exist). Per [RESTART_PLAN.md](../docs/RESTART_PLAN.md), it lives at `src/anndata_proteomics/parsing_rules/` so it ships with the package — confirm with the user during execution if there's any reason to put it at the project root instead. Default = inside `src/anndata_proteomics/parsing_rules/_schema/`.

## Pydantic model design

One model per TOML section, plus a top-level container. Use `pydantic.BaseModel` with `model_config = ConfigDict(extra='forbid')` so unknown keys fail loudly — this is a strict spec and silent typos cost more than the convenience of permissiveness.

```python
# Sketch — final names/details to confirm against toml_schema.md during implementation

InputShape = Literal["long", "wide"]
EncodingMode = Literal["numeric", "factor"]
DuplicateMode = Literal["error", "aggregate", "keep_first", "keep_all_as_raw_table"]


class Axis(BaseModel):
    obs_keys: list[str]
    var_keys: list[str]
    x_layer: str  # must match one Layer.name; cross-checked at the top level


class Columns(BaseModel):
    obs: dict[str, str]   # internal_name -> vendor_column (or "<sample>" for wide)
    var: dict[str, str]


class Layer(BaseModel):
    name: str
    encoding_mode: EncodingMode = "numeric"
    categories: dict[str, int] | None = None
    # long-only:
    source_column: str | None = None
    # wide-only:
    column_pattern: str | None = None

    @model_validator(mode="after")
    def _factor_requires_categories(self): ...
    # source_column XOR column_pattern is enforced at ParseRule level, where input_shape is known.


class Duplicates(BaseModel):
    mode: DuplicateMode = "error"


class SampleNameCleanup(BaseModel):
    pattern: str  # may be empty per the FragPipe wide example


class ParseRule(BaseModel):
    schema_version: str
    file_version: str
    software_name: str
    software_version: str | None = None
    input_shape: InputShape
    axis: Axis
    columns: Columns
    layers: list[Layer] = Field(min_length=1)
    duplicates: Duplicates = Duplicates()
    sample_name_cleanup: SampleNameCleanup | None = None

    @model_validator(mode="after")
    def _shape_layer_consistency(self):
        # input_shape == "long":  every layer.source_column set, every layer.column_pattern None
        # input_shape == "wide":  every layer.column_pattern set, every layer.source_column None
        ...

    @model_validator(mode="after")
    def _x_layer_exists(self):
        # axis.x_layer in {l.name for l in layers}
        ...

    @model_validator(mode="after")
    def _cleanup_only_for_wide(self):
        # sample_name_cleanup must be None when input_shape == "long"
        ...
```

### Conditional rules to enforce ([toml_schema.md](../docs/toml_schema.md) §"Formal validation")

1. `input_shape = "long"` → every layer must have `source_column`, none may have `column_pattern`.
2. `input_shape = "wide"` → every layer must have `column_pattern`, none may have `source_column`.
3. `encoding_mode = "factor"` → `categories` required (non-empty dict of str→int).
4. `axis.x_layer` must equal one of `layers[*].name`.
5. `sample_name_cleanup` is only valid when `input_shape = "wide"` (it has no meaning for long).

### What we deliberately do NOT validate at this layer

- Whether `axis.obs_keys` / `axis.var_keys` reference real entries in `columns.obs` / `columns.var`. The spec is ambiguous about whether these are LHS internal names or RHS vendor names (the DIA-NN example uses RHS vendor names for `var_keys` but the FragPipe example uses LHS for `obs_keys`). Defer cross-checking to the converter, where we also have the actual DataFrame columns to compare against.
- Whether `column_pattern` is a valid regex — pydantic doesn't compile it; the converter will, and that's the natural place for a clearer error.
- Whether vendor columns actually exist in the data — that's a converter-time check, not a schema-time check.

This separation keeps `schema.py` purely about TOML structural correctness; semantic checks against actual data live in `converters/`.

## JSON Schema export

`src/anndata_proteomics/rules/_export_schema.py` — a tiny module that:

```python
import json
from pathlib import Path
from anndata_proteomics.rules.schema import ParseRule

def main() -> None:
    schema = ParseRule.model_json_schema()
    out = Path(__file__).parent.parent / "parsing_rules" / "_schema" / "parse_rule.schema.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(schema, indent=2) + "\n")

if __name__ == "__main__":
    main()
```

Add a `[project.scripts]` entry to [pyproject.toml](../pyproject.toml):

```toml
[project.scripts]
export-rule-schema = "anndata_proteomics.rules._export_schema:main"
```

The committed `parse_rule.schema.json` becomes the contract for editors. We'll add a CI check (later, not this plan) that fails if the file is out of sync with the models.

## Tests — `tests/test_rule_models.py`

Test fixtures: copy the two examples in [toml_schema.md](../docs/toml_schema.md) §"Long example" and §"Wide example" verbatim, parse with `tomllib`, feed the dict to `ParseRule.model_validate(...)`. They must round-trip cleanly.

Negative cases (one focused test each):

- Long rule with a layer missing `source_column` → ValidationError mentioning the layer name.
- Long rule with a layer that has `column_pattern` set → ValidationError.
- Wide rule with a layer missing `column_pattern` → ValidationError.
- Wide rule with a layer that has `source_column` set → ValidationError.
- `encoding_mode = "factor"` without `categories` → ValidationError.
- `axis.x_layer = "DoesNotExist"` → ValidationError.
- `duplicates.mode = "wrong"` → ValidationError on the literal.
- Top-level unknown key (e.g. `foo = "bar"`) → ValidationError (because `extra='forbid'`).
- `sample_name_cleanup` set on a long rule → ValidationError.

One sanity test for the export script: call `ParseRule.model_json_schema()` and assert the result has the expected top-level `properties` keys (`schema_version`, `input_shape`, `axis`, `columns`, `layers`, `duplicates`).

## pyproject.toml additions

- Add `pydantic` to `[project] dependencies` (if not already pulled transitively — confirm during execution).
- Add the `export-rule-schema` console script (see above).
- `tomllib` is stdlib on Python 3.11+; keep `requires-python = ">=3.9"` only if we plan to use `tomli` as a fallback. **Recommendation: bump `requires-python` to `>=3.11`** since (a) we're starting fresh and (b) every supported Python is on 3.11+ now. Confirm with user during execution if they'd rather stay on 3.9 with `tomli`.

## Verification

```bash
cd /Users/wolski/projects/anndata_bridge/anndata_proteomics_bridge
source .venv/bin/activate
uv pip install -e .

# 1. Models import cleanly
.venv/bin/python -c "from anndata_proteomics.rules.schema import ParseRule; print(ParseRule.__name__)"

# 2. Tests pass
.venv/bin/python -m pytest tests/ -q

# 3. Round-trip the spec examples
.venv/bin/python - <<'PY'
import tomllib
from anndata_proteomics.rules.schema import ParseRule
LONG = """<paste long example from toml_schema.md>"""
WIDE = """<paste wide example from toml_schema.md>"""
for src in (LONG, WIDE):
    rule = ParseRule.model_validate(tomllib.loads(src))
    print(rule.software_name, rule.input_shape, len(rule.layers))
PY

# 4. JSON Schema export works and produces a non-trivial file
.venv/bin/python -m anndata_proteomics.rules._export_schema
ls -la src/anndata_proteomics/parsing_rules/_schema/parse_rule.schema.json
.venv/bin/python -c "import json; s = json.load(open('src/anndata_proteomics/parsing_rules/_schema/parse_rule.schema.json')); print(sorted(s['properties']))"

# 5. Generated schema lints (optional sanity check)
.venv/bin/python -c "import jsonschema; jsonschema.Draft202012Validator.check_schema(__import__('json').load(open('src/anndata_proteomics/parsing_rules/_schema/parse_rule.schema.json')))"
```

## After this plan lands

Next plan (RESTART_PLAN step 3): `rules/loader.py` — TOML file → validated `ParseRule`, plus `rules/validate.py` for batch-validation of packaged rules. Then step 4 (move first concrete TOMLs into `parsing_rules/<vendor>/`), then readers (step 5).

## Decisions resolved (2026-05-01)

1. **`parsing_rules/` location** → inside `src/anndata_proteomics/parsing_rules/` (ships with the package).
2. **Python floor** → bump to `>=3.13`. Use stdlib `tomllib`. Update `[tool.black] target-version` and `[tool.ruff] target-version` accordingly.
3. **Generated JSON Schema** → committed at `src/anndata_proteomics/parsing_rules/_schema/parse_rule.schema.json`.
