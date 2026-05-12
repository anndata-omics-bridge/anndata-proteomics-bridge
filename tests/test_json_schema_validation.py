"""Round-trip validate parsed TOMLs against the generated parse_rule.schema.json.

This is **structural-parity only** — pydantic remains the source of truth for
cross-field rules ("long → every layer has source_column", "factor encoding
requires categories", etc.), which JSON Schema cannot express. JSON Schema
covers only types, literals, required fields, and additionalProperties.
See docs/ARCHITECTURE.md.
"""

from __future__ import annotations

import json
import tomllib
from pathlib import Path

import jsonschema
import pytest

from anndata_proteomics.rules.registry import iter_packaged_rules, packaged_rules_root


SCHEMA_PATH = packaged_rules_root() / "_schema" / "parse_rule.schema.json"


# Minimal valid long-format TOML used as baseline for negative-case mutations.
_VALID_LONG_TOML = """
schema_version = "0.1"
file_version = "1"
software_name = "Fake"
input_shape = "long"
quantification_level = "ion"

[axis]
obs_keys = ["Run"]
var_keys = ["Foo"]
x_layer = "X"

[axis.duplicates]
mode = "error"

[columns.obs.select]
Run = "Run"

[columns.var.select]
Foo = "Foo"

[[layers]]
name = "X"
source_column = "Foo"
"""


def _load_schema() -> dict:
    return json.loads(SCHEMA_PATH.read_text())


def test_exported_schema_is_valid_draft_2020_12() -> None:
    """The generated parse_rule.schema.json must itself be a well-formed JSON Schema."""
    jsonschema.Draft202012Validator.check_schema(_load_schema())


@pytest.mark.parametrize(
    "toml_path",
    list(iter_packaged_rules()),
    ids=lambda p: f"{p.parent.name}/{p.name}",
)
def test_packaged_rule_passes_json_schema(toml_path: Path) -> None:
    """Every packaged TOML must validate against the generated JSON Schema."""
    data = tomllib.loads(toml_path.read_text())
    jsonschema.validate(instance=data, schema=_load_schema())


def test_baseline_toml_is_valid() -> None:
    """Sanity check: the baseline used by the negative tests is itself valid."""
    data = tomllib.loads(_VALID_LONG_TOML)
    jsonschema.validate(instance=data, schema=_load_schema())


def test_json_schema_rejects_missing_required_field() -> None:
    bad = _VALID_LONG_TOML.replace('quantification_level = "ion"\n', "")
    data = tomllib.loads(bad)
    with pytest.raises(jsonschema.ValidationError, match="quantification_level"):
        jsonschema.validate(instance=data, schema=_load_schema())


def test_json_schema_rejects_unknown_top_level_key() -> None:
    bad = _VALID_LONG_TOML + '\nfoo = "bar"\n'
    data = tomllib.loads(bad)
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(instance=data, schema=_load_schema())


def test_json_schema_rejects_invalid_literal() -> None:
    bad = _VALID_LONG_TOML.replace('mode = "error"', 'mode = "wrong"')
    data = tomllib.loads(bad)
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(instance=data, schema=_load_schema())
