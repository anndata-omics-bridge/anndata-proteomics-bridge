"""Validate effective parsing-rule JSON against the generated JSON Schema.

This is **structural-parity only** — pydantic remains the source of truth for
cross-field rules ("long → every layer has source", "factor encoding
requires categories", etc.), which JSON Schema cannot express. JSON Schema
covers only types, literals, required fields, and additionalProperties.
See docs/ARCHITECTURE.md.
"""

from __future__ import annotations

import json

import jsonschema
import pytest

from anndata_proteomics.rules.loader import load_rule
from anndata_proteomics.rules.registry import (
    RuleLocator,
    iter_packaged_rules,
    packaged_rules_root,
)
from anndata_proteomics.rules.schema import ParseRule, ParseRuleDocument

SCHEMA_ROOT = packaged_rules_root() / "_schema"
SCHEMA_PATH = SCHEMA_ROOT / "parse_rule.schema.json"


_VALID_LONG = {
    "schema_version": "0.1",
    "file_version": "1",
    "software_name": "Fake",
    "software_version": "^1$",
    "input_shape": "long",
    "quantification_level": "ion",
    "axis": {
        "obs_keys": ["Run"],
        "var_keys": ["Foo"],
        "x_layer": "X",
        "duplicates": {"mode": "error"},
    },
    "columns": {
        "obs": {"select": {"Run": "Run"}},
        "var": {"select": {"Foo": "Foo"}},
    },
    "layers": [{"name": "X", "source": "Foo"}],
}


def _load_schema() -> dict[str, object]:
    return json.loads(SCHEMA_PATH.read_text())


def _locator_id(locator: RuleLocator) -> str:
    return f"{locator.path.parent.name}/{locator.level}"


def test_exported_schema_is_valid_draft_2020_12() -> None:
    """The generated parse_rule.schema.json must itself be a well-formed JSON Schema."""
    jsonschema.Draft202012Validator.check_schema(_load_schema())


@pytest.mark.parametrize(
    ("filename", "model"),
    [
        ("parse_rule.schema.json", ParseRule),
        ("parse_rule_document.schema.json", ParseRuleDocument),
    ],
)
def test_packaged_json_schema_matches_models(
    filename: str,
    model: type[ParseRule] | type[ParseRuleDocument],
) -> None:
    """The packaged schemas are exact semantic exports of their Pydantic models."""
    expected = json.loads((SCHEMA_ROOT / filename).read_text(encoding="utf-8"))
    assert model.model_json_schema() == expected


@pytest.mark.parametrize(
    "locator",
    list(iter_packaged_rules()),
    ids=_locator_id,
)
def test_packaged_rule_passes_json_schema(locator: RuleLocator) -> None:
    """Every effective document level must validate against the JSON Schema."""
    data = load_rule(locator).model_dump(by_alias=True, mode="json")
    jsonschema.validate(instance=data, schema=_load_schema())


def test_baseline_json_is_valid() -> None:
    """Sanity check: the baseline used by the negative tests is itself valid."""
    jsonschema.validate(instance=_VALID_LONG, schema=_load_schema())


def test_json_schema_rejects_missing_required_field() -> None:
    data = {key: value for key, value in _VALID_LONG.items() if key != "quantification_level"}
    with pytest.raises(jsonschema.ValidationError, match="quantification_level"):
        jsonschema.validate(instance=data, schema=_load_schema())


def test_json_schema_rejects_unknown_top_level_key() -> None:
    data = {**_VALID_LONG, "foo": "bar"}
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(instance=data, schema=_load_schema())


def test_json_schema_rejects_invalid_literal() -> None:
    data = json.loads(json.dumps(_VALID_LONG))
    data["axis"]["duplicates"]["mode"] = "wrong"
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(instance=data, schema=_load_schema())
