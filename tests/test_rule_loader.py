"""Tests for self-contained parsing-rule documents."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from pydantic import ValidationError

from anndata_proteomics.rules.loader import (
    RuleDocumentError,
    load_packaged_rule,
    load_rule,
    load_rule_document,
    parse_rule_document,
)
from anndata_proteomics.rules.registry import find_rule
from anndata_proteomics.rules.schema import ParseRuleDocument, _merge_rule_dicts


def test_merge_scalar_level_wins() -> None:
    assert _merge_rule_dicts({"a": 1, "b": 2}, {"b": 3}) == {"a": 1, "b": 3}


def test_merge_objects_deep_merge() -> None:
    base = {"axis": {"obs_keys": ["Run"], "duplicates": {"mode": "error"}}}
    level = {"axis": {"var_keys": ["Ion"], "x_layer": "X"}}
    assert _merge_rule_dicts(base, level) == {
        "axis": {
            "obs_keys": ["Run"],
            "duplicates": {"mode": "error"},
            "var_keys": ["Ion"],
            "x_layer": "X",
        }
    }


def test_merge_object_arrays_append_and_scalar_arrays_replace() -> None:
    assert _merge_rule_dicts(
        {"compute": [{"name": "peptidoform"}]},
        {"compute": [{"name": "ion"}]},
    )["compute"] == [{"name": "peptidoform"}, {"name": "ion"}]
    assert _merge_rule_dicts({"var_keys": ["A"]}, {"var_keys": ["B"]}) == {"var_keys": ["B"]}


def _document() -> dict:
    return {
        "schema_version": "0.1",
        "file_version": "1",
        "software_name": "MyVendor",
        "software_version": "^1$",
        "base": {
            "input_shape": "long",
            "axis": {"obs_keys": ["Run"], "duplicates": {"mode": "error"}},
            "columns": {"obs": {"select": {"Run": "Run"}}},
        },
        "levels": {
            "ion": {
                "axis": {"var_keys": ["Ion"], "x_layer": "Intensity"},
                "columns": {"var": {"select": {"Ion": "Ion"}}},
                "layers": [{"name": "Intensity", "source": "Intensity"}],
            }
        },
    }


def _write_document(tmp_path: Path, data: dict | None = None) -> Path:
    path = tmp_path / "rules.json"
    path.write_text(json.dumps(data or _document()))
    return path


def test_parse_document_returns_pydantic_source_model() -> None:
    document = parse_rule_document(json.dumps(_document()))
    assert isinstance(document, ParseRuleDocument)
    assert list(document.levels) == ["ion"]


def test_effective_rule_merges_base_and_level() -> None:
    rule = parse_rule_document(json.dumps(_document())).effective_rule("ion")
    assert rule.software_name == "MyVendor"
    assert rule.axis.obs_keys == ["Run"]
    assert rule.axis.var_keys == ["Ion"]


def test_load_rule_path_without_level_works_for_single_level(tmp_path: Path) -> None:
    assert load_rule(_write_document(tmp_path)).quantification_level == "ion"


def test_load_rule_path_requires_level_for_multi_level_document(tmp_path: Path) -> None:
    data = _document()
    data["levels"]["protein"] = data["levels"]["ion"]
    path = _write_document(tmp_path, data)
    with pytest.raises(RuleDocumentError, match="select one explicitly"):
        load_rule(path)


def test_load_rule_locator_selects_level() -> None:
    locator = find_rule("diann", "protein", "1.9.2")
    assert load_rule(locator).quantification_level == "protein"


def test_load_rule_document_validates_every_level(tmp_path: Path) -> None:
    data = _document()
    data["levels"]["protein"] = {
        "axis": {"var_keys": ["Protein"], "x_layer": "missing"},
        "columns": {"var": {"select": {"Protein": "Protein"}}},
        "layers": [{"name": "Abundance", "source": "Abundance"}],
    }
    path = _write_document(tmp_path, data)
    with pytest.raises(ValidationError, match="x_layer") as error:
        load_rule_document(path)
    assert any("level: protein" in note for note in error.value.__notes__)


def test_old_extends_format_is_rejected() -> None:
    data = _document()
    data["$extends"] = "base.json"
    with pytest.raises(ValidationError, match="extra_forbidden"):
        parse_rule_document(json.dumps(data))


def test_malformed_json_attaches_source_path(tmp_path: Path) -> None:
    path = tmp_path / "bad.json"
    path.write_text('{"broken": }')
    with pytest.raises(json.JSONDecodeError) as error:
        load_rule_document(path)
    assert any(str(path) in note for note in error.value.__notes__)


def test_duplicate_key_and_nonstandard_constant_are_rejected() -> None:
    with pytest.raises(RuleDocumentError, match="duplicate JSON key"):
        parse_rule_document('{"base": {}, "base": {}}')
    with pytest.raises(RuleDocumentError, match="non-standard JSON"):
        parse_rule_document('{"schema_version": NaN}')


def test_load_packaged_rule_uses_existing_version_groups() -> None:
    assert load_packaged_rule("diann", "fragment", "1.9.2").software_version == "^1\\..*"
    with pytest.raises(ValueError, match="no packaged rule"):
        load_packaged_rule("diann", "fragment", "2.3.0")
