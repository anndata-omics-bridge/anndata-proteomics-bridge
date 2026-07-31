"""Tests for self-contained parsing-rule documents."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from pydantic import ValidationError

from anndata_proteomics.rules import loader
from anndata_proteomics.rules.loader import (
    RuleDocumentError,
    load_packaged_rule_for_version,
    load_rule,
    load_rule_document,
    load_single_level_rule,
    parse_rule_document,
)
from anndata_proteomics.rules.registry import RuleLocator, find_rule_for_version
from anndata_proteomics.rules.schema import ParseRuleDocument, RuleCompositionError


def _document() -> dict[str, object]:
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


def _write_document(
    tmp_path: Path,
    data: dict[str, object] | None = None,
) -> Path:
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


def test_effective_rule_composes_typed_axis_and_column_fragments() -> None:
    data = _document()
    base = data["base"]
    levels = data["levels"]
    assert isinstance(base, dict)
    assert isinstance(levels, dict)
    ion = levels["ion"]
    assert isinstance(ion, dict)
    base["axis"] = {
        "obs_keys": ["Run"],
        "var_keys": ["BaseIon"],
        "duplicates": {"mode": "error"},
    }
    base["columns"] = {
        "obs": {"select": {"Run": "Run"}},
        "var": {"select": {"BaseIon": "BaseIon"}},
    }
    ion["axis"] = {
        "var_keys": ["Ion"],
        "x_layer": "Intensity",
        "duplicates": {"mode": "aggregate"},
    }

    rule = parse_rule_document(json.dumps(data)).effective_rule("ion")

    assert rule.axis.obs_keys == ["Run"]
    assert rule.axis.var_keys == ["Ion"]
    assert rule.axis.duplicates.mode == "aggregate"
    assert rule.columns.var.select == {"BaseIon": "BaseIon", "Ion": "Ion"}


def test_effective_rule_appends_typed_layer_objects() -> None:
    data = _document()
    base = data["base"]
    assert isinstance(base, dict)
    base["layers"] = [{"name": "Quality", "source": "Quality"}]

    rule = parse_rule_document(json.dumps(data)).effective_rule("ion")

    assert [layer.name for layer in rule.layers] == ["Quality", "Intensity"]


@pytest.mark.parametrize(
    ("field", "expected", "expected_path"),
    [
        ("input_shape", "requires input_shape", ("input_shape",)),
        ("axis.var_keys", "requires axis.var_keys", ("axis", "var_keys")),
        ("layers", "requires at least one layer", ("layers",)),
    ],
)
def test_explicit_null_or_empty_required_fragments_raise_composition_errors(
    field: str,
    expected: str,
    expected_path: tuple[str, ...],
) -> None:
    data = _document()
    base = data["base"]
    levels = data["levels"]
    assert isinstance(base, dict)
    assert isinstance(levels, dict)
    ion = levels["ion"]
    assert isinstance(ion, dict)
    if field == "input_shape":
        base["input_shape"] = None
    elif field == "axis.var_keys":
        axis = ion["axis"]
        assert isinstance(axis, dict)
        axis["var_keys"] = None
    else:
        ion["layers"] = []

    with pytest.raises(RuleCompositionError, match=expected) as captured:
        parse_rule_document(json.dumps(data))
    assert captured.value.path == expected_path


def test_effective_rule_appends_typed_modification_map_entries() -> None:
    data = _document()
    base = data["base"]
    levels = data["levels"]
    assert isinstance(base, dict)
    assert isinstance(levels, dict)
    ion = levels["ion"]
    assert isinstance(ion, dict)
    base["modifications"] = {
        "parser": "token_regex",
        "source_column": "Ion",
        "token_pattern": "\\[([^]]+)\\]",
        "case_sensitive": True,
        "map": [{"token": "Oxidation", "accession": "UNIMOD:35"}],
    }
    ion["modifications"] = {
        "parser": "token_regex",
        "source_column": "Ion",
        "token_pattern": "\\[([^]]+)\\]",
        "unknown_policy": "error",
        "map": [{"token": "Acetyl", "accession": "UNIMOD:1"}],
    }

    modifications = parse_rule_document(json.dumps(data)).effective_rule("ion").modifications

    assert modifications is not None
    assert modifications.case_sensitive
    assert modifications.unknown_policy == "error"
    assert [entry.token for entry in modifications.map] == ["Oxidation", "Acetyl"]


def test_load_rule_path_without_level_works_for_single_level(tmp_path: Path) -> None:
    assert load_single_level_rule(_write_document(tmp_path)).quantification_level == "ion"


def test_load_rule_path_requires_level_for_multi_level_document(tmp_path: Path) -> None:
    data = _document()
    levels = data["levels"]
    assert isinstance(levels, dict)
    levels["protein"] = levels["ion"]
    path = _write_document(tmp_path, data)
    with pytest.raises(RuleDocumentError, match="select one explicitly"):
        load_single_level_rule(path)


def test_load_rule_locator_selects_level() -> None:
    locator = find_rule_for_version("diann", "protein", "1.9.2")
    assert load_rule(locator).quantification_level == "protein"


def test_load_rule_document_validates_every_level(tmp_path: Path) -> None:
    data = _document()
    levels = data["levels"]
    assert isinstance(levels, dict)
    levels["protein"] = {
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
    assert (
        load_packaged_rule_for_version("diann", "fragment", "1.9.2").software_version == "^1\\..*"
    )
    with pytest.raises(ValueError, match="no packaged rule"):
        load_packaged_rule_for_version("diann", "fragment", "2.3.0")


def test_resolve_without_version_accepts_equivalent_documents(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    path = find_rule_for_version("diann", "ion", "1.9.2").path
    monkeypatch.setattr(loader, "document_paths_for_software", lambda _software: (path, path))

    locator = loader.resolve_rule_locator_without_version("diann", "ion")

    assert isinstance(locator, RuleLocator)
    assert locator.path == path
