"""Focused validation and storage edge paths for rules and summaries."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from anndata_proteomics.params.registry import get_parser
from anndata_proteomics.readers import summary
from anndata_proteomics.rules import loader
from anndata_proteomics.rules.schema import (
    ParseRule,
    ParseRuleDocument,
)


def _rule_document() -> dict[str, Any]:
    return {
        "schema_version": "0.1",
        "file_version": "test",
        "software_name": "Synthetic",
        "software_version": ".*",
        "input_shape": "long",
        "quantification_level": "ion",
        "axis": {
            "obs_keys": ["Run"],
            "var_keys": ["Feature"],
            "x_layer": "Intensity",
        },
        "columns": {
            "obs": {"select": {"Run": "Run"}},
            "var": {
                "select": {
                    "Feature": "Feature",
                    "Other": "Other",
                }
            },
        },
        "layers": [{"name": "Intensity", "source": "Intensity"}],
    }


def _source_document(*, layer_source: str = "Intensity") -> dict[str, Any]:
    rule = _rule_document()
    base = {
        key: value
        for key, value in rule.items()
        if key
        not in {
            "schema_version",
            "file_version",
            "software_name",
            "software_version",
            "quantification_level",
        }
    }
    base["layers"][0]["source"] = layer_source
    return {
        "schema_version": "0.1",
        "file_version": "test",
        "software_name": "Synthetic",
        "software_version": ".*",
        "base": base,
        "levels": {"ion": {}},
    }


def test_rule_schema_rejects_remaining_invalid_contracts() -> None:
    document = _rule_document()
    document["input_shape"] = "wide"
    document["columns"]["obs"] = {"select": {"Run": "<sample>"}}
    document["layers"][0]["source"] = "["
    with pytest.raises(ValueError, match="valid regex"):
        ParseRule.model_validate(document)

    for axis, message in [
        ({"obs_keys": ["Missing"]}, "axis.obs_keys"),
        ({"var_keys": ["Missing"]}, "axis.var_keys"),
    ]:
        document = _rule_document()
        document["axis"].update(axis)
        with pytest.raises(ValueError, match=message):
            ParseRule.model_validate(document)

    document = _rule_document()
    document["fragments"] = {
        "label_strategy": "column",
        "label_column": "Fragment.Info",
        "label_output": "Fragment",
        "value_columns": ["Fragment.Intensity"],
        "delimiter": ";",
    }
    with pytest.raises(ValueError, match="only valid"):
        ParseRule.model_validate(document)

    document = _rule_document()
    document["columns"]["var"]["compute"] = [
        {
            "name": "Combined",
            "from": ["Feature", "Missing"],
            "how": "coalesce",
        }
    ]
    with pytest.raises(ValueError, match="undeclared"):
        ParseRule.model_validate(document)

    document = _rule_document()
    document["modifications"] = {
        "source_column": "Modified",
        "parser": "token_regex",
        "token_pattern": r"\(([^)]+)\)",
        "map": [{"token": "ox", "accession": "UNIMOD:35"}],
    }
    document["columns"]["var"]["compute"] = [
        {
            "name": "ProForma_peptidoform",
            "from": ["Feature", "Other"],
            "how": "proforma_sequence",
        }
    ]
    with pytest.raises(ValueError, match="exactly one"):
        ParseRule.model_validate(document)

    document = _rule_document()
    document["quantification_level"] = "protein"
    document["columns"]["var"]["compute"] = [
        {
            "name": "ProForma_ion",
            "from": ["Feature", "Other"],
            "how": "proforma_ion",
        }
    ]
    with pytest.raises(ValueError, match="only for ion or fragment"):
        ParseRule.model_validate(document)

    document = _rule_document()
    document["columns"]["var"]["compute"] = [
        {
            "name": "ProForma_fragment",
            "from": ["Feature", "Other"],
            "how": "proforma_fragment",
        }
    ]
    with pytest.raises(ValueError, match="only for fragment"):
        ParseRule.model_validate(document)

    document = _rule_document()
    document["quantification_level"] = "fragment"
    document["columns"]["var"]["compute"] = [
        {
            "name": "ProForma_fragment",
            "from": ["Feature"],
            "how": "proforma_fragment",
        }
    ]
    with pytest.raises(ValueError, match="exactly two"):
        ParseRule.model_validate(document)

    document = _rule_document()
    document["quantification_level"] = "fragment"
    document["columns"]["var"]["compute"] = [
        {
            "name": "ProForma_fragment",
            "from": ["Feature", "Other"],
            "how": "proforma_fragment",
        }
    ]
    with pytest.raises(ValueError, match="axis.var_keys"):
        ParseRule.model_validate(document)

    document = _rule_document()
    document["columns"]["obs"]["compute"] = [
        {
            "name": "Combined",
            "from": ["Run", "Other"],
            "how": "coalesce",
        }
    ]
    with pytest.raises(ValueError, match="only for columns.var"):
        ParseRule.model_validate(document)


def test_rule_schema_accepts_a_valid_fragment_compute() -> None:
    document = _rule_document()
    document["quantification_level"] = "fragment"
    document["axis"]["var_keys"] = ["ProForma_fragment"]
    document["columns"]["var"]["compute"] = [
        {
            "name": "ProForma_fragment",
            "from": ["Feature", "Other"],
            "how": "proforma_fragment",
        }
    ]
    assert ParseRule.model_validate(document).quantification_level == "fragment"


def test_rule_document_missing_level_and_loader_paths(
    tmp_path: Path,
) -> None:
    parsed = ParseRuleDocument.model_validate(_source_document())
    with pytest.raises(KeyError):
        parsed.effective_rule("protein")

    with pytest.raises(loader.RuleDocumentError, match="root must be"):
        loader.parse_rule_source("[]")
    missing = tmp_path / "missing.json"
    with pytest.raises(FileNotFoundError):
        loader.read_rule_document(missing)
    with pytest.raises(FileNotFoundError):
        loader.load_rule_document(missing)

    path = tmp_path / "rules.json"
    path.write_text(json.dumps(_source_document()), encoding="utf-8")
    assert loader.read_rule_document(path)["software_name"] == "Synthetic"
    assert set(loader.load_rules(path)) == {"ion"}
    with pytest.raises(loader.RuleDocumentError, match="has no level"):
        loader.load_rule(path, "protein")
    with pytest.raises(ValueError, match="invalid software_version regex"):
        loader.software_version_matches("[", "1")


def test_parameter_registry_rejects_unknown_software() -> None:
    with pytest.raises(KeyError, match="no parameter parser registered"):
        get_parser("unknown")


def test_rule_locator_rejects_ambiguous_non_equivalent_documents(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    paths: list[Path] = []
    for index, source in enumerate(("Intensity", "Abundance")):
        path = tmp_path / f"rules-{index}.json"
        path.write_text(
            json.dumps(_source_document(layer_source=source)),
            encoding="utf-8",
        )
        paths.append(path)
    monkeypatch.setattr(
        loader,
        "document_paths_for_software",
        lambda _software: tuple(paths),
    )
    assert loader.resolve_rule_locator("synthetic", "ion", None) is None


def test_summary_json_compatibility_helper() -> None:
    converted = summary._to_json_compatible(
        {
            "array": np.asarray([np.int64(1)]),
            "tuple": (np.float64(2.0),),
            "bytes": b"text",
        }
    )
    assert converted == {
        "array": [1],
        "tuple": [2.0],
        "bytes": "text",
    }
