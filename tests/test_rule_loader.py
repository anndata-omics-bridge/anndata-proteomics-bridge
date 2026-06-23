"""Tests for rules/loader.py."""

from __future__ import annotations

import tomllib
from pathlib import Path

import pytest
from pydantic import ValidationError

from anndata_proteomics.rules.loader import load_packaged_rule, load_rule
from anndata_proteomics.rules.registry import find_rule


def test_load_rule_happy() -> None:
    rule = load_rule(find_rule("diann", "ion"))
    assert rule.software_name == "DIA-NN"
    assert rule.quantification_level == "ion"
    assert rule.input_shape == "long"


def test_load_rule_missing_path(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        load_rule(tmp_path / "nope.toml")


def test_load_rule_malformed_toml(tmp_path: Path) -> None:
    bad = tmp_path / "bad.toml"
    bad.write_text("this is not = valid toml [[[")
    with pytest.raises(tomllib.TOMLDecodeError):
        load_rule(bad)


def test_load_rule_pydantic_invalid_attaches_path(tmp_path: Path) -> None:
    # Valid TOML but missing required quantification_level field.
    bad = tmp_path / "missing_field.toml"
    bad.write_text(
        """
schema_version = "0.1"
file_version = "1"
software_name = "Fake"
input_shape = "long"

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
source = "Foo"
"""
    )
    with pytest.raises(ValidationError) as exc_info:
        load_rule(bad)
    notes = getattr(exc_info.value, "__notes__", [])
    assert any(str(bad) in n for n in notes), f"path not attached as note: {notes}"


def test_load_packaged_rule_matches_direct_load() -> None:
    direct = load_rule(find_rule("wombat", "peptidoform"))
    via_packaged = load_packaged_rule("wombat", "peptidoform")
    assert direct.model_dump() == via_packaged.model_dump()
