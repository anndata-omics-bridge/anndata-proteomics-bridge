"""Tests for rules/loader.py."""

from __future__ import annotations

import tomllib
from pathlib import Path

import pytest
from pydantic import ValidationError

from anndata_proteomics.rules.loader import (
    _merge_rule_dicts,
    load_packaged_rule,
    load_rule,
)
from anndata_proteomics.rules.registry import find_rule


# --- merge engine (Finding 1: convention-based base/leaf inheritance) ---------


def test_merge_scalar_child_wins() -> None:
    assert _merge_rule_dicts({"a": 1, "b": 2}, {"b": 3}) == {"a": 1, "b": 3}


def test_merge_tables_deep_merge_child_wins() -> None:
    base = {"axis": {"obs_keys": ["Run"], "duplicates": {"mode": "error"}}}
    leaf = {"axis": {"var_keys": ["Ion"], "x_layer": "X"}}
    assert _merge_rule_dicts(base, leaf) == {
        "axis": {
            "obs_keys": ["Run"],
            "duplicates": {"mode": "error"},
            "var_keys": ["Ion"],
            "x_layer": "X",
        }
    }


def test_merge_nested_scalar_override() -> None:
    # protein leaves override the base's axis.duplicates.mode ("error" -> "keep_first").
    base = {"axis": {"obs_keys": ["Run"], "duplicates": {"mode": "error"}}}
    leaf = {"axis": {"var_keys": ["PG"], "x_layer": "X", "duplicates": {"mode": "keep_first"}}}
    merged = _merge_rule_dicts(base, leaf)
    assert merged["axis"]["duplicates"]["mode"] == "keep_first"
    assert merged["axis"]["obs_keys"] == ["Run"]
    assert merged["axis"]["var_keys"] == ["PG"]


def test_merge_arrays_of_tables_append_base_first() -> None:
    base = {"compute": [{"name": "peptidoform"}]}
    leaf = {"compute": [{"name": "ion"}]}
    # base entries come first so dependency order (peptidoform before ion) is preserved.
    assert _merge_rule_dicts(base, leaf)["compute"] == [{"name": "peptidoform"}, {"name": "ion"}]


def test_merge_scalar_lists_replaced_not_appended() -> None:
    # var_keys is a scalar list — the leaf replaces the base, it is not concatenated.
    assert _merge_rule_dicts({"var_keys": ["A"]}, {"var_keys": ["B"]}) == {"var_keys": ["B"]}


def test_merge_subobject_inherited_whole_when_leaf_omits_it() -> None:
    base = {"modifications": {"parser": "token_regex", "source_column": "Mod"}}
    assert _merge_rule_dicts(base, {"quantification_level": "protein"}) == {
        "modifications": {"parser": "token_regex", "source_column": "Mod"},
        "quantification_level": "protein",
    }


_BASE_TOML = """
schema_version = "0.1"
software_name = "MyVendor"
input_shape = "long"

[axis]
obs_keys = ["Run"]

[axis.duplicates]
mode = "error"

[columns.obs.select]
Run = "Run"

[columns.var.select]
Modified_Sequence = "Modified.Sequence"

[[columns.var.compute]]
name = "ProForma_peptidoform"
from = ["Modified_Sequence"]
how = "proforma_sequence"

[modifications]
source_column = "Modified.Sequence"
parser = "token_regex"
token_pattern = "\\\\(([^()]*)\\\\)"

[[modifications.map]]
token = "UniMod:35"
accession = "UNIMOD:35"
"""

_LEAF_ION_TOML = """
file_version = "1"
software_version = "1.0"
quantification_level = "ion"

[axis]
var_keys = ["ProForma_ion"]
x_layer = "Intensity"

[columns.var.select]
Precursor_Charge = "Precursor.Charge"

[[columns.var.compute]]
name = "ProForma_ion"
from = ["ProForma_peptidoform", "Precursor_Charge"]
how = "proforma_ion"

[[layers]]
name = "Intensity"
source = "Intensity"
"""


def _make_vendor(tmp_path: Path, version_folder: str | None = None) -> Path:
    """Write a base + ion leaf under tmp_path/myvendor/, return the leaf path."""
    vendor = tmp_path / "myvendor"
    vendor.mkdir()
    (vendor / "myvendor.toml").write_text(_BASE_TOML)
    leaf_dir = vendor / version_folder if version_folder else vendor
    leaf_dir.mkdir(exist_ok=True)
    leaf = leaf_dir / "parse_myvendor_ion.toml"
    leaf.write_text(_LEAF_ION_TOML)
    return leaf


def test_load_rule_merges_vendor_base(tmp_path: Path) -> None:
    rule = load_rule(_make_vendor(tmp_path))
    # scalars: base + leaf
    assert rule.software_name == "MyVendor"
    assert rule.quantification_level == "ion"
    # axis deep-merge: obs_keys from base, var_keys/x_layer from leaf
    assert rule.axis.obs_keys == ["Run"]
    assert rule.axis.var_keys == ["ProForma_ion"]
    # var.select deep-merge: base + leaf columns
    assert "Modified_Sequence" in rule.columns.var.select
    assert "Precursor_Charge" in rule.columns.var.select
    # compute append, base first (dependency order preserved)
    assert [c.name for c in rule.columns.var.compute] == ["ProForma_peptidoform", "ProForma_ion"]
    # sub-object inherited from base
    assert rule.modifications is not None
    assert rule.modifications.source_column == "Modified.Sequence"


def test_load_rule_merges_base_across_version_folder(tmp_path: Path) -> None:
    rule = load_rule(_make_vendor(tmp_path, version_folder="v1"))
    assert rule.software_name == "MyVendor"
    assert rule.modifications is not None


def test_load_rule_no_base_loads_standalone(tmp_path: Path) -> None:
    # A leaf with no sibling <vendor>.toml must still load on its own (backward compatible).
    rule = load_rule(find_rule("wombat", "peptidoform"))
    assert rule.software_name == "WOMBAT"


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
