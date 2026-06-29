"""Tests for the canonical-modification registry."""

from __future__ import annotations

import tomllib

import pytest
from pydantic import ValidationError

from anndata_proteomics.modifications.pipeline import _to_runtime_rule
from anndata_proteomics.modifications.unimod_registry import (
    UnimodRegistry,
    load_registry,
    resolve,
)
from anndata_proteomics.rules.schema import ParseRule


BASE = """
schema_version = "0.1"
file_version = "1"
software_name = "T"
software_version = "1.0"
input_shape = "long"
quantification_level = "ion"

[axis]
obs_keys = ["Run"]
var_keys = ["ProForma_peptidoform", "Charge"]
x_layer = "Intensity"

[axis.duplicates]
mode = "error"

[columns.obs.select]
Run = "Run"

[columns.var.select]
Modified_Sequence = "Modified Sequence"
Charge = "Charge"

[[columns.var.compute]]
name = "ProForma_peptidoform"
from = ["Modified_Sequence"]
how = "proforma_sequence"

[[layers]]
name = "Intensity"
source = "Intensity"
"""


def test_registry_loads_with_required_accessions():
    registry = load_registry()
    for accession in (
        "UNIMOD:1",
        "UNIMOD:4",
        "UNIMOD:21",
        "UNIMOD:35",
        "UNIMOD:27",
        "UNIMOD:28",
        "UNIMOD:121",
    ):
        assert accession in registry, f"{accession} should be in the registry"


def test_resolve_returns_canonical_record():
    entry = resolve("UNIMOD:35")
    assert entry.name == "Oxidation"
    assert entry.target == ["M"]
    assert entry.position == "Anywhere"
    assert entry.mass_delta == pytest.approx(15.9949)


def test_resolve_phospho_targets_s_t_y():
    entry = resolve("UNIMOD:21")
    assert entry.name == "Phospho"
    assert entry.target == ["S", "T", "Y"]


def test_resolve_glygly_on_lysine():
    entry = resolve("UNIMOD:121")
    assert entry.target == ["K"]
    assert entry.mass_delta == pytest.approx(114.04293)


def test_resolve_unknown_accession_raises():
    with pytest.raises(KeyError, match="UNIMOD:99999"):
        resolve("UNIMOD:99999")


def test_registry_rejects_duplicate_accession(tmp_path, monkeypatch):
    # Point load_registry at a temp TOML with a repeated accession and confirm it raises.
    from anndata_proteomics.modifications import unimod_registry as reg

    dup = tmp_path / "dup_registry.toml"
    dup.write_text(
        """
[[entries]]
accession = "UNIMOD:35"
name = "A"
target = ["M"]
position = "Anywhere"
mass_delta = 1.0

[[entries]]
accession = "UNIMOD:35"
name = "B"
target = ["M"]
position = "Anywhere"
mass_delta = 2.0
"""
    )
    monkeypatch.setattr(reg, "_REGISTRY_TOML", dup)
    reg.load_registry.cache_clear()
    try:
        with pytest.raises(ValueError, match="duplicate accession"):
            reg.load_registry()
    finally:
        reg.load_registry.cache_clear()  # drop the temp result so other tests see the real registry


def test_registry_entry_extras_forbidden():
    with pytest.raises(ValidationError):
        UnimodRegistry(
            entries=[
                {
                    "accession": "UNIMOD:35",
                    "name": "Oxidation",
                    "target": ["M"],
                    "position": "Anywhere",
                    "mass_delta": 15.9949,
                    "vendor_specific": "nope",
                }
            ]
        )


def test_runtime_rule_resolves_canonical_fields_from_accession():
    rule_toml = (
        BASE
        + """
[modifications]
source_column = "Modified Sequence"
parser = "token_regex"
token_pattern = "\\\\[([^\\\\]]+)\\\\]"
token_position = "after_residue"

[[modifications.map]]
token = "15.9949"
accession = "UNIMOD:35"
"""
    )
    rule = ParseRule(**tomllib.loads(rule_toml))
    runtime = _to_runtime_rule(rule.modifications)
    assert runtime.entries[0].name == "Oxidation"
    assert runtime.entries[0].target == ["M"]
    assert runtime.entries[0].position == "Anywhere"
    assert runtime.entries[0].mass_delta == pytest.approx(15.9949)


def test_runtime_rule_errors_on_unknown_accession():
    rule_toml = (
        BASE
        + """
[modifications]
source_column = "Modified Sequence"
parser = "token_regex"
token_pattern = "\\\\[([^\\\\]]+)\\\\]"

[[modifications.map]]
token = "nope"
accession = "UNIMOD:99999"
"""
    )
    rule = ParseRule(**tomllib.loads(rule_toml))
    with pytest.raises(KeyError, match="UNIMOD:99999"):
        _to_runtime_rule(rule.modifications)
