"""Tests for the optional [modifications] block on ParseRule."""

from __future__ import annotations

import tomllib

import pytest
from pydantic import ValidationError

from anndata_proteomics.rules.schema import ParseRule


BASE = """
schema_version = "0.1"
file_version = "1"
software_name = "DIA-NN"
software_version = "1.9.1"
input_shape = "long"
quantification_level = "ion"

[axis]
obs_keys = ["Run"]
var_keys = ["Modified.Sequence", "Precursor.Charge"]
x_layer = "Precursor_Normalised"

[columns.obs]
Run = "Run"

[columns.var]
Modified_Sequence = "Modified.Sequence"
Precursor_Charge = "Precursor.Charge"

[[layers]]
name = "Precursor_Normalised"
source_column = "Precursor.Normalised"

[duplicates]
mode = "error"
"""


def _parse(toml: str) -> ParseRule:
    return ParseRule(**tomllib.loads(toml))


def test_rule_without_modifications_still_validates():
    rule = _parse(BASE)
    assert rule.modifications is None


def test_rule_with_token_regex_modifications():
    extra = """
[modifications]
source_column = "Modified.Sequence"
parser = "token_regex"
token_pattern = "\\\\(([^()]*)\\\\)"
token_position = "after_residue"
unknown_policy = "preserve"

[[modifications.map]]
token = "UniMod:35"
name = "Oxidation"
accession = "UNIMOD:35"
target = "M"
position = "Anywhere"
mass_delta = 15.9949
"""
    rule = _parse(BASE + extra)
    assert rule.modifications is not None
    assert rule.modifications.parser == "token_regex"
    assert rule.modifications.map[0].accession == "UNIMOD:35"


def test_token_regex_requires_token_pattern():
    extra = """
[modifications]
source_column = "Modified.Sequence"
parser = "token_regex"

[[modifications.map]]
token = "ox"
name = "Oxidation"
"""
    with pytest.raises(ValidationError, match="token_pattern"):
        _parse(BASE + extra)


def test_token_regex_requires_map_entries():
    extra = """
[modifications]
source_column = "Modified.Sequence"
parser = "token_regex"
token_pattern = "\\\\[([^]]+)\\\\]"
"""
    with pytest.raises(ValidationError, match="map"):
        _parse(BASE + extra)


def test_already_proforma_rejects_token_pattern():
    extra = """
[modifications]
source_column = "ProForma"
parser = "already_proforma"
token_pattern = "anything"
"""
    with pytest.raises(ValidationError, match="token_pattern"):
        _parse(BASE + extra)


def test_already_proforma_rejects_map():
    extra = """
[modifications]
source_column = "ProForma"
parser = "already_proforma"

[[modifications.map]]
token = "x"
name = "X"
"""
    with pytest.raises(ValidationError, match="map"):
        _parse(BASE + extra)


def test_modifications_extra_keys_forbidden():
    extra = """
[modifications]
source_column = "Modified.Sequence"
parser = "token_regex"
token_pattern = "\\\\[([^]]+)\\\\]"
weird_field = true

[[modifications.map]]
token = "ox"
name = "Oxidation"
"""
    with pytest.raises(ValidationError):
        _parse(BASE + extra)
