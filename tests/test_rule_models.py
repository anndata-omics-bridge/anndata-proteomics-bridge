"""Tests for the pydantic ParseRule schema."""

from __future__ import annotations

import copy
from typing import Any

import pytest
from pydantic import ValidationError

from anndata_proteomics.rules.schema import ParseRule


LONG_EXAMPLE: dict[str, Any] = {
    "schema_version": "0.1",
    "file_version": "1",
    "software_name": "DIA-NN",
    "software_version": "1.9.1",
    "input_shape": "long",
    "quantification_level": "ion",
    "axis": {
        "obs_keys": ["Run"],
        "var_keys": ["ProForma_ion"],
        "x_layer": "Precursor_Normalised",
        "duplicates": {"mode": "error"},
    },
    "columns": {
        "obs": {"select": {"File_Name": "File.Name", "Run": "Run"}},
        "var": {
            "select": {
                "Modified_Sequence": "Modified.Sequence",
                "Protein_Ids": "Protein.Ids",
                "Precursor_Charge": "Precursor.Charge",
                "Genes": "Genes",
            },
            "compute": [
                {
                    "name": "ProForma_peptidoform",
                    "from": ["Modified_Sequence"],
                    "how": "proforma_sequence",
                },
                {
                    "name": "ProForma_ion",
                    "from": ["ProForma_peptidoform", "Precursor_Charge"],
                    "how": "proforma_ion",
                },
            ],
        },
    },
    "layers": [
        {"name": "Precursor_Normalised", "source": "Precursor.Normalised"},
        {"name": "Q_Value", "source": "Q.Value"},
        {"name": "RT", "source": "RT"},
        {"name": "Ms1_Area", "source": "Ms1.Area"},
    ],
    "modifications": {
        "source_column": "Modified.Sequence",
        "parser": "token_regex",
        "token_pattern": r"\(([^()]*)\)",
        "map": [{"token": "UniMod:35", "accession": "UNIMOD:35"}],
    },
}


WIDE_EXAMPLE: dict[str, Any] = {
    "schema_version": "0.1",
    "file_version": "1",
    "software_name": "FragPipe",
    "software_version": "23.0",
    "input_shape": "wide",
    "quantification_level": "ion",
    "axis": {
        "obs_keys": ["sample"],
        "var_keys": ["Modified_Sequence", "Charge"],
        "x_layer": "Intensity",
        "duplicates": {"mode": "error"},
    },
    "columns": {
        "obs": {"select": {"sample": "<sample>"}},
        "var": {
            "select": {
                "Peptide_Sequence": "Peptide Sequence",
                "Modified_Sequence": "Modified Sequence",
                "Charge": "Charge",
                "Protein_ID": "Protein ID",
                "Gene": "Gene",
            }
        },
    },
    "layers": [
        {"name": "Intensity", "source": r"^(?P<sample>.+) Intensity$"},
        {"name": "Spectral_Count", "source": r"^(?P<sample>.+) Spectral Count$"},
        {
            "name": "Match_Type",
            "source": r"^(?P<sample>.+) Match Type$",
            "encoding_mode": "factor",
            "categories": {"unmatched": 0, "MS/MS": 1, "MBR": 2},
        },
        {
            "name": "Localization",
            "source": r"^(?P<sample>.+) Localization$",
            "encoding_mode": "factor",
            "categories": {"Localized": 1, "Ambiguous": 0},
        },
    ],
    "sample_name_cleanup": {"pattern": ""},
}


def _parse(document: dict[str, Any]) -> ParseRule:
    return ParseRule.model_validate(document)


def test_long_example_validates():
    rule = _parse(LONG_EXAMPLE)
    assert rule.input_shape == "long"
    assert rule.software_name == "DIA-NN"
    assert len(rule.layers) == 4
    assert rule.axis.x_layer == "Precursor_Normalised"
    assert rule.axis.duplicates.mode == "error"


def test_wide_example_validates():
    rule = _parse(WIDE_EXAMPLE)
    assert rule.input_shape == "wide"
    match_type = next(layer for layer in rule.layers if layer.name == "Match_Type")
    assert match_type.encoding_mode == "factor"
    assert match_type.categories == {"unmatched": 0, "MS/MS": 1, "MBR": 2}


def test_long_layer_missing_source():
    bad = copy.deepcopy(LONG_EXAMPLE)
    bad["layers"][2].pop("source")
    with pytest.raises(ValidationError, match="source"):
        _parse(bad)


def test_layer_column_pattern_is_unknown_field():
    # column_pattern was removed from the model; it is now just an extra (forbidden) key.
    bad = copy.deepcopy(LONG_EXAMPLE)
    bad["layers"][2]["column_pattern"] = "^.+ RT$"
    with pytest.raises(ValidationError, match="column_pattern"):
        _parse(bad)


def test_wide_layer_missing_source():
    bad = copy.deepcopy(WIDE_EXAMPLE)
    bad["layers"][0].pop("source")
    with pytest.raises(ValidationError, match="source"):
        _parse(bad)


def test_layer_source_column_is_unknown_field():
    # source_column was removed from the layer model; it is now just an extra (forbidden) key.
    bad = copy.deepcopy(WIDE_EXAMPLE)
    bad["layers"][0]["source_column"] = "Intensity"
    with pytest.raises(ValidationError, match="source_column"):
        _parse(bad)


def test_wide_source_requires_sample_group():
    bad = copy.deepcopy(WIDE_EXAMPLE)
    bad["layers"][0]["source"] = "^.+ Intensity$"
    with pytest.raises(ValidationError, match="sample"):
        _parse(bad)


def test_factor_requires_categories():
    bad = copy.deepcopy(WIDE_EXAMPLE)
    bad["layers"][2].pop("categories")
    with pytest.raises(ValidationError, match="categories"):
        _parse(bad)


def test_x_layer_must_exist():
    bad = copy.deepcopy(LONG_EXAMPLE)
    bad["axis"]["x_layer"] = "DoesNotExist"
    with pytest.raises(ValidationError, match="x_layer"):
        _parse(bad)


def test_layer_required_defaults_to_x_layer_only():
    rule = _parse(LONG_EXAMPLE)
    by_name = {layer.name: layer for layer in rule.layers}
    assert rule.layer_required(by_name["Precursor_Normalised"]) is True  # x_layer always required
    assert rule.layer_required(by_name["Q_Value"]) is False  # optional by default


def test_invalid_duplicates_mode():
    bad = copy.deepcopy(LONG_EXAMPLE)
    bad["axis"]["duplicates"]["mode"] = "wrong"
    with pytest.raises(ValidationError):
        _parse(bad)


def test_top_level_duplicates_rejected():
    bad = {**copy.deepcopy(LONG_EXAMPLE), "duplicates": {"mode": "error"}}
    with pytest.raises(ValidationError, match="duplicates"):
        _parse(bad)


def test_unknown_top_level_key_rejected():
    bad = {**copy.deepcopy(LONG_EXAMPLE), "foo": "bar"}
    with pytest.raises(ValidationError, match="foo"):
        _parse(bad)


def test_sample_name_cleanup_rejected_for_long():
    bad = {**copy.deepcopy(LONG_EXAMPLE), "sample_name_cleanup": {"pattern": "(.+)"}}
    with pytest.raises(ValidationError, match="sample_name_cleanup"):
        _parse(bad)


def test_json_schema_export_has_expected_top_level_properties():
    schema = ParseRule.model_json_schema()
    expected = {
        "schema_version",
        "file_version",
        "software_name",
        "software_version",
        "input_shape",
        "quantification_level",
        "axis",
        "columns",
        "layers",
        "sample_name_cleanup",
        "modifications",
        "fragments",
    }
    assert set(schema["properties"]) == expected


def test_invalid_quantification_level():
    bad = copy.deepcopy(LONG_EXAMPLE)
    bad["quantification_level"] = "wrong"
    with pytest.raises(ValidationError):
        _parse(bad)


def test_fragment_level_can_be_native_row_level_without_fragments_block():
    good = copy.deepcopy(LONG_EXAMPLE)
    good["quantification_level"] = "fragment"
    rule = _parse(good)
    assert rule.quantification_level == "fragment"
    assert rule.fragments is None


def test_missing_quantification_level():
    bad = copy.deepcopy(LONG_EXAMPLE)
    bad.pop("quantification_level")
    with pytest.raises(ValidationError, match="quantification_level"):
        _parse(bad)


def test_missing_software_version():
    bad = copy.deepcopy(LONG_EXAMPLE)
    bad.pop("software_version")
    with pytest.raises(ValidationError, match="software_version"):
        _parse(bad)


def test_proforma_ion_requires_two_sources():
    bad = copy.deepcopy(LONG_EXAMPLE)
    bad["columns"]["var"]["compute"][1]["from"] = ["ProForma_peptidoform"]
    with pytest.raises(ValidationError, match="exactly two"):
        _parse(bad)


def test_proforma_ion_must_be_var_axis_key():
    bad = copy.deepcopy(LONG_EXAMPLE)
    bad["axis"]["var_keys"] = ["ProForma_peptidoform"]
    with pytest.raises(ValidationError, match="axis.var_keys"):
        _parse(bad)


def test_proforma_compute_names_are_pinned():
    bad = copy.deepcopy(LONG_EXAMPLE)
    bad["columns"]["var"]["compute"][0]["name"] = "MyPeptidoform"
    with pytest.raises(ValidationError, match="ProForma_peptidoform"):
        _parse(bad)


def test_apb_derived_columns_cannot_be_selected():
    bad = copy.deepcopy(LONG_EXAMPLE)
    bad["columns"]["var"]["select"]["Bad"] = "proforma_sequence"
    with pytest.raises(ValidationError, match="derived"):
        _parse(bad)


def test_proforma_sequence_compute_requires_modifications():
    bad = copy.deepcopy(LONG_EXAMPLE)
    bad.pop("modifications")
    with pytest.raises(ValidationError, match="modifications"):
        _parse(bad)
