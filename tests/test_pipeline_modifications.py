"""End-to-end pipeline tests for [modifications] + params_path integration."""

from __future__ import annotations

import json
import tomllib
from pathlib import Path

import pandas as pd
import pytest

from anndata_proteomics.converters.assemble import convert
from anndata_proteomics.rules.schema import ParseRule


RULE_TOML = """
schema_version = "0.1"
file_version = "1"
software_name = "Sage"
input_shape = "long"
quantification_level = "ion"

[axis]
obs_keys = ["Run"]
var_keys = ["ProForma"]
x_layer = "Intensity"

[axis.duplicates]
mode = "error"

[columns.obs.select]
Run = "Run"

[columns.var.select]
Vendor_Sequence = "Modified.Sequence"
Precursor_Charge = "Precursor.Charge"

[[columns.var.compute]]
name = "Peptidoform"
from = ["Vendor_Sequence"]
how = "proforma_sequence"

[[columns.var.compute]]
name = "Stripped"
from = ["Vendor_Sequence"]
how = "stripped_sequence"

[[columns.var.compute]]
name = "ProForma"
from = ["Peptidoform", "Precursor_Charge"]
how = "proforma_ion"

[[layers]]
name = "Intensity"
source_column = "Intensity"

[modifications]
source_column = "Modified.Sequence"
parser = "token_regex"
token_pattern = "\\\\[([^]]+)\\\\]"
token_position = "after_residue"
unknown_policy = "preserve"
output_column = "proforma_sequence"

[[modifications.map]]
token = "15.9949"
accession = "UNIMOD:35"

[[modifications.map]]
token = "57.0215"
accession = "UNIMOD:4"
"""


def _make_rule() -> ParseRule:
    return ParseRule(**tomllib.loads(RULE_TOML))


def _make_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Run": ["run1", "run1", "run2", "run2"],
            "Modified.Sequence": [
                "PEPM[15.9949]TIDE",
                "PEPC[57.0215]TIDE",
                "PEPM[15.9949]TIDE",
                "PEPC[57.0215]TIDE",
            ],
            "Precursor.Charge": [2, 2, 2, 2],
            "Intensity": [100.0, 200.0, 150.0, 250.0],
        }
    )


def test_convert_adds_proforma_column_to_var():
    adata = convert(_make_df(), _make_rule())
    assert "ProForma" in adata.var.columns
    proforma_values = sorted(adata.var["ProForma"].tolist())
    assert "PEPM[UNIMOD:35]TIDE/2" in proforma_values
    assert "PEPC[UNIMOD:4]TIDE/2" in proforma_values
    peptidoform_values = sorted(adata.var["Peptidoform"].tolist())
    assert "PEPM[UNIMOD:35]TIDE" in peptidoform_values


def test_convert_var_indexed_by_proforma():
    adata = convert(_make_df(), _make_rule())
    assert sorted(adata.var_names) == [
        "PEPC[UNIMOD:4]TIDE/2",
        "PEPM[UNIMOD:35]TIDE/2",
    ]


def test_convert_normalizes_float_charge_in_proforma_ion():
    df = _make_df()
    df["Precursor.Charge"] = [2.0, 2.0, 2.0, 2.0]
    adata = convert(df, _make_rule())
    assert sorted(adata.var_names) == [
        "PEPC[UNIMOD:4]TIDE/2",
        "PEPM[UNIMOD:35]TIDE/2",
    ]


def test_convert_rejects_non_positive_charge_for_proforma_ion():
    df = _make_df()
    df.loc[0, "Precursor.Charge"] = 0
    with pytest.raises(ValueError, match="positive"):
        convert(df, _make_rule())


def test_convert_rejects_non_integral_charge_for_proforma_ion():
    df = _make_df()
    df["Precursor.Charge"] = df["Precursor.Charge"].astype(float)
    df.loc[0, "Precursor.Charge"] = 2.5
    with pytest.raises(ValueError, match="integer"):
        convert(df, _make_rule())


def test_convert_with_params_path_attaches_search_parameters(tmp_path):
    proteobench_params = Path(
        "/Users/wolski/projects/anndata_bridge/ProteoBench/test/params/sage_parameterfile.json"
    )
    if not proteobench_params.exists():
        pytest.skip("ProteoBench fixture missing")
    adata = convert(_make_df(), _make_rule(), params_path=proteobench_params)
    uns = adata.uns["anndata_proteomics"]
    assert "search_parameters" in uns
    parsed = json.loads(uns["search_parameters"])
    assert parsed["software_name"] == "Sage"
    assert parsed["software_version"] == "0.14.6"
    assert uns["search_parameters_path"] == str(proteobench_params)


def test_convert_with_params_path_for_unknown_software_keeps_path_only(tmp_path):
    rule_toml = RULE_TOML.replace('software_name = "Sage"', 'software_name = "UnknownTool"')
    rule = ParseRule(**tomllib.loads(rule_toml))
    fake = tmp_path / "fake_params.txt"
    fake.write_text("dummy")
    adata = convert(_make_df(), rule, params_path=fake)
    uns = adata.uns["anndata_proteomics"]
    assert "search_parameters" not in uns
    assert uns["search_parameters_path"] == str(fake)
