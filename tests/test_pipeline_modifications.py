"""End-to-end pipeline tests for [modifications] + params_path integration."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from conversion_support import convert_to_anndata

from anndata_proteomics.adapters.anndata import conversion as conversion_adapter
from anndata_proteomics.converters import pipeline as conversion_pipeline
from anndata_proteomics.rules.schema import ParseRule

RULE: dict[str, Any] = {
    "schema_version": "0.1",
    "file_version": "1",
    "software_name": "Sage",
    "software_version": "1.0",
    "input_shape": "long",
    "quantification_level": "ion",
    "axis": {
        "obs_keys": ["Run"],
        "var_keys": ["ProForma_ion"],
        "x_layer": "Intensity",
        "duplicates": {"mode": "error"},
    },
    "columns": {
        "obs": {"select": {"Run": "Run"}},
        "var": {
            "select": {
                "Vendor_Sequence": "Modified.Sequence",
                "Precursor_Charge": "Precursor.Charge",
            },
            "compute": [
                {
                    "name": "ProForma_peptidoform",
                    "from": ["Vendor_Sequence"],
                    "how": "proforma_sequence",
                },
                {
                    "name": "ProForma_peptide",
                    "from": ["Vendor_Sequence"],
                    "how": "stripped_sequence",
                },
                {
                    "name": "ProForma_ion",
                    "from": ["ProForma_peptidoform", "Precursor_Charge"],
                    "how": "proforma_ion",
                },
            ],
        },
    },
    "layers": [{"name": "Intensity", "source": "Intensity"}],
    "modifications": {
        "source_column": "Modified.Sequence",
        "parser": "token_regex",
        "token_pattern": r"\[([^]]+)\]",
        "token_position": "after_residue",
        "unknown_policy": "preserve",
        "output_column": "proforma_sequence",
        "map": [
            {"token": "15.9949", "accession": "UNIMOD:35"},
            {"token": "57.0215", "accession": "UNIMOD:4"},
        ],
    },
}


def _make_rule() -> ParseRule:
    return ParseRule.model_validate(RULE)


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
    adata = convert_to_anndata(_make_df(), _make_rule())
    assert "ProForma_ion" in adata.var.columns
    proforma_values = sorted(adata.var["ProForma_ion"].tolist())
    assert "PEPM[UNIMOD:35]TIDE/2" in proforma_values
    assert "PEPC[UNIMOD:4]TIDE/2" in proforma_values
    peptidoform_values = sorted(adata.var["ProForma_peptidoform"].tolist())
    assert "PEPM[UNIMOD:35]TIDE" in peptidoform_values


def test_convert_var_indexed_by_proforma():
    adata = convert_to_anndata(_make_df(), _make_rule())
    assert sorted(adata.var_names) == [
        "PEPC[UNIMOD:4]TIDE/2",
        "PEPM[UNIMOD:35]TIDE/2",
    ]


def test_convert_normalizes_float_charge_in_proforma_ion():
    df = _make_df()
    df["Precursor.Charge"] = [2.0, 2.0, 2.0, 2.0]
    adata = convert_to_anndata(df, _make_rule())
    assert sorted(adata.var_names) == [
        "PEPC[UNIMOD:4]TIDE/2",
        "PEPM[UNIMOD:35]TIDE/2",
    ]


def test_convert_rejects_non_positive_charge_for_proforma_ion():
    df = _make_df()
    df.loc[0, "Precursor.Charge"] = 0
    with pytest.raises(ValueError, match="positive"):
        convert_to_anndata(df, _make_rule())


def test_convert_rejects_non_integral_charge_for_proforma_ion():
    df = _make_df()
    df["Precursor.Charge"] = df["Precursor.Charge"].astype(float)
    df.loc[0, "Precursor.Charge"] = 2.5
    with pytest.raises(ValueError, match="integer"):
        convert_to_anndata(df, _make_rule())


def test_convert_with_params_path_attaches_search_parameters(tmp_path: Path) -> None:
    proteobench_params = Path(__file__).resolve().parent / "params" / "sage_parameterfile.json"
    if not proteobench_params.exists():
        pytest.skip("ProteoBench fixture missing")
    adata = convert_to_anndata(_make_df(), _make_rule())
    resolution = conversion_pipeline.resolve_parameters(proteobench_params, "sage")
    conversion_adapter.write_parameter_resolution(
        adata,
        resolution,
    )
    uns = adata.uns["anndata_proteomics"]
    assert "search_parameters" in uns
    parsed = json.loads(uns["search_parameters"])
    assert parsed["software_name"] == "Sage"
    assert parsed["software_version"] == "0.14.6"
    assert uns["search_parameters_path"] == str(proteobench_params)
