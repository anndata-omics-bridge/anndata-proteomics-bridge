"""Tests for converters/wide.py."""

from __future__ import annotations

import pandas as pd

from anndata_proteomics.converters.wide import convert_wide
from anndata_proteomics.rules.schema import ParseRule


def test_convert_wide_happy_path() -> None:
    df = pd.DataFrame(
        {
            "Modified Sequence": ["PEP1", "PEP2", "PEP3"],
            "Charge": [2, 2, 3],
            "S1 Intensity": [10.0, 20.0, 30.0],
            "S2 Intensity": [11.0, 21.0, 31.0],
            "S1 Spectral Count": [1, 2, 3],
            "S2 Spectral Count": [4, 5, 6],
        }
    )
    rule = ParseRule.model_validate(
        {
            "schema_version": "0.1",
            "file_version": "1",
            "software_name": "Synthetic",
            "software_version": "1.0",
            "input_shape": "wide",
            "quantification_level": "ion",
            "axis": {
                "obs_keys": ["sample"],
                "var_keys": ["Modified Sequence", "Charge"],
                "x_layer": "Intensity",
                "duplicates": {"mode": "keep_first"},
            },
            "columns": {
                "obs": {"select": {"sample": "<sample>"}},
                "var": {
                    "select": {
                        "Modified Sequence": "Modified Sequence",
                        "Charge": "Charge",
                    }
                },
            },
            "layers": [
                {
                    "name": "Intensity",
                    "source": "^(?P<sample>S\\d+) Intensity$",
                },
                {
                    "name": "Spectral_Count",
                    "source": "^(?P<sample>S\\d+) Spectral Count$",
                },
            ],
        }
    )
    pieces = convert_wide(df, rule)
    assert pieces.X.shape == (2, 3)  # 2 samples × 3 features
    assert list(pieces.obs.index) == ["S1", "S2"]
    assert pieces.var.shape == (3, 2)
    assert "Spectral_Count" in pieces.layers
    assert pieces.layers["Spectral_Count"].shape == (2, 3)


def test_convert_wide_factor_layer() -> None:
    df = pd.DataFrame(
        {
            "Modified Sequence": ["PEP1", "PEP2"],
            "Charge": [2, 2],
            "S1 Intensity": [10.0, 20.0],
            "S2 Intensity": [11.0, 21.0],
            "S1 Type": ["MS/MS", "MBR"],
            "S2 Type": ["MS/MS", "unmatched"],
        }
    )
    rule = ParseRule.model_validate(
        {
            "schema_version": "0.1",
            "file_version": "1",
            "software_name": "Synthetic",
            "software_version": "1.0",
            "input_shape": "wide",
            "quantification_level": "ion",
            "axis": {
                "obs_keys": ["sample"],
                "var_keys": ["Modified Sequence", "Charge"],
                "x_layer": "Intensity",
                "duplicates": {"mode": "keep_first"},
            },
            "columns": {
                "obs": {"select": {"sample": "<sample>"}},
                "var": {
                    "select": {
                        "Modified Sequence": "Modified Sequence",
                        "Charge": "Charge",
                    }
                },
            },
            "layers": [
                {"name": "Intensity", "source": "^(?P<sample>S\\d+) Intensity$"},
                {
                    "name": "Type",
                    "source": "^(?P<sample>S\\d+) Type$",
                    "encoding_mode": "factor",
                    "categories": {"unmatched": 0, "MS/MS": 1, "MBR": 2},
                },
            ],
        }
    )
    pieces = convert_wide(df, rule)
    type_layer = pieces.layers["Type"]
    # S1: ['MS/MS', 'MBR'] → [1, 2];  S2: ['MS/MS', 'unmatched'] → [1, 0]
    assert type_layer[0, 0] == 1
    assert type_layer[0, 1] == 2
    assert type_layer[1, 0] == 1
    assert type_layer[1, 1] == 0
