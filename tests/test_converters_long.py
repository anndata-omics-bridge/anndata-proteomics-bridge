"""Tests for converters/long.py."""

from __future__ import annotations

import pandas as pd

from anndata_proteomics.converters.long import convert_long
from anndata_proteomics.rules.schema import ParseRule


def _build_long_rule() -> ParseRule:
    return ParseRule.model_validate(
        {
            "schema_version": "0.1",
            "file_version": "1",
            "software_name": "Synthetic",
            "input_shape": "long",
            "quantification_level": "ion",
            "axis": {
                "obs_keys": ["Run"],
                "var_keys": ["Sequence", "Charge"],
                "x_layer": "Intensity",
            },
            "columns": {
                "obs": {"Run": "Run", "Condition": "Condition"},
                "var": {"Sequence": "Sequence", "Charge": "Charge"},
            },
            "layers": [
                {"name": "Intensity", "source_column": "Intensity"},
                {"name": "Score", "source_column": "Score"},
            ],
            "duplicates": {"mode": "first"},
        }
    )


def test_convert_long_happy_path() -> None:
    df = pd.DataFrame(
        {
            "Run": ["S1", "S1", "S2", "S2"],
            "Condition": ["A", "A", "B", "B"],
            "Sequence": ["PEP1", "PEP2", "PEP1", "PEP2"],
            "Charge": [2, 2, 2, 2],
            "Intensity": [10.0, 20.0, 30.0, 40.0],
            "Score": [0.9, 0.8, 0.7, 0.6],
        }
    )
    rule = ParseRule.model_validate(
        {
            "schema_version": "0.1",
            "file_version": "1",
            "software_name": "Synthetic",
            "input_shape": "long",
            "quantification_level": "ion",
            "axis": {
                "obs_keys": ["Run"],
                "var_keys": ["Sequence", "Charge"],
                "x_layer": "Intensity",
            },
            "columns": {
                "obs": {"Run": "Run", "Condition": "Condition"},
                "var": {"Sequence": "Sequence", "Charge": "Charge"},
            },
            "layers": [
                {"name": "Intensity", "source_column": "Intensity"},
                {"name": "Score", "source_column": "Score"},
            ],
            "duplicates": {"mode": "keep_first"},
        }
    )
    pieces = convert_long(df, rule)
    assert pieces.X.shape == (2, 2)  # 2 samples × 2 features
    assert list(pieces.obs.index) == ["S1", "S2"]
    assert list(pieces.obs.columns) == ["Run", "Condition"]
    assert pieces.var.shape[0] == 2
    assert "Score" in pieces.layers
