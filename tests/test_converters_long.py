"""Tests for converters/long.py."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from anndata_proteomics.converters.long import convert_long
from anndata_proteomics.rules.schema import ParseRule


def _build_long_rule(*, score_required: bool = False) -> ParseRule:
    return ParseRule.model_validate(
        {
            "schema_version": "0.1",
            "file_version": "1",
            "software_name": "Synthetic",
            "software_version": "1.0",
            "input_shape": "long",
            "quantification_level": "ion",
            "axis": {
                "obs_keys": ["Run"],
                "var_keys": ["Sequence", "Charge"],
                "x_layer": "Intensity",
                "duplicates": {"mode": "keep_first"},
            },
            "columns": {
                "obs": {"select": {"Run": "Run", "Condition": "Condition"}},
                "var": {"select": {"Sequence": "Sequence", "Charge": "Charge"}},
            },
            "layers": [
                {"name": "Intensity", "source": "Intensity"},
                {"name": "Score", "source": "Score", "required": score_required},
            ],
        }
    )


def _df_without_score() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Run": ["S1", "S1", "S2", "S2"],
            "Condition": ["A", "A", "B", "B"],
            "Sequence": ["PEP1", "PEP2", "PEP1", "PEP2"],
            "Charge": [2, 2, 2, 2],
            "Intensity": [10.0, 20.0, 30.0, 40.0],
        }
    )


def test_convert_long_skips_absent_optional_layer() -> None:
    # Layers are optional by default -> absent source is skipped, conversion still succeeds.
    pieces = convert_long(_df_without_score(), _build_long_rule())
    assert pieces.X.shape == (2, 2)
    assert "Intensity" in pieces.layers
    assert "Score" not in pieces.layers


def test_convert_long_raises_on_absent_required_layer() -> None:
    with pytest.raises(KeyError, match="Score"):
        convert_long(_df_without_score(), _build_long_rule(score_required=True))


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
            "software_version": "1.0",
            "input_shape": "long",
            "quantification_level": "ion",
            "axis": {
                "obs_keys": ["Run"],
                "var_keys": ["Sequence", "Charge"],
                "x_layer": "Intensity",
                "duplicates": {"mode": "keep_first"},
            },
            "columns": {
                "obs": {"select": {"Run": "Run", "Condition": "Condition"}},
                "var": {"select": {"Sequence": "Sequence", "Charge": "Charge"}},
            },
            "layers": [
                {"name": "Intensity", "source": "Intensity"},
                {"name": "Score", "source": "Score"},
            ],
        }
    )
    pieces = convert_long(df, rule)
    assert pieces.X.shape == (2, 2)  # 2 samples × 2 features
    assert list(pieces.obs.index) == ["S1", "S2"]
    assert list(pieces.obs.columns) == ["Run", "Condition"]
    assert pieces.var.shape[0] == 2
    assert "Score" in pieces.layers


def test_convert_long_replaces_declared_numeric_missing_values() -> None:
    df = _df_without_score()
    df.loc[1, "Intensity"] = 0.0
    document = _build_long_rule().model_dump(by_alias=True)
    document["layers"][0]["missing_values"] = [0]

    intensity = convert_long(df, ParseRule.model_validate(document)).layers["Intensity"]
    assert intensity[0, 0] == 10.0
    assert np.isnan(intensity[0, 1])
    assert intensity[1, 0] == 30.0
    assert intensity[1, 1] == 40.0


def test_convert_long_duplicate_modes_are_honored() -> None:
    df = pd.DataFrame(
        {
            "Run": ["S1", "S1"],
            "Condition": ["A", "A"],
            "Sequence": ["PEP1", "PEP1"],
            "Charge": [2, 2],
            "Intensity": [10.0, 20.0],
        }
    )
    document = _build_long_rule().model_dump(by_alias=True)

    document["axis"]["duplicates"]["mode"] = "error"
    with pytest.raises(ValueError, match="duplicate observation-feature keys"):
        convert_long(df, ParseRule.model_validate(document))

    document["axis"]["duplicates"]["mode"] = "keep_first"
    assert convert_long(df, ParseRule.model_validate(document)).X[0, 0] == 10.0

    document["axis"]["duplicates"]["mode"] = "aggregate"
    assert convert_long(df, ParseRule.model_validate(document)).X[0, 0] == 30.0
