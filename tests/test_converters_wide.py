"""Tests for converters/wide.py."""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
import pandas as pd
import pytest

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


def test_convert_wide_replaces_declared_numeric_missing_values() -> None:
    df = pd.DataFrame(
        {
            "Modified Sequence": ["PEP1", "PEP2"],
            "Charge": [2, 2],
            "S1 Intensity": [0.0, 20.0],
            "S2 Intensity": [11.0, 0.0],
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
                    "missing_values": [0],
                }
            ],
        }
    )

    intensity = convert_wide(df, rule).layers["Intensity"]
    assert np.isnan(intensity[0, 0])
    assert intensity[0, 1] == 20.0
    assert intensity[1, 0] == 11.0
    assert np.isnan(intensity[1, 1])


def test_convert_wide_duplicate_modes_are_honored() -> None:
    df = pd.DataFrame(
        {
            "Modified Sequence": ["PEP1", "PEP1"],
            "Charge": [2, 2],
            "S1 Intensity": [10.0, 20.0],
        }
    )
    document: dict[str, Any] = {
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
            "duplicates": {"mode": "error"},
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
        "layers": [{"name": "Intensity", "source": "^(?P<sample>S\\d+) Intensity$"}],
    }

    with pytest.raises(ValueError, match="duplicate observation-feature keys"):
        convert_wide(df, ParseRule.model_validate(document))

    document["axis"]["duplicates"]["mode"] = "keep_first"
    assert convert_wide(df, ParseRule.model_validate(document)).X[0, 0] == 10.0

    document["axis"]["duplicates"]["mode"] = "aggregate"
    assert convert_wide(df, ParseRule.model_validate(document)).X[0, 0] == 30.0


def test_x_layer_alone_defines_observation_axis(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """An optional-layer-only token must not become an observation."""
    df = pd.DataFrame(
        {
            "Feature": ["F1"],
            "S1 Intensity": [10.0],
            "S1 Auxiliary": [11.0],
            "Summary Auxiliary": [99.0],
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
                "var_keys": ["Feature"],
                "x_layer": "Intensity",
            },
            "columns": {
                "obs": {"select": {"sample": "<sample>"}},
                "var": {"select": {"Feature": "Feature"}},
            },
            "layers": [
                {
                    "name": "Intensity",
                    "source": r"^(?P<sample>S\d+) Intensity$",
                },
                {
                    "name": "Auxiliary",
                    "source": r"^(?P<sample>.+) Auxiliary$",
                },
            ],
        }
    )

    with caplog.at_level(logging.WARNING):
        pieces = convert_wide(df, rule)

    assert pieces.obs.index.tolist() == ["S1"]
    assert pieces.layers["Auxiliary"].shape == (1, 1)
    assert pieces.layers["Auxiliary"][0, 0] == 11.0
    assert "Summary" in caplog.text


def _structured_layer_rule(**layer_overrides: Any) -> ParseRule:
    """A single-sample wide rule whose auxiliary layer holds structured strings."""
    return ParseRule.model_validate(
        {
            "schema_version": "0.1",
            "file_version": "1",
            "software_name": "Synthetic",
            "software_version": "1.0",
            "input_shape": "wide",
            "quantification_level": "ion",
            "axis": {
                "obs_keys": ["sample"],
                "var_keys": ["Feature"],
                "x_layer": "Intensity",
                "duplicates": {"mode": "error"},
            },
            "columns": {
                "obs": {"select": {"sample": "<sample>"}},
                "var": {"select": {"Feature": "Feature"}},
            },
            "layers": [
                {"name": "Intensity", "source": r"^(?P<sample>S\d+) Intensity$"},
                {
                    "name": "Score",
                    "source": r"^(?P<sample>S\d+) Score$",
                    **layer_overrides,
                },
            ],
        }
    )


def test_convert_wide_value_pattern_extracts_numbers_from_structured_cells() -> None:
    df = pd.DataFrame(
        {
            "Feature": ["F1", "F2", "F3"],
            "S1 Intensity": [10.0, 20.0, 30.0],
            "S1 Score": ["site:12.5", None, "site:7.0"],
        }
    )

    pieces = convert_wide(
        df,
        _structured_layer_rule(value_pattern={"mode": "regex", "pattern": r":(-?\d+(?:\.\d+)?)$"}),
    )

    score = pieces.layers["Score"]
    assert score[0, 0] == 12.5
    assert np.isnan(score[0, 1])
    assert score[0, 2] == 7.0


def test_convert_wide_single_column_layer_does_not_fill_down_the_feature_axis() -> None:
    """A one-column-per-sample layer must never coalesce values between features.

    `_gather_layer_matrix` bfills across a sample's columns; with a single column that
    must be a no-op. pandas 2.3 got this wrong for nullable dtypes and filled along the
    feature axis instead, so `coerce_numeric` pins the result to plain float64.
    """
    df = pd.DataFrame(
        {
            "Feature": ["F1", "F2", "F3", "F4"],
            "S1 Intensity": [10.0, 20.0, 30.0, 40.0],
            "S1 Score": [None, None, "site:9.0", None],
        }
    )

    pieces = convert_wide(
        df,
        _structured_layer_rule(value_pattern={"mode": "regex", "pattern": r":(-?\d+(?:\.\d+)?)$"}),
    )

    score = pieces.layers["Score"]
    assert score[0, 2] == 9.0
    assert np.isnan(score[0, [0, 1, 3]]).all()


def test_convert_wide_warns_when_a_matched_layer_is_entirely_missing(
    caplog: pytest.LogCaptureFixture,
) -> None:
    df = pd.DataFrame(
        {
            "Feature": ["F1", "F2"],
            "S1 Intensity": [10.0, 20.0],
            "S1 Score": ["not-a-number", "also-not"],
        }
    )

    with caplog.at_level(logging.WARNING):
        pieces = convert_wide(df, _structured_layer_rule())

    assert np.isnan(pieces.layers["Score"]).all()
    assert "Score" in caplog.text
    assert "missing after numeric coercion" in caplog.text
