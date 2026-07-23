"""Tests for generic declarative computed columns."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd

from anndata_proteomics.converters.assemble import _compute_column, convert
from anndata_proteomics.rules.schema import ColumnCompute, ParseRule


def test_coalesce_uses_first_non_null_value_without_treating_empty_as_missing():
    frame = pd.DataFrame(
        {
            "primary": pd.Series(["P1", None, "", None], dtype="string"),
            "fallback": pd.Categorical(["F1", "F2", "F3", None]),
        }
    )
    compute = ColumnCompute.model_validate(
        {
            "name": "Proteins",
            "from": ["primary", "fallback"],
            "how": "coalesce",
        }
    )

    result = _compute_column(frame, compute)

    assert result.iloc[:3].tolist() == ["P1", "F2", ""]
    assert pd.isna(result.iloc[3])


def test_join_nonempty_skips_null_and_empty_values_in_source_order():
    frame = pd.DataFrame(
        {
            "leading": pd.Series(["P1", None, "", None], dtype="string"),
            "mapped": pd.Categorical(["P2;P3", "P4", "P5", None]),
        }
    )
    compute = ColumnCompute.model_validate(
        {
            "name": "Proteins",
            "from": ["leading", "mapped"],
            "how": "join_nonempty",
            "separator": ",",
        }
    )

    result = _compute_column(frame, compute)

    assert result.iloc[:3].tolist() == ["P1,P2;P3", "P4", "P5"]
    assert pd.isna(result.iloc[3])


def test_generic_string_computes_survive_h5ad_round_trip(tmp_path: Path) -> None:
    frame = pd.DataFrame(
        {
            "Feature": ["F1", "F2", "F3"],
            "Primary": pd.Categorical(["P1", None, ""]),
            "Fallback": pd.Series(["M1", "M2", None], dtype="string"),
            "S1 Intensity": [1.0, 2.0, 3.0],
        }
    )
    rule = ParseRule.model_validate(
        {
            "schema_version": "0.1",
            "file_version": "1",
            "software_name": "Synthetic",
            "software_version": "1",
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
                "var": {
                    "select": {
                        "Feature": "Feature",
                        "Primary": "Primary",
                        "Fallback": "Fallback",
                    },
                    "compute": [
                        {
                            "name": "Complete",
                            "from": ["Primary", "Fallback"],
                            "how": "coalesce",
                        },
                        {
                            "name": "Joined",
                            "from": ["Primary", "Fallback"],
                            "how": "join_nonempty",
                            "separator": ",",
                        },
                    ],
                },
            },
            "layers": [
                {
                    "name": "Intensity",
                    "source": "^(?P<sample>S\\d+) Intensity$",
                }
            ],
        }
    )

    output = tmp_path / "computed.h5ad"
    convert(frame, rule).write_h5ad(output)
    restored = ad.read_h5ad(output)

    assert restored.var["Complete"].tolist() == ["P1", "M2", ""]
    assert restored.var["Joined"].iloc[:2].tolist() == ["P1,M1", "M2"]
    assert pd.isna(restored.var["Joined"].iloc[2])
