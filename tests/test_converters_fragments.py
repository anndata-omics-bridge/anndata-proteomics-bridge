"""Unit tests for converters/_fragments.explode_fragments."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from anndata_proteomics.converters._fragments import explode_fragments
from anndata_proteomics.rules.schema import Fragments

_FRAGS = Fragments(
    label_column="Fragment.Info",
    value_columns=["Fragment.Quant.Raw", "Fragment.Quant.Corrected"],
    delimiter=";",
)


def _df() -> pd.DataFrame:
    # r1 has 2 fragments, r2 has 1; both lists carry DIA-NN's trailing ';' terminator.
    return pd.DataFrame(
        {
            "Run": ["r1", "r2"],
            "Precursor.Id": ["PEP2", "OTHER3"],
            "Fragment.Info": ["b4-unknown^1/327.16;y5-unknown^1/600.3;", "b2-unknown^1/200.1;"],
            "Fragment.Quant.Raw": ["100;200;", "9;"],
            "Fragment.Quant.Corrected": ["110;210;", "8;"],
        }
    )


def test_explodes_one_row_per_fragment_and_drops_trailing_empty() -> None:
    out = explode_fragments(_df(), _FRAGS)
    assert len(out) == 3  # 2 + 1, trailing empty terminator dropped
    assert list(out["fragment_label"]) == ["b4-unknown^1", "y5-unknown^1", "b2-unknown^1"]
    # parallel lists stay aligned and the precursor scalars replicate
    assert list(out["Fragment.Quant.Raw"]) == ["100", "200", "9"]
    assert list(out["Fragment.Quant.Corrected"]) == ["110", "210", "8"]
    assert list(out["Run"]) == ["r1", "r1", "r2"]
    assert list(out["Precursor.Id"]) == ["PEP2", "PEP2", "OTHER3"]


def test_label_is_token_before_slash() -> None:
    out = explode_fragments(_df(), _FRAGS)
    # the m/z after '/' is dropped from the label
    assert "/" not in "".join(out["fragment_label"])


def test_length_mismatch_raises() -> None:
    bad = pd.DataFrame(
        {
            "Run": ["r1"],
            "Fragment.Info": ["b4^1/1;y5^1/2;"],  # 2 fragments
            "Fragment.Quant.Raw": ["10;"],  # 1 value
            "Fragment.Quant.Corrected": ["10;20;"],
        }
    )
    with pytest.raises(ValueError, match="matching element counts"):
        explode_fragments(bad, _FRAGS)


def test_missing_column_raises() -> None:
    df = _df().drop(columns=["Fragment.Quant.Corrected"])
    with pytest.raises(KeyError, match="Fragment.Quant.Corrected"):
        explode_fragments(df, _FRAGS)


def test_precursor_with_no_fragments_is_dropped() -> None:
    df = _df()
    df.loc[1, ["Fragment.Info", "Fragment.Quant.Raw", "Fragment.Quant.Corrected"]] = [
        np.nan,
        np.nan,
        np.nan,
    ]
    out = explode_fragments(df, _FRAGS)
    assert len(out) == 2  # only r1's two fragments survive
    assert set(out["Run"]) == {"r1"}
