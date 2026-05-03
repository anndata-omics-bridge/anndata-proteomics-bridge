"""Tests for converters/factors.py."""

from __future__ import annotations

import pandas as pd

from anndata_proteomics.converters.factors import UNKNOWN_CODE, encode_factor


def test_encode_factor_known_values() -> None:
    s = pd.Series(["A", "B", "A", "C"])
    out = encode_factor(s, {"A": 1, "B": 2, "C": 3})
    assert list(out) == [1, 2, 1, 3]


def test_encode_factor_unknown_value_uses_default() -> None:
    s = pd.Series(["A", "X", "B"])
    out = encode_factor(s, {"A": 1, "B": 2})
    assert list(out) == [1, UNKNOWN_CODE, 2]


def test_encode_factor_nan_uses_default() -> None:
    s = pd.Series(["A", None, "B"])
    out = encode_factor(s, {"A": 1, "B": 2})
    assert list(out) == [1, UNKNOWN_CODE, 2]


def test_encode_factor_custom_default() -> None:
    s = pd.Series(["A", "X"])
    out = encode_factor(s, {"A": 1}, default=0)
    assert list(out) == [1, 0]


def test_encode_factor_returns_int64() -> None:
    s = pd.Series(["A"])
    out = encode_factor(s, {"A": 1})
    assert out.dtype == "int64"
