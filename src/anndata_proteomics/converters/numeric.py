"""Normalize vendor-reported numeric layer values."""

from __future__ import annotations

import pandas as pd


def coerce_numeric(series: pd.Series, missing_values: list[float]) -> pd.Series:
    """Coerce numeric values and replace rule-declared missing sentinels with NaN."""
    values = pd.to_numeric(series, errors="coerce")
    if missing_values:
        values = values.mask(values.isin(missing_values))
    return values
