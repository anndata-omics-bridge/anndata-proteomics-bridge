"""Shared axis-frame and key-index helpers for the long/wide converters."""

from __future__ import annotations

import pandas as pd

KEY_SEPARATOR = "_"


def join_keys(row: pd.Series) -> str:
    """Join a row of key values into a single string index token."""
    return KEY_SEPARATOR.join(str(v) for v in row)


def build_index(df: pd.DataFrame, keys: list[str]) -> pd.Series:
    """Build a string index from one or more key columns.

    Vectorised string concatenation (not a row-wise apply) so it stays cheap on the
    full, un-deduplicated frame the long converter scatters from.
    """
    if len(keys) == 1:
        return df[keys[0]].astype(str)
    joined = df[keys[0]].astype(str)
    for key in keys[1:]:
        joined = joined + KEY_SEPARATOR + df[key].astype(str)
    return joined


def build_axis_frame(
    df: pd.DataFrame, keys: list[str], output_columns: list[str]
) -> pd.DataFrame:
    """Take first occurrence per key tuple for already-materialized output columns."""
    needed_cols = list(dict.fromkeys(list(keys) + output_columns))
    block = df[needed_cols].drop_duplicates(subset=keys).copy()
    out = block[output_columns].copy()
    out.index = build_index(block, keys).values
    out.index.name = KEY_SEPARATOR.join(keys)
    return out
