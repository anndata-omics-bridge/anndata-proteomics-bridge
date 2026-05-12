"""Convert a long-format DataFrame into AnnData pieces using a ParseRule."""

from __future__ import annotations

import numpy as np
import pandas as pd

from anndata_proteomics.converters._pieces import ConversionPieces
from anndata_proteomics.converters.factors import encode_factor
from anndata_proteomics.rules.schema import ParseRule

_KEY_SEPARATOR = "_"


def _join_keys(row: pd.Series) -> str:
    return _KEY_SEPARATOR.join(str(v) for v in row)


def _build_index(df: pd.DataFrame, keys: list[str]) -> pd.Series:
    """Build a string index from one or more key columns."""
    if len(keys) == 1:
        return df[keys[0]].astype(str)
    return df[keys].apply(_join_keys, axis=1)


def _build_axis_frame(
    df: pd.DataFrame, keys: list[str], output_columns: list[str]
) -> pd.DataFrame:
    """Take first occurrence per key tuple for already-materialized output columns."""
    needed_cols = list(dict.fromkeys(list(keys) + output_columns))
    block = df[needed_cols].drop_duplicates(subset=keys).copy()
    out = block[output_columns].copy()
    out.index = _build_index(block, keys).values
    out.index.name = _KEY_SEPARATOR.join(keys)
    return out


def _aggfunc_for(rule: ParseRule) -> str:
    mode = rule.axis.duplicates.mode
    if mode == "aggregate":
        return "sum"
    if mode == "keep_all_as_raw_table":
        raise NotImplementedError(
            "duplicates.mode='keep_all_as_raw_table' is not yet supported"
        )
    return "first"


def _pivot_layer(
    df: pd.DataFrame,
    obs_keys: list[str],
    var_keys: list[str],
    values: pd.Series,
    aggfunc: str,
) -> pd.DataFrame:
    work = df[obs_keys + var_keys].copy()
    work["__value__"] = values.values
    pivot = work.pivot_table(
        index=obs_keys, columns=var_keys, values="__value__", aggfunc=aggfunc
    )
    pivot.index = (
        pivot.index.astype(str)
        if len(obs_keys) == 1
        else pivot.index.to_frame().apply(_join_keys, axis=1).values
    )
    pivot.columns = (
        pivot.columns.astype(str)
        if len(var_keys) == 1
        else pivot.columns.to_frame().apply(_join_keys, axis=1).values
    )
    return pivot


def convert_long(df: pd.DataFrame, rule: ParseRule) -> ConversionPieces:
    """Convert a long DataFrame to AnnData pieces using a long ParseRule."""
    if rule.input_shape != "long":
        raise ValueError(f"convert_long called with {rule.input_shape!r} rule")

    obs_keys = list(rule.axis.obs_keys)
    var_keys = list(rule.axis.var_keys)

    obs_df = _build_axis_frame(df, obs_keys, rule.columns.obs.names)
    var_df = _build_axis_frame(df, var_keys, rule.columns.var.names)

    aggfunc = _aggfunc_for(rule)

    layers: dict[str, np.ndarray] = {}
    for layer in rule.layers:
        values = df[layer.source_column]
        if layer.encoding_mode == "factor":
            values = encode_factor(values, layer.categories or {})
        else:
            # Vendors sometimes use sentinels like "-" for missing in otherwise-numeric
            # columns; coerce so they become NaN rather than blowing up the pivot.
            values = pd.to_numeric(values, errors="coerce")

        pivot = _pivot_layer(df, obs_keys, var_keys, values, aggfunc)
        matrix = pivot.reindex(index=obs_df.index, columns=var_df.index)
        layers[layer.name] = matrix.values

    X = layers[rule.axis.x_layer]
    return ConversionPieces(X=X, obs=obs_df, var=var_df, layers=layers)
