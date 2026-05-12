"""Convert a wide-format DataFrame into AnnData pieces using a ParseRule."""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

from anndata_proteomics.converters._pieces import ConversionPieces
from anndata_proteomics.converters.factors import encode_factor
from anndata_proteomics.rules.schema import Layer, ParseRule

_KEY_SEPARATOR = "_"
_SAMPLE_PLACEHOLDER = "<sample>"


def _join_keys(row: pd.Series) -> str:
    return _KEY_SEPARATOR.join(str(v) for v in row)


def _matching_columns(headers: list[str], pattern: str) -> list[tuple[str, str]]:
    """Return [(column, sample_token), ...] for columns matching `pattern`."""
    compiled = re.compile(pattern)
    out: list[tuple[str, str]] = []
    for h in headers:
        m = compiled.match(h)
        if m is None:
            continue
        try:
            sample = m.group("sample")
        except (IndexError, KeyError):
            sample = m.group(0)
        out.append((h, sample))
    return out


def _build_var_frame(df: pd.DataFrame, rule: ParseRule) -> pd.DataFrame:
    var_keys = list(rule.axis.var_keys)
    output_columns = rule.columns.var.names
    needed_cols = list(dict.fromkeys(var_keys + output_columns))
    block = df[needed_cols].drop_duplicates(subset=var_keys).copy()
    out = block[output_columns].copy()
    if len(var_keys) == 1:
        out.index = block[var_keys[0]].astype(str).values
    else:
        out.index = block[var_keys].apply(_join_keys, axis=1).values
    out.index.name = _KEY_SEPARATOR.join(var_keys)
    return out


def _gather_layer_matrix(
    df: pd.DataFrame, layer: Layer, sample_order: list[str], var_index: pd.Index
) -> np.ndarray:
    """Build (n_obs × n_var) matrix for a single wide layer."""
    matches = _matching_columns(list(df.columns), layer.column_pattern or "")
    sample_to_col = {sample: col for col, sample in matches}

    n_obs = len(sample_order)
    n_var = len(var_index)
    matrix = np.full((n_obs, n_var), np.nan, dtype="float64")

    for i, sample in enumerate(sample_order):
        col = sample_to_col.get(sample)
        if col is None:
            continue
        series = df[col]
        if layer.encoding_mode == "factor":
            series = encode_factor(series, layer.categories or {})
        else:
            # Coerce sentinels like "-" / empty strings to NaN so they fit the float matrix.
            series = pd.to_numeric(series, errors="coerce")
        matrix[i, :] = series.values[: n_var]
    return matrix


def _apply_sample_cleanup(samples: list[str], rule: ParseRule) -> list[str]:
    """Apply optional sample_name_cleanup pattern to sample tokens."""
    if rule.sample_name_cleanup is None or not rule.sample_name_cleanup.pattern:
        return samples
    pattern = re.compile(rule.sample_name_cleanup.pattern)
    out: list[str] = []
    for s in samples:
        m = pattern.search(s)
        out.append(m.group(1) if m and m.groups() else (m.group(0) if m else s))
    return out


def convert_wide(df: pd.DataFrame, rule: ParseRule) -> ConversionPieces:
    """Convert a wide DataFrame to AnnData pieces using a wide ParseRule."""
    if rule.input_shape != "wide":
        raise ValueError(f"convert_wide called with {rule.input_shape!r} rule")

    headers = list(df.columns)

    # Sample tokens — collect across all layers (union), preserve insertion order.
    sample_order: list[str] = []
    seen: set[str] = set()
    for layer in rule.layers:
        for _, sample in _matching_columns(headers, layer.column_pattern or ""):
            if sample not in seen:
                sample_order.append(sample)
                seen.add(sample)

    if not sample_order:
        raise ValueError(
            f"no columns matched any layer pattern for rule {rule.software_name!r}; "
            f"layers: {[layer.column_pattern for layer in rule.layers]}"
        )

    var_df = _build_var_frame(df, rule)

    layers: dict[str, np.ndarray] = {}
    for layer in rule.layers:
        layers[layer.name] = _gather_layer_matrix(df, layer, sample_order, var_df.index)

    obs_names = _apply_sample_cleanup(sample_order, rule)
    obs_index = pd.Index(obs_names, name="sample")
    obs_data: dict[str, list[str]] = {}
    for out_name, source in rule.columns.obs.select.items():
        if source == _SAMPLE_PLACEHOLDER:
            obs_data[out_name] = list(obs_names)
        else:
            raise ValueError(
                f"wide rule columns.obs entry {out_name!r} = {source!r}: "
                f"only the {_SAMPLE_PLACEHOLDER!r} placeholder is supported in wide shape"
            )
    obs_df = pd.DataFrame(obs_data, index=obs_index)

    X = layers[rule.axis.x_layer]
    return ConversionPieces(X=X, obs=obs_df, var=var_df, layers=layers)
