"""Convert a long-format DataFrame into AnnData pieces using a ParseRule.

Each layer is built by scattering the long values into a dense (obs × var) matrix via
integer category codes, rather than ``DataFrame.pivot_table``. pivot_table materialises a
huge transient for high-cardinality var axes (the fragment level fans one report out to
millions of rows × hundreds of thousands of features and peaks at many GB); the scatter is
O(nnz + obs·var) and matches pivot_table's semantics exactly.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from anndata_proteomics.converters._axis import build_axis_frame, build_index
from anndata_proteomics.converters._pieces import ConversionPieces
from anndata_proteomics.converters.factors import encode_factor
from anndata_proteomics.rules.schema import ParseRule


def _aggfunc_for(rule: ParseRule) -> str:
    mode = rule.axis.duplicates.mode
    if mode == "aggregate":
        return "sum"
    if mode == "keep_all_as_raw_table":
        raise NotImplementedError("duplicates.mode='keep_all_as_raw_table' is not yet supported")
    return "first"


def _build_matrix(
    obs_codes: np.ndarray,
    var_codes: np.ndarray,
    values: np.ndarray,
    key_ok: np.ndarray,
    n_obs: int,
    n_var: int,
    aggfunc: str,
) -> np.ndarray:
    """Scatter ``values`` into a dense (n_obs × n_var) matrix.

    Mirrors ``pivot_table`` semantics: rows with a null axis key are dropped (as the
    pivot's groupby drops NaN keys); a cell with no contributing row stays NaN.
    - ``"first"``: keep the first non-null value in row order.
    - ``"sum"``: sum non-null values; a cell that has rows but only null values is 0.0
      (matching ``GroupBy.sum``), while a cell with no rows stays NaN.
    """
    matrix = np.full((n_obs, n_var), np.nan, dtype="float64")
    if aggfunc == "sum":
        finite = key_ok & ~np.isnan(values)
        totals = np.zeros((n_obs, n_var), dtype="float64")
        np.add.at(totals, (obs_codes[finite], var_codes[finite]), values[finite])
        present = np.zeros((n_obs, n_var), dtype=bool)
        present[obs_codes[key_ok], var_codes[key_ok]] = True
        matrix[present] = totals[present]
    else:  # "first" non-null: assign in reverse so the lowest-index value wins
        keep = key_ok & ~np.isnan(values)
        oc, vc, vv = obs_codes[keep], var_codes[keep], values[keep]
        matrix[oc[::-1], vc[::-1]] = vv[::-1]
    return matrix


def convert_long(df: pd.DataFrame, rule: ParseRule) -> ConversionPieces:
    """Convert a long DataFrame to AnnData pieces using a long ParseRule."""
    if rule.input_shape != "long":
        raise ValueError(f"convert_long called with {rule.input_shape!r} rule")

    obs_keys = list(rule.axis.obs_keys)
    var_keys = list(rule.axis.var_keys)

    obs_df = build_axis_frame(df, obs_keys, rule.columns.obs.names)
    var_df = build_axis_frame(df, var_keys, rule.columns.var.names)

    # Map every input row to its position in the obs/var axes. build_axis_frame keeps the
    # first occurrence per key, so the Categorical codes index directly into obs_df/var_df.
    obs_codes = pd.Categorical(build_index(df, obs_keys), categories=obs_df.index).codes
    var_codes = pd.Categorical(build_index(df, var_keys), categories=var_df.index).codes
    key_ok = df[obs_keys + var_keys].notna().all(axis=1).to_numpy()

    aggfunc = _aggfunc_for(rule)
    n_obs, n_var = len(obs_df), len(var_df)

    layers: dict[str, np.ndarray] = {}
    for layer in rule.layers:
        values = df[layer.source]
        if layer.encoding_mode == "factor":
            values = encode_factor(values, layer.categories)
        else:
            # Vendors sometimes use sentinels like "-" for missing in otherwise-numeric
            # columns; coerce so they become NaN rather than blowing up the scatter.
            values = pd.to_numeric(values, errors="coerce")

        layers[layer.name] = _build_matrix(
            obs_codes,
            var_codes,
            np.asarray(values, dtype="float64"),
            key_ok,
            n_obs,
            n_var,
            aggfunc,
        )

    X = layers[rule.axis.x_layer]
    return ConversionPieces(X=X, obs=obs_df, var=var_df, layers=layers)
