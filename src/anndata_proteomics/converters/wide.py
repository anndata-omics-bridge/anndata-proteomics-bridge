"""Convert a wide-format DataFrame into AnnData pieces using a ParseRule."""

from __future__ import annotations

import logging
import re

import numpy as np
import pandas as pd

from anndata_proteomics.converters._axis import build_axis_frame, build_index
from anndata_proteomics.converters._pieces import ConversionPieces
from anndata_proteomics.converters.factors import encode_factor
from anndata_proteomics.converters.numeric import coerce_numeric
from anndata_proteomics.rules.schema import Layer, ParseRule

logger = logging.getLogger(__name__)
_SAMPLE_PLACEHOLDER = "<sample>"


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


def _gather_layer_matrix(
    df: pd.DataFrame,
    layer: Layer,
    sample_order: list[str],
    var_index: pd.Index,
    var_keys: list[str],
    duplicate_mode: str,
) -> np.ndarray:
    """Build (n_obs × n_var) matrix for a single wide layer."""
    matches = _matching_columns(list(df.columns), layer.source)
    sample_to_columns: dict[str, list[str]] = {}
    for column, sample in matches:
        sample_to_columns.setdefault(sample, []).append(column)

    n_obs = len(sample_order)
    n_var = len(var_index)
    matrix = np.full((n_obs, n_var), np.nan, dtype="float64")
    feature_index = build_index(df, var_keys)

    for i, sample in enumerate(sample_order):
        columns = sample_to_columns.get(sample, [])
        if not columns:
            continue
        if duplicate_mode == "error" and len(columns) > 1:
            raise ValueError(
                "duplicate observation-feature keys are not allowed when "
                f"axis.duplicates.mode='error'; layer {layer.name!r} has multiple "
                f"columns for sample {sample!r}: {columns}"
            )
        values = [_coerce_layer_series(df[column], layer) for column in columns]
        combined = pd.concat(values, axis=1)
        if duplicate_mode == "aggregate":
            series = combined.sum(axis=1)
        else:
            series = combined.bfill(axis=1).iloc[:, 0]
        series.index = feature_index
        grouped = series.groupby(level=0, sort=False)
        if duplicate_mode == "aggregate":
            series = grouped.sum()
        else:
            series = grouped.first()
        matrix[i, :] = series.reindex(var_index).to_numpy(dtype="float64")
    return matrix


def _coerce_layer_series(series: pd.Series, layer: Layer) -> pd.Series:
    if layer.encoding_mode == "factor":
        return encode_factor(series, layer.categories)
    return coerce_numeric(series, layer.missing_values)


def _raise_on_duplicate_features(df: pd.DataFrame, rule: ParseRule) -> None:
    """Reject repeated feature keys in a wide table when requested by the rule."""
    if rule.axis.duplicates.mode != "error":
        return
    keys = list(rule.axis.var_keys)
    valid = df[keys].notna().all(axis=1)
    duplicated = df.loc[valid, keys].duplicated(keep=False)
    if not duplicated.any():
        return
    examples = (
        df.loc[valid, keys].loc[duplicated].drop_duplicates().head(5).to_dict(orient="records")
    )
    raise ValueError(
        "duplicate observation-feature keys are not allowed when "
        f"axis.duplicates.mode='error'; examples: {examples}"
    )


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
    _raise_on_duplicate_features(df, rule)

    headers = list(df.columns)

    # The x-layer defines the observation axis. Optional auxiliary layers may expose
    # summary columns or malformed tokens; those must not expand the run axis.
    x_layer = next(layer for layer in rule.layers if layer.name == rule.axis.x_layer)
    sample_order = list(
        dict.fromkeys(sample for _, sample in _matching_columns(headers, x_layer.source))
    )
    sample_set = set(sample_order)

    if not sample_order:
        raise ValueError(
            f"no columns matched any layer pattern for rule {rule.software_name!r}; "
            f"layers: {[layer.source for layer in rule.layers]}"
        )

    var_df = build_axis_frame(df, list(rule.axis.var_keys), rule.columns.var.names)

    layers: dict[str, np.ndarray] = {}
    for layer in rule.layers:
        layer_matches = _matching_columns(headers, layer.source)
        extra_samples = list(
            dict.fromkeys(sample for _, sample in layer_matches if sample not in sample_set)
        )
        if extra_samples:
            logger.warning(
                "ignoring layer %r sample token(s) outside x-layer axis: %s",
                layer.name,
                extra_samples,
            )
        axis_matches = [
            (column, sample) for column, sample in layer_matches if sample in sample_set
        ]
        if not rule.layer_required(layer) and not axis_matches:
            logger.info(
                "skipping optional layer %r: no x-layer samples matched %r",
                layer.name,
                layer.source,
            )
            continue
        layers[layer.name] = _gather_layer_matrix(
            df,
            layer,
            sample_order,
            var_df.index,
            list(rule.axis.var_keys),
            rule.axis.duplicates.mode,
        )

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
