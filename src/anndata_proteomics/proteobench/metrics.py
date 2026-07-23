"""ProteoBench-compatible aggregate HYE score metrics."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import rankdata

PROTEOBENCH_COMPATIBILITY_VERSION = "0.17.0"
PROTEOBENCH_SOURCE_REVISION = "fc95e712ca0466485814d3895087a048cfc0d2b0"

_TOP_LEVEL_METRICS = (
    "median_abs_epsilon_global",
    "mean_abs_epsilon_global",
    "median_abs_epsilon_eq_species",
    "mean_abs_epsilon_eq_species",
    "median_abs_epsilon_precision_global",
    "mean_abs_epsilon_precision_global",
    "median_abs_epsilon_precision_eq_species",
    "mean_abs_epsilon_precision_eq_species",
    "nr_feature",
)


def build_scores(
    intermediate: pd.DataFrame,
    intermediate_hash: str,
    *,
    default_cutoff: int = 3,
    max_nr_observed: int = 6,
) -> dict[str, Any]:
    """Compute the compatible score-only ProteoBench JSON document."""
    if default_cutoff < 1 or default_cutoff > max_nr_observed:
        raise ValueError(f"default cutoff {default_cutoff} must be within 1..{max_nr_observed}")
    results = {
        str(cutoff): _metrics_at_cutoff(intermediate, cutoff)
        for cutoff in range(1, max_nr_observed + 1)
    }
    selected = results[str(default_cutoff)]
    payload: dict[str, Any] = {
        "intermediate_hash": intermediate_hash,
        "results": results,
        "proteobench_version": PROTEOBENCH_COMPATIBILITY_VERSION,
    }
    payload.update({key: selected[key] for key in _TOP_LEVEL_METRICS})
    return _json_compatible(payload)


def compute_roc_auc(frame: pd.DataFrame) -> float:
    """Compute binary ROC-AUC from absolute fold changes using average tie ranks."""
    required = {"species", "log2_A_vs_B", "log2_expectedRatio"}
    if frame.empty or not required <= set(frame.columns):
        return np.nan
    species_ratios = frame[["species", "log2_expectedRatio"]].drop_duplicates()
    if species_ratios.empty:
        return np.nan
    unchanged_index = species_ratios["log2_expectedRatio"].abs().idxmin()
    unchanged = species_ratios.loc[unchanged_index, "species"]

    y_true = (frame["species"] != unchanged).to_numpy(dtype=np.int8)
    y_score = frame["log2_A_vs_B"].abs().to_numpy(dtype=np.float64)
    valid = ~np.isnan(y_score)
    y_true = y_true[valid]
    y_score = y_score[valid]
    if len(y_true) < 2 or len(np.unique(y_true)) < 2:
        return np.nan

    ranks = rankdata(y_score, method="average")
    positives = y_true == 1
    n_positive = int(np.count_nonzero(positives))
    n_negative = len(y_true) - n_positive
    rank_sum = float(ranks[positives].sum())
    return (rank_sum - n_positive * (n_positive + 1) / 2) / (n_positive * n_negative)


def _metrics_at_cutoff(frame: pd.DataFrame, cutoff: int) -> dict[str, Any]:
    selected = frame[frame["nr_observed"] >= cutoff]
    metrics = {
        **_epsilon_metrics(frame, cutoff, aggregation="median"),
        **_epsilon_metrics(frame, cutoff, aggregation="mean"),
        **_precision_metrics(frame, cutoff, aggregation="median"),
        **_precision_metrics(frame, cutoff, aggregation="mean"),
        **_cv_metrics(selected),
        "variance_epsilon_global": selected["epsilon"].var() if len(selected) else 0.0,
        "nr_feature": len(selected),
        "roc_auc": compute_roc_auc(selected),
    }
    return metrics


def _epsilon_metrics(
    frame: pd.DataFrame,
    cutoff: int,
    *,
    aggregation: str,
) -> dict[str, float]:
    selected = frame[frame["nr_observed"] >= cutoff]
    aggregate = _absolute_aggregate(aggregation)
    per_species = selected.groupby("species")["epsilon"].apply(aggregate)
    return {
        f"{aggregation}_abs_epsilon_global": aggregate(selected["epsilon"]),
        f"{aggregation}_abs_epsilon_eq_species": per_species.mean(),
        **{f"{aggregation}_abs_epsilon_{species}": value for species, value in per_species.items()},
    }


def _precision_metrics(
    frame: pd.DataFrame,
    cutoff: int,
    *,
    aggregation: str,
) -> dict[str, float]:
    selected = frame[frame["nr_observed"] >= cutoff]
    center_name = "median" if aggregation == "median" else "mean"
    grouped = selected.groupby("species")["log2_A_vs_B"]
    centers = grouped.transform(center_name)
    precision = selected["log2_A_vs_B"] - centers
    aggregate = _absolute_aggregate(aggregation)

    precision_frame = selected[["species"]].copy()
    precision_frame["epsilon_precision"] = precision
    per_species = precision_frame.groupby("species")["epsilon_precision"].apply(aggregate)
    empirical_centers = grouped.agg(center_name)
    return {
        **{
            f"{aggregation}_log2_empirical_{species}": value
            for species, value in empirical_centers.items()
        },
        f"{aggregation}_abs_epsilon_precision_global": aggregate(precision),
        f"{aggregation}_abs_epsilon_precision_eq_species": per_species.mean(),
        **{
            f"{aggregation}_abs_epsilon_precision_{species}": value
            for species, value in per_species.items()
        },
    }


def _cv_metrics(selected: pd.DataFrame) -> dict[str, float]:
    quantiles = selected[["CV_A", "CV_B"]].quantile([0.5, 0.75, 0.9, 0.95])
    averages = quantiles.mean(axis=1)
    return {
        "CV_median": averages.loc[0.50],
        "CV_q75": averages.loc[0.75],
        "CV_q90": averages.loc[0.90],
        "CV_q95": averages.loc[0.95],
    }


def _absolute_aggregate(name: str):
    if name == "median":
        return lambda values: values.abs().median()
    if name == "mean":
        return lambda values: values.abs().mean()
    raise ValueError(f"unsupported aggregation {name!r}")


def _json_compatible(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_compatible(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_compatible(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return number if np.isfinite(number) else None
    if isinstance(value, np.bool_):
        return bool(value)
    return value
