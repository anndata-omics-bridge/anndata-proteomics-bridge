"""ProteoBench-compatible aggregate HYE score metrics."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict, Field, RootModel, model_validator
from scipy.stats import rankdata

PROTEOBENCH_COMPATIBILITY_VERSION = "0.17.0"
PROTEOBENCH_SOURCE_REVISION = "fc95e712ca0466485814d3895087a048cfc0d2b0"

type ScoreMetric = int | float
type MetricAggregate = Callable[[pd.Series], float]


class _ScoreModel(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True, ser_json_inf_nan="null")


class ScoreConfig(_ScoreModel):
    """Cutoff range used to build one ProteoBench score document."""

    default_cutoff: int = Field(default=3, ge=1)
    max_nr_observed: int = Field(default=6, ge=1)

    @model_validator(mode="after")
    def _validate_cutoff(self) -> ScoreConfig:
        if self.default_cutoff > self.max_nr_observed:
            raise ValueError("default_cutoff must not exceed max_nr_observed")
        return self


class CutoffScores(RootModel[dict[str, ScoreMetric]]):
    """All compatible score metrics calculated at one observation cutoff."""

    model_config = ConfigDict(frozen=True, ser_json_inf_nan="null")


class ProteoBenchScores(_ScoreModel):
    """Typed storage contract for one ProteoBench score document."""

    intermediate_hash: str
    results: dict[str, CutoffScores]
    proteobench_version: str
    median_abs_epsilon_global: float
    mean_abs_epsilon_global: float
    median_abs_epsilon_eq_species: float
    mean_abs_epsilon_eq_species: float
    median_abs_epsilon_precision_global: float
    mean_abs_epsilon_precision_global: float
    median_abs_epsilon_precision_eq_species: float
    mean_abs_epsilon_precision_eq_species: float
    nr_feature: int


def build_scores(
    intermediate: pd.DataFrame,
    intermediate_hash: str,
    config: ScoreConfig,
) -> ProteoBenchScores:
    """Compute the compatible score-only ProteoBench storage model."""
    results = {
        str(cutoff): _metrics_at_cutoff(intermediate, cutoff)
        for cutoff in range(1, config.max_nr_observed + 1)
    }
    selected = results[str(config.default_cutoff)]
    selected_metrics = selected.root
    return ProteoBenchScores(
        intermediate_hash=intermediate_hash,
        results=results,
        proteobench_version=PROTEOBENCH_COMPATIBILITY_VERSION,
        median_abs_epsilon_global=_required_float_metric(
            selected_metrics, "median_abs_epsilon_global"
        ),
        mean_abs_epsilon_global=_required_float_metric(selected_metrics, "mean_abs_epsilon_global"),
        median_abs_epsilon_eq_species=_required_float_metric(
            selected_metrics, "median_abs_epsilon_eq_species"
        ),
        mean_abs_epsilon_eq_species=_required_float_metric(
            selected_metrics, "mean_abs_epsilon_eq_species"
        ),
        median_abs_epsilon_precision_global=_required_float_metric(
            selected_metrics, "median_abs_epsilon_precision_global"
        ),
        mean_abs_epsilon_precision_global=_required_float_metric(
            selected_metrics, "mean_abs_epsilon_precision_global"
        ),
        median_abs_epsilon_precision_eq_species=_required_float_metric(
            selected_metrics, "median_abs_epsilon_precision_eq_species"
        ),
        mean_abs_epsilon_precision_eq_species=_required_float_metric(
            selected_metrics, "mean_abs_epsilon_precision_eq_species"
        ),
        nr_feature=_required_int_metric(selected_metrics, "nr_feature"),
    )


def compute_roc_auc(frame: pd.DataFrame) -> float:
    """Compute binary ROC-AUC from absolute fold changes using average tie ranks."""
    required = {"species", "log2_A_vs_B", "log2_expectedRatio"}
    if frame.empty or not required <= set(frame.columns):
        return np.nan
    species_ratios = frame[["species", "log2_expectedRatio"]].drop_duplicates()
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


def _metrics_at_cutoff(frame: pd.DataFrame, cutoff: int) -> CutoffScores:
    selected = frame[frame["nr_observed"] >= cutoff]
    metrics: dict[str, ScoreMetric] = {
        **_epsilon_metrics(frame, cutoff, aggregation="median"),
        **_epsilon_metrics(frame, cutoff, aggregation="mean"),
        **_precision_metrics(frame, cutoff, aggregation="median"),
        **_precision_metrics(frame, cutoff, aggregation="mean"),
        **_cv_metrics(selected),
        "variance_epsilon_global": float(selected["epsilon"].var()) if len(selected) else 0.0,
        "nr_feature": len(selected),
        "roc_auc": compute_roc_auc(selected),
    }
    return CutoffScores(metrics)


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
        f"{aggregation}_abs_epsilon_eq_species": float(per_species.mean()),
        **{
            f"{aggregation}_abs_epsilon_{species}": float(value)
            for species, value in per_species.items()
        },
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
            f"{aggregation}_log2_empirical_{species}": float(value)
            for species, value in empirical_centers.items()
        },
        f"{aggregation}_abs_epsilon_precision_global": aggregate(precision),
        f"{aggregation}_abs_epsilon_precision_eq_species": float(per_species.mean()),
        **{
            f"{aggregation}_abs_epsilon_precision_{species}": float(value)
            for species, value in per_species.items()
        },
    }


def _cv_metrics(selected: pd.DataFrame) -> dict[str, float]:
    quantiles = selected[["CV_A", "CV_B"]].quantile([0.5, 0.75, 0.9, 0.95])
    averages = quantiles.mean(axis=1)
    return {
        "CV_median": float(averages.loc[0.50]),
        "CV_q75": float(averages.loc[0.75]),
        "CV_q90": float(averages.loc[0.90]),
        "CV_q95": float(averages.loc[0.95]),
    }


def _absolute_aggregate(name: str) -> MetricAggregate:
    if name == "median":
        return _median_absolute
    if name == "mean":
        return _mean_absolute
    raise ValueError(f"unsupported aggregation {name!r}")


def _median_absolute(values: pd.Series) -> float:
    return float(values.abs().median())


def _mean_absolute(values: pd.Series) -> float:
    return float(values.abs().mean())


def _required_float_metric(metrics: dict[str, ScoreMetric], name: str) -> float:
    value = metrics[name]
    if isinstance(value, bool):
        raise TypeError(f"score metric {name!r} is not numeric")
    return float(value)


def _required_int_metric(metrics: dict[str, ScoreMetric], name: str) -> int:
    value = metrics[name]
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"score metric {name!r} is not an integer")
    return value
