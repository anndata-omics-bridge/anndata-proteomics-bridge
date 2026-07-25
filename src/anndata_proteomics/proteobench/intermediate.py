"""Matrix-native ProteoBench HYE intermediate calculations."""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from anndata_proteomics._matrix_types import is_sparse_matrix
from anndata_proteomics.proteobench.config import ModuleSettings, SampleSettings
from anndata_proteomics.proteobench.mapping import (
    map_reported_proteins,
    render_proteobench_features,
)
from anndata_proteomics.proteobench.resolve import ResolvedRoles
from anndata_proteomics.rules.schema import ParseRule

_DEFAULT_RUN_CLEANUP = re.compile(
    r"(?:\.mzML\.gz|\.mzML|\.mgf|\.raw|\.RAW|\.d|\.wiff|_uncalibrated|"
    r" Intensity| Normalized Area)$"
)
_CHUNK_SIZE = 50_000


@dataclass(frozen=True)
class RunDesign:
    """Module sample design aligned to the target observation axis."""

    conditions: np.ndarray
    raw_files: tuple[str, ...]
    sample_names: tuple[str, ...]


@dataclass(frozen=True)
class IntermediateResult:
    """Feature-aligned storage table and ProteoBench-compatible legacy table."""

    varm: pd.DataFrame
    legacy: pd.DataFrame
    intermediate_hash: str
    protein_mapping: dict[str, Any]


def align_runs(
    target: Any,
    rule: ParseRule,
    roles: ResolvedRoles,
    module_settings: ModuleSettings,
) -> RunDesign:
    """Align converted observations to module samples without requiring annotation."""
    obs_column = roles.raw_file or rule.axis.obs_keys[0]
    observed = (
        target.obs[obs_column].astype(str).tolist()
        if obs_column in target.obs.columns
        else target.obs_names.astype(str).tolist()
    )
    cleanup = _DEFAULT_RUN_CLEANUP
    by_raw = {
        _clean_run_name(sample.raw_file, cleanup): sample for sample in module_settings.samples
    }
    by_sample_name = {sample.sample_name: sample for sample in module_settings.samples}
    wide = _wide_run_mapping(rule, module_settings, cleanup)

    matched: list[tuple[SampleSettings, str]] = []
    for value in observed:
        clean_value = _clean_run_name(value, cleanup)
        pair = wide.get(value) or wide.get(clean_value)
        if pair is None:
            sample = by_raw.get(clean_value) or by_sample_name.get(value)
            pair = (sample, _clean_run_name(sample.raw_file, cleanup)) if sample else None
        if pair is None:
            raise ValueError(
                f"converted run {value!r} does not match any [[samples]] entry in the module TOML"
            )
        matched.append(pair)

    samples = [pair[0] for pair in matched]
    if len({sample.raw_file for sample in samples}) != len(samples):
        raise ValueError("converted observations do not map one-to-one to module samples")
    expected = {sample.raw_file for sample in module_settings.samples}
    actual = {sample.raw_file for sample in samples}
    if actual != expected:
        missing = sorted(expected - actual)
        extra = sorted(actual - expected)
        raise ValueError(f"module sample alignment is incomplete; missing={missing}, extra={extra}")

    conditions = np.asarray([sample.condition for sample in samples], dtype=object)
    sample_names = tuple(sample.sample_name for sample in samples)
    _validate_existing_annotations(target, conditions, sample_names)
    raw_files = tuple(pair[1] for pair in matched)
    if len(raw_files) != len(set(raw_files)):
        raise ValueError(f"ProteoBench raw-file names are not unique after cleanup: {raw_files}")
    return RunDesign(conditions=conditions, raw_files=raw_files, sample_names=sample_names)


def compute_intermediate(  # noqa: C901, PLR0915 - scoring pipeline orchestration
    target: Any,
    module_settings: ModuleSettings,
    roles: ResolvedRoles,
    design: RunDesign,
) -> IntermediateResult:
    """Compute feature statistics and assemble the legacy ProteoBench table."""
    if not target.var_names.is_unique:
        raise ValueError("ProteoBench scoring requires unique var_names")

    feature_ids = target.var[roles.feature].astype(str).to_numpy()
    if len(set(feature_ids)) != len(feature_ids):
        raise ValueError(
            f"ProteoBench feature identifiers in var[{roles.feature!r}] are not unique"
        )

    reported_proteins = target.var[roles.proteins].astype("string").fillna("")
    mapping_result = map_reported_proteins(reported_proteins)
    proteins = mapping_result.proteins
    compatibility_features = render_proteobench_features(target.var[roles.feature])
    compatibility_feature_ids = compatibility_features.to_numpy(dtype=str)
    species_flags = {
        species: proteins.str.contains(flag, regex=True, na=False).to_numpy(dtype=bool)
        for flag, species in module_settings.species_mapper.items()
    }
    unique = np.sum(np.vstack(list(species_flags.values())), axis=0, dtype=np.int64)
    contaminants = _contaminants(proteins)
    decoys = np.zeros(target.n_vars, dtype=bool)

    matrix = target.X
    source_dtype = np.float32 if _is_float32_backed(matrix) else np.float64
    conditions = sorted(set(design.conditions.tolist()))
    stats: dict[str, np.ndarray] = {}
    condition_counts: dict[str, np.ndarray] = {}
    for condition in conditions:
        rows = np.flatnonzero(design.conditions == condition)
        values = _condition_statistics(matrix, rows, source_dtype)
        condition_counts[condition] = values.pop("count")
        for metric, metric_values in values.items():
            stats[f"{metric}_{condition}"] = metric_values

    nr_observed = np.zeros(target.n_vars, dtype=np.int64)
    for counts in condition_counts.values():
        nr_observed += counts
    stats["log2_A_vs_B"] = stats["log_Intensity_mean_A"] - stats["log_Intensity_mean_B"]

    multi_species = unique > module_settings.general.min_count_multispec
    pre_unique = (nr_observed > 0) & ~contaminants & ~decoys & ~multi_species
    included = pre_unique & (unique == 1)

    species_values = np.full(target.n_vars, "", dtype=object)
    expected = np.full(target.n_vars, np.nan, dtype=np.float64)
    for species, ratio in module_settings.species_expected_ratio.items():
        selected = included & species_flags[species]
        species_values[selected] = species
        expected[selected] = np.log2(ratio.a_vs_b)

    epsilon = stats["log2_A_vs_B"] - expected
    centers = _empirical_centers(stats["log2_A_vs_B"], species_values, included)

    varm = pd.DataFrame(index=target.var_names.copy())
    for metric in (
        "log_Intensity_mean",
        "log_Intensity_std",
        "Intensity_mean",
        "Intensity_std",
        "CV",
    ):
        for condition in conditions:
            varm[f"{metric}_{condition}"] = stats[f"{metric}_{condition}"]
    varm["log2_A_vs_B"] = stats["log2_A_vs_B"]
    varm["nr_observed"] = nr_observed
    for species, flags in species_flags.items():
        varm[species] = flags
    varm["unique"] = unique
    varm["species"] = species_values
    varm["log2_expectedRatio"] = expected
    varm["epsilon"] = epsilon
    varm["log2_empirical_median"] = centers["median"]
    varm["log2_empirical_mean"] = centers["mean"]
    varm["epsilon_precision_median"] = stats["log2_A_vs_B"] - centers["median"]
    varm["epsilon_precision_mean"] = stats["log2_A_vs_B"] - centers["mean"]
    varm["included"] = included

    legacy = _compute_legacy_intermediate(
        target,
        feature_ids=compatibility_feature_ids,
        species_flags=species_flags,
        contaminants=contaminants,
        decoys=decoys,
        module_settings=module_settings,
        design=design,
        source_dtype=source_dtype,
        conditions=conditions,
    )
    digest = hashlib.sha1(legacy.to_string().encode("utf-8")).hexdigest()
    return IntermediateResult(
        varm=varm,
        legacy=legacy,
        intermediate_hash=digest,
        protein_mapping={
            "species_mapper": dict(module_settings.species_mapper),
            "accession_mapper": {
                "asset": "ProteoBench mapper.csv",
                "sha256": mapping_result.mapper_sha256,
                "entries": mapping_result.mapper_entries,
                "matched_token_occurrences": mapping_result.matched_token_occurrences,
                "unmatched_token_occurrences": mapping_result.unmatched_token_occurrences,
            },
        },
    )


def _compute_legacy_intermediate(  # noqa: PLR0913 - legacy scoring contract
    target: Any,
    *,
    feature_ids: np.ndarray,
    species_flags: dict[str, np.ndarray],
    contaminants: np.ndarray,
    decoys: np.ndarray,
    module_settings: ModuleSettings,
    design: RunDesign,
    source_dtype: type[np.floating[Any]],
    conditions: list[str],
) -> pd.DataFrame:
    """Reproduce ProteoBench's legacy feature grouping and intermediate."""
    unique_features, group_codes = np.unique(feature_ids, return_inverse=True)
    canonical_unique = np.sum(np.vstack(list(species_flags.values())), axis=0, dtype=np.int64)
    eligible = (
        ~contaminants & ~decoys & (canonical_unique <= module_settings.general.min_count_multispec)
    )
    matrix = _collapse_positive_matrix(
        target.X,
        group_codes,
        len(unique_features),
        eligible,
        source_dtype,
    )

    grouped_flags: dict[str, np.ndarray] = {}
    for species_name, flags in species_flags.items():
        grouped = np.zeros(len(unique_features), dtype=bool)
        np.logical_or.at(grouped, group_codes[eligible], flags[eligible])
        grouped_flags[species_name] = grouped
    unique = np.sum(np.vstack(list(grouped_flags.values())), axis=0, dtype=np.int64)

    stats: dict[str, np.ndarray] = {}
    condition_counts: dict[str, np.ndarray] = {}
    for condition in conditions:
        rows = np.flatnonzero(design.conditions == condition)
        values = _condition_statistics(matrix, rows, source_dtype)
        condition_counts[condition] = values.pop("count")
        for metric, metric_values in values.items():
            stats[f"{metric}_{condition}"] = metric_values
    nr_observed = np.zeros(len(unique_features), dtype=np.int64)
    for counts in condition_counts.values():
        nr_observed += counts
    stats["log2_A_vs_B"] = stats["log_Intensity_mean_A"] - stats["log_Intensity_mean_B"]

    pre_unique = nr_observed > 0
    included = pre_unique & (unique == 1)
    species_values = np.full(len(unique_features), "", dtype=object)
    expected = np.full(len(unique_features), np.nan, dtype=np.float64)
    for species_name, ratio in module_settings.species_expected_ratio.items():
        selected = included & grouped_flags[species_name]
        species_values[selected] = species_name
        expected[selected] = np.log2(ratio.a_vs_b)
    epsilon = stats["log2_A_vs_B"] - expected
    centers = _empirical_centers(stats["log2_A_vs_B"], species_values, included)

    derived = pd.DataFrame(index=pd.Index(unique_features))
    for metric in (
        "log_Intensity_mean",
        "log_Intensity_std",
        "Intensity_mean",
        "Intensity_std",
        "CV",
    ):
        for condition in conditions:
            derived[f"{metric}_{condition}"] = stats[f"{metric}_{condition}"]
    derived["log2_A_vs_B"] = stats["log2_A_vs_B"]
    derived["nr_observed"] = nr_observed
    for species_name, flags in grouped_flags.items():
        derived[species_name] = flags
    derived["unique"] = unique
    derived["species"] = species_values
    derived["log2_expectedRatio"] = expected
    derived["epsilon"] = epsilon
    derived["log2_empirical_median"] = centers["median"]
    derived["log2_empirical_mean"] = centers["mean"]
    derived["epsilon_precision_median"] = stats["log2_A_vs_B"] - centers["median"]
    derived["epsilon_precision_mean"] = stats["log2_A_vs_B"] - centers["mean"]

    return assemble_legacy_intermediate(
        derived,
        matrix,
        feature_ids=unique_features,
        pre_unique=pre_unique,
        included=included,
        design=design,
        level=module_settings.general.level,
        source_dtype=source_dtype,
        species=list(module_settings.species_expected_ratio),
        conditions=conditions,
    )


def assemble_legacy_intermediate(  # noqa: PLR0913 - legacy scoring contract
    derived: pd.DataFrame,
    matrix: Any,
    *,
    feature_ids: np.ndarray,
    pre_unique: np.ndarray,
    included: np.ndarray,
    design: RunDesign,
    level: str,
    source_dtype: type[np.floating[Any]],
    species: list[str],
    conditions: list[str],
) -> pd.DataFrame:
    """Reconstruct ProteoBench's ``result_performance.csv`` representation."""
    candidate_order = np.flatnonzero(pre_unique)
    candidate_order = candidate_order[np.argsort(feature_ids[candidate_order], kind="stable")]
    legacy_index = {
        feature_index: position for position, feature_index in enumerate(candidate_order)
    }

    selected = np.flatnonzero(included)
    selected = selected[np.argsort(feature_ids[selected], kind="stable")]
    index = [legacy_index[feature_index] for feature_index in selected]
    legacy = pd.DataFrame(index=index)
    precursor_column = "precursor ion" if level == "ion" else "peptidoform"
    legacy[precursor_column] = feature_ids[selected]

    for metric in (
        "log_Intensity_mean",
        "log_Intensity_std",
        "Intensity_mean",
        "Intensity_std",
        "CV",
    ):
        for condition in conditions:
            column = f"{metric}_{condition}"
            legacy[column] = derived[column].to_numpy()[selected]
    legacy["log2_A_vs_B"] = derived["log2_A_vs_B"].to_numpy()[selected]

    raw_order = sorted(range(len(design.raw_files)), key=lambda row: design.raw_files[row])
    for row in raw_order:
        values = _matrix_row(matrix, row, source_dtype)
        values[~np.isfinite(values) | (values <= 0)] = np.nan
        legacy[design.raw_files[row]] = values[selected]

    legacy["nr_observed"] = derived["nr_observed"].to_numpy()[selected]
    for species_name in species:
        legacy[species_name] = derived[species_name].to_numpy()[selected]
    for column in (
        "unique",
        "species",
        "log2_expectedRatio",
        "epsilon",
        "log2_empirical_median",
        "log2_empirical_mean",
        "epsilon_precision_median",
        "epsilon_precision_mean",
    ):
        legacy[column] = derived[column].to_numpy()[selected]
    return legacy


def _collapse_positive_matrix(
    matrix: Any,
    group_codes: np.ndarray,
    n_groups: int,
    eligible: np.ndarray,
    dtype: type[np.floating[Any]],
) -> np.ndarray:
    """Sum positive canonical-feature values into ProteoBench feature groups."""
    collapsed = np.full((matrix.shape[0], n_groups), np.nan, dtype=dtype)
    for row in range(matrix.shape[0]):
        values = _matrix_row(matrix, row, dtype)
        valid = eligible & np.isfinite(values) & (values > 0)
        if not np.any(valid):
            continue
        totals = np.zeros(n_groups, dtype=dtype)
        np.add.at(totals, group_codes[valid], values[valid])
        present = np.zeros(n_groups, dtype=bool)
        np.logical_or.at(present, group_codes[valid], True)
        collapsed[row, present] = totals[present]
    return collapsed


def _condition_statistics(
    matrix: Any,
    rows: np.ndarray,
    source_dtype: type[np.floating[Any]],
) -> dict[str, np.ndarray]:
    n_vars = matrix.shape[1]
    result = {
        "log_Intensity_mean": np.full(n_vars, np.nan, dtype=np.float64),
        "log_Intensity_std": np.full(n_vars, np.nan, dtype=np.float64),
        "Intensity_mean": np.full(n_vars, np.nan, dtype=np.float64),
        "Intensity_std": np.full(n_vars, np.nan, dtype=np.float64),
        "CV": np.full(n_vars, np.nan, dtype=np.float64),
        "count": np.zeros(n_vars, dtype=np.int64),
    }
    for start in range(0, n_vars, _CHUNK_SIZE):
        stop = min(start + _CHUNK_SIZE, n_vars)
        block = _matrix_block(matrix, rows, start, stop, source_dtype)
        valid = np.isfinite(block) & (block > 0)
        count = np.count_nonzero(valid, axis=0)
        intensity_mean_native, intensity_std = _mean_and_sample_std(block, valid, count)
        logged = np.full_like(block, np.nan)
        np.log2(block, out=logged, where=valid)
        log_mean_native, log_std = _mean_and_sample_std(logged, valid, count)
        cv = np.divide(
            intensity_std,
            intensity_mean_native,
            out=np.full_like(intensity_std, np.nan),
            where=np.isfinite(intensity_mean_native) & (intensity_mean_native != 0),
        )
        section = slice(start, stop)
        result["Intensity_mean"][section] = intensity_mean_native.astype(np.float64)
        result["Intensity_std"][section] = intensity_std
        result["log_Intensity_mean"][section] = log_mean_native.astype(np.float64)
        result["log_Intensity_std"][section] = log_std
        result["CV"][section] = cv
        result["count"][section] = count
    return result


def _mean_and_sample_std(
    values: np.ndarray,
    valid: np.ndarray,
    count: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    native_sum = np.sum(np.where(valid, values, 0), axis=0, dtype=values.dtype)
    mean = np.divide(
        native_sum,
        count,
        out=np.full(values.shape[1], np.nan, dtype=values.dtype),
        where=count > 0,
    )
    centered = np.where(valid, values.astype(np.float64) - mean.astype(np.float64), 0.0)
    squared = np.sum(centered * centered, axis=0, dtype=np.float64)
    std = np.sqrt(
        np.divide(
            squared,
            count - 1,
            out=np.full(values.shape[1], np.nan, dtype=np.float64),
            where=count > 1,
        )
    )
    return mean, std


def _empirical_centers(
    fold_change: np.ndarray,
    species: np.ndarray,
    included: np.ndarray,
) -> dict[str, np.ndarray]:
    median = np.full(len(fold_change), np.nan, dtype=np.float64)
    mean = np.full(len(fold_change), np.nan, dtype=np.float64)
    frame = pd.DataFrame(
        {"fold_change": fold_change[included], "species": species[included]},
        index=np.flatnonzero(included),
    )
    if not frame.empty:
        median[frame.index] = frame.groupby("species")["fold_change"].transform("median")
        mean[frame.index] = frame.groupby("species")["fold_change"].transform("mean")
    return {"median": median, "mean": mean}


def _wide_run_mapping(
    rule: ParseRule,
    module_settings: ModuleSettings,
    cleanup: re.Pattern[str],
) -> dict[str, tuple[SampleSettings, str]]:
    if rule.input_shape != "wide":
        return {}
    x_layer = next(layer for layer in rule.layers if layer.name == rule.axis.x_layer)
    pattern = re.compile(x_layer.source)
    result: dict[str, tuple[SampleSettings, str]] = {}
    for sample in module_settings.samples:
        raw_column = sample.raw_file
        match = pattern.match(raw_column)
        observed_name = match.group("sample") if match else _clean_run_name(raw_column, cleanup)
        result[observed_name] = (sample, _clean_run_name(raw_column, cleanup, strip_path=False))
    return result


def _clean_run_name(
    value: str,
    cleanup: re.Pattern[str],
    *,
    strip_path: bool = True,
) -> str:
    name = str(value)
    if strip_path:
        name = name.replace("\\", "/").rsplit("/", 1)[-1]
    return cleanup.sub("", name)


def _validate_existing_annotations(
    target: Any,
    conditions: np.ndarray,
    sample_names: tuple[str, ...],
) -> None:
    for column, expected in (
        ("condition", conditions),
        ("sample_name", np.asarray(sample_names, dtype=object)),
    ):
        if column not in target.obs.columns:
            continue
        actual = target.obs[column].astype(str).to_numpy()
        if not np.array_equal(actual, expected.astype(str)):
            raise ValueError(f"existing obs[{column!r}] disagrees with the required module TOML")


def _contaminants(proteins: pd.Series) -> np.ndarray:
    return proteins.str.contains("Cont_", regex=False, na=False).to_numpy(dtype=bool)


def _is_float32_backed(matrix: Any) -> bool:
    values = matrix.data if is_sparse_matrix(matrix) else np.asarray(matrix).ravel()
    finite = np.asarray(values)[np.isfinite(values)]
    if not finite.size:
        return False
    return np.array_equal(finite, finite.astype(np.float32).astype(finite.dtype))


def _matrix_block(
    matrix: Any,
    rows: np.ndarray,
    start: int,
    stop: int,
    dtype: type[np.floating[Any]],
) -> np.ndarray:
    block = matrix[rows, start:stop]
    if is_sparse_matrix(block):
        block = block.toarray()
    return np.asarray(block, dtype=dtype)


def _matrix_row(
    matrix: Any,
    row: int,
    dtype: type[np.floating[Any]],
) -> np.ndarray:
    values = matrix[row, :]
    if is_sparse_matrix(values):
        values = values.toarray()
    return np.asarray(values, dtype=dtype).reshape(-1).copy()
