"""Matrix-native ProteoBench HYE intermediate calculations."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData

from anndata_proteomics._matrix_types import QuantMatrix, is_sparse_matrix
from anndata_proteomics.proteobench.config import ExpectedRatio, ModuleSettings
from anndata_proteomics.proteobench.mapping import (
    map_reported_proteins,
    render_proteobench_features,
)
from anndata_proteomics.proteobench.resolve import ResolvedRoles

_CHUNK_SIZE = 50_000
_CONDITION_METRICS = (
    "log_Intensity_mean",
    "log_Intensity_std",
    "Intensity_mean",
    "Intensity_std",
    "CV",
)
# ProteoBench's own intermediate names its feature column after the module level; levels it has no
# module for keep the plain level name.
_LEGACY_FEATURE_COLUMN = {"ion": "precursor ion", "peptidoform": "peptidoform"}


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


@dataclass(frozen=True)
class _SpeciesStatistics:
    """Species-dependent statistics shared by both output identities."""

    species: np.ndarray
    expected: np.ndarray
    epsilon: np.ndarray
    empirical_median: np.ndarray
    empirical_mean: np.ndarray


@dataclass(frozen=True)
class _LegacyComputation:
    target: AnnData
    feature_ids: np.ndarray
    species_flags: dict[str, np.ndarray]
    contaminants: np.ndarray
    decoys: np.ndarray
    module_settings: ModuleSettings
    design: RunDesign
    source_dtype: type[np.floating[Any]]
    conditions: list[str]
    level: str


@dataclass(frozen=True)
class _LegacyAssembly:
    derived: pd.DataFrame
    matrix: QuantMatrix
    feature_ids: np.ndarray
    pre_unique: np.ndarray
    included: np.ndarray
    design: RunDesign
    level: str
    source_dtype: type[np.floating[Any]]
    species: list[str]
    conditions: list[str]


def align_runs(
    target: AnnData,
    module_settings: ModuleSettings,
) -> RunDesign:
    """Validate and align the sample design added by ``apb annotate``."""
    obs = target.obs
    if not isinstance(obs, pd.DataFrame):
        raise TypeError("ProteoBench scoring requires an in-memory obs DataFrame")
    required = ("sample_name", "condition")
    missing = [column for column in required if column not in obs.columns]
    if missing:
        raise ValueError(
            "ProteoBench scoring requires prior sample annotation; "
            f"missing obs columns: {missing}. Run 'apb annotate' first."
        )
    if obs.loc[:, list(required)].isna().any(axis=None):
        raise ValueError(
            "ProteoBench scoring requires complete sample annotation; "
            "obs['sample_name'] and obs['condition'] must not contain missing values"
        )

    observed_names = obs["sample_name"].astype(str).tolist()
    observed_conditions = obs["condition"].astype(str).tolist()
    by_sample_name = {sample.sample_name: sample for sample in module_settings.samples}

    matched = []
    for sample_name, condition in zip(observed_names, observed_conditions, strict=True):
        sample = by_sample_name.get(sample_name)
        if sample is None:
            raise ValueError(
                f"annotated sample_name {sample_name!r} does not match any "
                "[[samples]].sample_name entry in the module TOML"
            )
        if condition != sample.condition:
            raise ValueError(
                f"annotated condition {condition!r} for sample_name {sample_name!r} "
                f"does not match module condition {sample.condition!r}"
            )
        matched.append(sample)

    if len(observed_names) != len(set(observed_names)):
        raise ValueError("annotated sample_name values do not map one-to-one to module samples")
    expected = set(by_sample_name)
    actual = set(observed_names)
    if actual != expected:
        missing = sorted(expected - actual)
        extra = sorted(actual - expected)
        raise ValueError(
            "annotated sample alignment is incomplete; "
            f"missing sample_name values={missing}, extra={extra}"
        )

    conditions = np.asarray(observed_conditions, dtype=object)
    sample_names = tuple(observed_names)
    raw_files = tuple(sample.raw_file for sample in matched)
    return RunDesign(conditions=conditions, raw_files=raw_files, sample_names=sample_names)


def compute_intermediate(
    target: AnnData,
    module_settings: ModuleSettings,
    roles: ResolvedRoles,
    design: RunDesign,
    *,
    level: str,
) -> IntermediateResult:
    """Compute feature statistics and assemble the legacy ProteoBench table."""
    if not target.var_names.is_unique:
        raise ValueError("ProteoBench scoring requires unique var_names")

    # The feature axis is `var_names` at every level (the rule's joined `axis.var_keys`), and the
    # uniqueness check above is therefore also the feature-identifier uniqueness check.
    features = target.var_names.to_series()

    reported_proteins = target.var[roles.proteins].astype("string").fillna("")
    mapping_result = map_reported_proteins(reported_proteins)
    proteins = mapping_result.proteins
    compatibility_features = render_proteobench_features(features)
    compatibility_feature_ids = compatibility_features.to_numpy(dtype=str)
    species_flags = {
        species: proteins.str.contains(flag, regex=True, na=False).to_numpy(dtype=bool)
        for flag, species in module_settings.species_mapper.items()
    }
    unique = np.sum(np.vstack(list(species_flags.values())), axis=0, dtype=np.int64)
    contaminants = _contaminants(proteins)
    decoys = np.zeros(target.n_vars, dtype=bool)

    matrix = _require_quant_matrix(target)
    source_dtype = np.float32 if _is_float32_backed(matrix) else np.float64
    conditions = sorted(set(design.conditions.tolist()))
    stats, nr_observed = _derive_condition_statistics(
        matrix,
        design,
        conditions,
        source_dtype,
    )

    multi_species = unique > module_settings.general.min_count_multispec
    pre_unique = (nr_observed > 0) & ~contaminants & ~decoys & ~multi_species
    included = pre_unique & (unique == 1)

    species_stats = _derive_species_statistics(
        stats["log2_A_vs_B"],
        species_flags,
        included,
        module_settings.species_expected_ratio,
    )
    varm = _assemble_feature_statistics(
        target.var_names.copy(),
        stats,
        nr_observed,
        species_flags,
        unique,
        species_stats,
        conditions,
    )
    varm["included"] = included

    legacy = _compute_legacy_intermediate(
        _LegacyComputation(
            target=target,
            feature_ids=compatibility_feature_ids,
            species_flags=species_flags,
            contaminants=contaminants,
            decoys=decoys,
            module_settings=module_settings,
            design=design,
            source_dtype=source_dtype,
            conditions=conditions,
            level=level,
        )
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


def _compute_legacy_intermediate(inputs: _LegacyComputation) -> pd.DataFrame:
    """Reproduce ProteoBench's legacy feature grouping and intermediate."""
    unique_features, group_codes = np.unique(inputs.feature_ids, return_inverse=True)
    canonical_unique = np.sum(
        np.vstack(list(inputs.species_flags.values())),
        axis=0,
        dtype=np.int64,
    )
    eligible = (
        ~inputs.contaminants
        & ~inputs.decoys
        & (canonical_unique <= inputs.module_settings.general.min_count_multispec)
    )
    matrix = _collapse_positive_matrix(
        _require_quant_matrix(inputs.target),
        group_codes,
        len(unique_features),
        eligible,
        inputs.source_dtype,
    )

    grouped_flags: dict[str, np.ndarray] = {}
    for species_name, flags in inputs.species_flags.items():
        grouped = np.zeros(len(unique_features), dtype=bool)
        np.logical_or.at(grouped, group_codes[eligible], flags[eligible])
        grouped_flags[species_name] = grouped
    unique = np.sum(np.vstack(list(grouped_flags.values())), axis=0, dtype=np.int64)

    stats, nr_observed = _derive_condition_statistics(
        matrix,
        inputs.design,
        inputs.conditions,
        inputs.source_dtype,
    )

    pre_unique = nr_observed > 0
    included = pre_unique & (unique == 1)
    species_stats = _derive_species_statistics(
        stats["log2_A_vs_B"],
        grouped_flags,
        included,
        inputs.module_settings.species_expected_ratio,
    )
    derived = _assemble_feature_statistics(
        pd.Index(unique_features),
        stats,
        nr_observed,
        grouped_flags,
        unique,
        species_stats,
        inputs.conditions,
    )

    return _assemble_legacy_intermediate(
        _LegacyAssembly(
            derived=derived,
            matrix=matrix,
            feature_ids=unique_features,
            pre_unique=pre_unique,
            included=included,
            design=inputs.design,
            level=inputs.level,
            source_dtype=inputs.source_dtype,
            species=list(inputs.module_settings.species_expected_ratio),
            conditions=inputs.conditions,
        )
    )


def _assemble_legacy_intermediate(inputs: _LegacyAssembly) -> pd.DataFrame:
    """Reconstruct ProteoBench's ``result_performance.csv`` representation."""
    candidate_order = np.flatnonzero(inputs.pre_unique)
    candidate_order = candidate_order[
        np.argsort(inputs.feature_ids[candidate_order], kind="stable")
    ]
    legacy_index = {
        feature_index: position for position, feature_index in enumerate(candidate_order)
    }

    selected = np.flatnonzero(inputs.included)
    selected = selected[np.argsort(inputs.feature_ids[selected], kind="stable")]
    index = [legacy_index[feature_index] for feature_index in selected]
    legacy = pd.DataFrame(index=index)
    legacy[_LEGACY_FEATURE_COLUMN.get(inputs.level, inputs.level)] = inputs.feature_ids[selected]

    for metric in _CONDITION_METRICS:
        for condition in inputs.conditions:
            column = f"{metric}_{condition}"
            legacy[column] = inputs.derived[column].to_numpy()[selected]
    legacy["log2_A_vs_B"] = inputs.derived["log2_A_vs_B"].to_numpy()[selected]

    raw_order = sorted(
        range(len(inputs.design.raw_files)),
        key=lambda row: inputs.design.raw_files[row],
    )
    for row in raw_order:
        values = _matrix_row(inputs.matrix, row, inputs.source_dtype)
        values[~np.isfinite(values) | (values <= 0)] = np.nan
        legacy[inputs.design.raw_files[row]] = values[selected]

    legacy["nr_observed"] = inputs.derived["nr_observed"].to_numpy()[selected]
    for species_name in inputs.species:
        legacy[species_name] = inputs.derived[species_name].to_numpy()[selected]
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
        legacy[column] = inputs.derived[column].to_numpy()[selected]
    return legacy


def _derive_condition_statistics(
    matrix: QuantMatrix,
    design: RunDesign,
    conditions: list[str],
    source_dtype: type[np.floating[Any]],
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """Derive per-condition statistics without imposing a feature identity."""
    statistics: dict[str, np.ndarray] = {}
    condition_counts: dict[str, np.ndarray] = {}
    for condition in conditions:
        rows = np.flatnonzero(design.conditions == condition)
        values = _condition_statistics(matrix, rows, source_dtype)
        condition_counts[condition] = values.pop("count")
        for metric, metric_values in values.items():
            statistics[f"{metric}_{condition}"] = metric_values

    nr_observed = np.zeros(matrix.shape[1], dtype=np.int64)
    for counts in condition_counts.values():
        nr_observed += counts
    statistics["log2_A_vs_B"] = (
        statistics["log_Intensity_mean_A"] - statistics["log_Intensity_mean_B"]
    )
    return statistics, nr_observed


def _derive_species_statistics(
    fold_change: np.ndarray,
    species_flags: dict[str, np.ndarray],
    included: np.ndarray,
    expected_ratios: dict[str, ExpectedRatio],
) -> _SpeciesStatistics:
    """Derive species statistics after the caller applies its identity-specific mask."""
    species_values = np.full(len(fold_change), "", dtype=object)
    expected = np.full(len(fold_change), np.nan, dtype=np.float64)
    for species_name, ratio in expected_ratios.items():
        selected = included & species_flags[species_name]
        species_values[selected] = species_name
        expected[selected] = np.log2(ratio.a_vs_b)

    centers = _empirical_centers(fold_change, species_values, included)
    return _SpeciesStatistics(
        species=species_values,
        expected=expected,
        epsilon=fold_change - expected,
        empirical_median=centers["median"],
        empirical_mean=centers["mean"],
    )


def _assemble_feature_statistics(
    index: pd.Index,
    condition_statistics: dict[str, np.ndarray],
    nr_observed: np.ndarray,
    species_flags: dict[str, np.ndarray],
    unique: np.ndarray,
    species_statistics: _SpeciesStatistics,
    conditions: list[str],
) -> pd.DataFrame:
    """Assemble statistics common to canonical and compatibility feature identities."""
    frame = pd.DataFrame(index=index)
    for metric in _CONDITION_METRICS:
        for condition in conditions:
            column = f"{metric}_{condition}"
            frame[column] = condition_statistics[column]
    fold_change = condition_statistics["log2_A_vs_B"]
    frame["log2_A_vs_B"] = fold_change
    frame["nr_observed"] = nr_observed
    for species_name, flags in species_flags.items():
        frame[species_name] = flags
    frame["unique"] = unique
    frame["species"] = species_statistics.species
    frame["log2_expectedRatio"] = species_statistics.expected
    frame["epsilon"] = species_statistics.epsilon
    frame["log2_empirical_median"] = species_statistics.empirical_median
    frame["log2_empirical_mean"] = species_statistics.empirical_mean
    frame["epsilon_precision_median"] = fold_change - species_statistics.empirical_median
    frame["epsilon_precision_mean"] = fold_change - species_statistics.empirical_mean
    return frame


def _require_quant_matrix(target: AnnData) -> QuantMatrix:
    """Return ``X`` as in-memory values, rejecting an absent or backed matrix."""
    matrix = target.X
    if matrix is None:
        raise ValueError("ProteoBench scoring requires quantitative values in X")
    if not isinstance(matrix, np.ndarray) and not is_sparse_matrix(matrix):
        raise TypeError("ProteoBench scoring requires an in-memory X; load the object first")
    return matrix


def _collapse_positive_matrix(
    matrix: QuantMatrix,
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
    matrix: QuantMatrix,
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


def _contaminants(proteins: pd.Series) -> np.ndarray:
    return proteins.str.contains("Cont_", regex=False, na=False).to_numpy(dtype=bool)


def _is_float32_backed(matrix: QuantMatrix) -> bool:
    values = matrix.data if is_sparse_matrix(matrix) else np.asarray(matrix).ravel()
    finite = np.asarray(values)[np.isfinite(values)]
    if not finite.size:
        return False
    return np.array_equal(finite, finite.astype(np.float32).astype(finite.dtype))


def _matrix_block(
    matrix: Any,  # SciPy slicing cannot be described structurally; see QuantMatrix.
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
    matrix: Any,  # SciPy slicing cannot be described structurally; see QuantMatrix.
    row: int,
    dtype: type[np.floating[Any]],
) -> np.ndarray:
    values = matrix[row, :]
    if is_sparse_matrix(values):
        values = values.toarray()
    return np.asarray(values, dtype=dtype).reshape(-1).copy()
