"""Matrix-native ProteoBench HYE intermediate calculations."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import overload

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from pydantic import BaseModel, ConfigDict

from anndata_proteomics.proteobench.config import ExpectedRatio, ModuleSettings
from anndata_proteomics.proteobench.contracts import (
    FloatArray,
    FloatDType,
    QuantMatrix,
)
from anndata_proteomics.proteobench.mapping import (
    map_reported_proteins,
    render_proteobench_features,
)
from anndata_proteomics.rules.schema import QuantificationLevel

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

    conditions: NDArray[np.str_]
    raw_files: tuple[str, ...]
    sample_names: tuple[str, ...]


class _ResultModel(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)


class AccessionMappingProvenance(_ResultModel):
    """Use of the bundled ProteoBench accession mapper."""

    asset: str
    sha256: str
    entries: int
    matched_token_occurrences: int
    unmatched_token_occurrences: int


class ProteinMappingProvenance(_ResultModel):
    """Protein normalization and species-mapping provenance."""

    species_mapper: dict[str, str]
    accession_mapper: AccessionMappingProvenance


@dataclass(frozen=True)
class IntermediateResult:
    """Feature-aligned storage table and ProteoBench-compatible legacy table."""

    varm: pd.DataFrame
    legacy: pd.DataFrame
    intermediate_hash: str
    protein_mapping: ProteinMappingProvenance


@dataclass(frozen=True)
class _SpeciesStatistics:
    """Species-dependent statistics shared by both output identities."""

    species: NDArray[np.object_]
    expected: NDArray[np.float64]
    epsilon: NDArray[np.float64]
    empirical_median: NDArray[np.float64]
    empirical_mean: NDArray[np.float64]


@dataclass(frozen=True)
class _ConditionStatistics:
    """Per-feature values and valid-observation counts for one condition."""

    values: dict[str, NDArray[np.float64]]
    count: NDArray[np.int64]


@dataclass(frozen=True)
class _LegacyComputation:
    matrix: QuantMatrix
    feature_ids: NDArray[np.str_]
    species_flags: dict[str, NDArray[np.bool_]]
    contaminants: NDArray[np.bool_]
    decoys: NDArray[np.bool_]
    module_settings: ModuleSettings
    design: RunDesign
    source_dtype: FloatDType
    conditions: tuple[str, ...]
    level: QuantificationLevel


@dataclass(frozen=True)
class _LegacyAssembly:
    derived: pd.DataFrame
    matrix: QuantMatrix
    feature_ids: NDArray[np.str_]
    pre_unique: NDArray[np.bool_]
    included: NDArray[np.bool_]
    design: RunDesign
    level: QuantificationLevel
    source_dtype: FloatDType
    species: tuple[str, ...]
    conditions: tuple[str, ...]


def align_runs(
    observations: pd.DataFrame,
    module_settings: ModuleSettings,
) -> RunDesign:
    """Validate and align the sample design added by ``apb annotate``."""
    required = ("sample_name", "condition")
    missing = [column for column in required if column not in observations.columns]
    if missing:
        raise ValueError(
            "ProteoBench scoring requires prior sample annotation; "
            f"missing obs columns: {missing}. Run 'apb annotate' first."
        )
    if observations.loc[:, list(required)].isna().any(axis=None):
        raise ValueError(
            "ProteoBench scoring requires complete sample annotation; "
            "obs['sample_name'] and obs['condition'] must not contain missing values"
        )

    observed_names = observations["sample_name"].astype(str).tolist()
    observed_conditions = observations["condition"].astype(str).tolist()
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

    conditions = np.asarray(observed_conditions, dtype=np.str_)
    sample_names = tuple(observed_names)
    raw_files = tuple(sample.raw_file for sample in matched)
    return RunDesign(conditions=conditions, raw_files=raw_files, sample_names=sample_names)


def compute_intermediate(
    matrix: QuantMatrix,
    feature_ids: pd.Index,
    reported_proteins: pd.Series,
    module_settings: ModuleSettings,
    design: RunDesign,
    level: QuantificationLevel,
) -> IntermediateResult:
    """Compute feature statistics and assemble the legacy ProteoBench table."""
    if not feature_ids.is_unique:
        raise ValueError("ProteoBench scoring requires unique var_names")
    if matrix.ndim != 2:
        raise ValueError("ProteoBench scoring requires a two-dimensional quantification matrix")
    if matrix.shape[0] != len(design.conditions):
        raise ValueError("quantification rows and aligned sample design have different lengths")
    if matrix.shape[1] != len(feature_ids):
        raise ValueError("quantification columns and feature identifiers have different lengths")
    if len(reported_proteins) != len(feature_ids):
        raise ValueError("reported proteins and feature identifiers have different lengths")

    features = feature_ids.to_series()
    normalized_proteins = reported_proteins.astype("string")
    mapping_result = map_reported_proteins(normalized_proteins)
    proteins = mapping_result.proteins
    compatibility_features = render_proteobench_features(features)
    compatibility_feature_ids = compatibility_features.to_numpy(dtype=str)
    species_flags = {
        species: proteins.str.contains(flag, regex=True, na=False).to_numpy(dtype=bool)
        for flag, species in module_settings.species_mapper.items()
    }
    unique = np.sum(np.vstack(list(species_flags.values())), axis=0, dtype=np.int64)
    contaminants = _contaminants(proteins)
    decoys = np.zeros(len(feature_ids), dtype=np.bool_)

    source_dtype = np.float32 if _is_float32_backed(matrix) else np.float64
    conditions = tuple(sorted(set(design.conditions.tolist())))
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
        feature_ids.copy(),
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
            matrix=matrix,
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
        protein_mapping=ProteinMappingProvenance(
            species_mapper=dict(module_settings.species_mapper),
            accession_mapper=AccessionMappingProvenance(
                asset="ProteoBench mapper.csv",
                sha256=mapping_result.mapper_sha256,
                entries=mapping_result.mapper_entries,
                matched_token_occurrences=mapping_result.matched_token_occurrences,
                unmatched_token_occurrences=mapping_result.unmatched_token_occurrences,
            ),
        ),
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
        inputs.matrix,
        group_codes,
        len(unique_features),
        eligible,
        inputs.source_dtype,
    )

    grouped_flags: dict[str, NDArray[np.bool_]] = {}
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
            species=tuple(inputs.module_settings.species_expected_ratio),
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
    conditions: tuple[str, ...],
    source_dtype: FloatDType,
) -> tuple[dict[str, NDArray[np.float64]], NDArray[np.int64]]:
    """Derive per-condition statistics without imposing a feature identity."""
    statistics: dict[str, NDArray[np.float64]] = {}
    condition_counts: dict[str, NDArray[np.int64]] = {}
    for condition in conditions:
        rows = np.flatnonzero(design.conditions == condition)
        condition_statistics = _condition_statistics(matrix, rows, source_dtype)
        condition_counts[condition] = condition_statistics.count
        for metric, metric_values in condition_statistics.values.items():
            statistics[f"{metric}_{condition}"] = metric_values

    nr_observed = np.zeros(matrix.shape[1], dtype=np.int64)
    for counts in condition_counts.values():
        nr_observed += counts
    statistics["log2_A_vs_B"] = (
        statistics["log_Intensity_mean_A"] - statistics["log_Intensity_mean_B"]
    )
    return statistics, nr_observed


def _derive_species_statistics(
    fold_change: NDArray[np.float64],
    species_flags: dict[str, NDArray[np.bool_]],
    included: NDArray[np.bool_],
    expected_ratios: dict[str, ExpectedRatio],
) -> _SpeciesStatistics:
    """Derive species statistics after the caller applies its identity-specific mask."""
    # Species inclusion is controlled exclusively by the explicit boolean ``included``
    # mask. ``pd.NA`` is the table's missing-value encoding for excluded features; it is
    # serialized as an empty CSV cell only at the legacy compatibility boundary.
    species_values = np.full(len(fold_change), pd.NA, dtype=object)
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
    condition_statistics: dict[str, NDArray[np.float64]],
    nr_observed: NDArray[np.int64],
    species_flags: dict[str, NDArray[np.bool_]],
    unique: NDArray[np.int64],
    species_statistics: _SpeciesStatistics,
    conditions: tuple[str, ...],
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


def _collapse_positive_matrix(
    matrix: QuantMatrix,
    group_codes: NDArray[np.intp],
    n_groups: int,
    eligible: NDArray[np.bool_],
    dtype: FloatDType,
) -> FloatArray:
    """Sum positive canonical-feature values into ProteoBench feature groups."""
    if dtype is np.float32:
        collapsed = np.full((matrix.shape[0], n_groups), np.nan, dtype=np.float32)
    else:
        collapsed = np.full((matrix.shape[0], n_groups), np.nan, dtype=np.float64)
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
    rows: NDArray[np.intp],
    source_dtype: FloatDType,
) -> _ConditionStatistics:
    n_vars = matrix.shape[1]
    result: dict[str, NDArray[np.float64]] = {
        "log_Intensity_mean": np.full(n_vars, np.nan, dtype=np.float64),
        "log_Intensity_std": np.full(n_vars, np.nan, dtype=np.float64),
        "Intensity_mean": np.full(n_vars, np.nan, dtype=np.float64),
        "Intensity_std": np.full(n_vars, np.nan, dtype=np.float64),
        "CV": np.full(n_vars, np.nan, dtype=np.float64),
    }
    counts = np.zeros(n_vars, dtype=np.int64)
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
        counts[section] = count
    return _ConditionStatistics(values=result, count=counts)


def _mean_and_sample_std(
    values: FloatArray,
    valid: NDArray[np.bool_],
    count: NDArray[np.int64],
) -> tuple[FloatArray, NDArray[np.float64]]:
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
    fold_change: NDArray[np.float64],
    species: NDArray[np.object_],
    included: NDArray[np.bool_],
) -> dict[str, NDArray[np.float64]]:
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


def _contaminants(proteins: pd.Series) -> NDArray[np.bool_]:
    return proteins.str.contains("Cont_", regex=False, na=False).to_numpy(dtype=bool)


def _is_float32_backed(matrix: QuantMatrix) -> bool:
    values = np.asarray(matrix).ravel() if isinstance(matrix, np.ndarray) else matrix.data
    finite = np.asarray(values)[np.isfinite(values)]
    if not finite.size:
        return False
    return np.array_equal(finite, finite.astype(np.float32).astype(finite.dtype))


def _matrix_block(
    matrix: QuantMatrix,
    rows: NDArray[np.intp],
    start: int,
    stop: int,
    dtype: FloatDType,
) -> FloatArray:
    if isinstance(matrix, np.ndarray):
        block = matrix[rows, start:stop]
    else:
        block = matrix[rows, start:stop].toarray()
    return _as_float_array(block, dtype)


def _matrix_row(
    matrix: QuantMatrix,
    row: int,
    dtype: FloatDType,
) -> FloatArray:
    if isinstance(matrix, np.ndarray):
        values = matrix[row, :]
    else:
        values = matrix[row, :].toarray()
    return _as_float_array(values, dtype).reshape(-1).copy()


@overload
def _as_float_array(values: FloatArray, dtype: type[np.float32]) -> NDArray[np.float32]: ...


@overload
def _as_float_array(values: FloatArray, dtype: type[np.float64]) -> NDArray[np.float64]: ...


def _as_float_array(values: FloatArray, dtype: FloatDType) -> FloatArray:
    """Materialize values in one of the two supported computation precisions."""
    if dtype is np.float32:
        return np.asarray(values, dtype=np.float32)
    return np.asarray(values, dtype=np.float64)
