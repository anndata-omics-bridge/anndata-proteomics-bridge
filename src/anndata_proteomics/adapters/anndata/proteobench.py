"""AnnData and MuData storage adapter for ProteoBench scoring."""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass
from typing import cast

import numpy as np
import pandas as pd
from anndata import AnnData
from mudata import MuData
from numpy.typing import NDArray
from scipy.sparse import csc_array, csc_matrix, csr_array, csr_matrix

from anndata_proteomics.adapters.anndata.namespace import read_namespace, update_namespace
from anndata_proteomics.adapters.anndata.rules import require_stored_rule
from anndata_proteomics.proteobench.contracts import QuantMatrix, SparseQuantMatrix
from anndata_proteomics.proteobench.intermediate import ProteinMappingProvenance
from anndata_proteomics.proteobench.metrics import (
    PROTEOBENCH_COMPATIBILITY_VERSION,
    PROTEOBENCH_SOURCE_REVISION,
    ProteoBenchScores,
)
from anndata_proteomics.proteobench.pipeline import ProteoBenchResult
from anndata_proteomics.rules.schema import ParseRule
from anndata_proteomics.serialization import JsonObject, JsonValue, to_json_compatible
from anndata_proteomics.workflows.proteobench import ProteoBenchLevelInput

_STORAGE_SCHEMA_VERSION = "0.1"
_STORAGE_KEY = "proteobench"
_SPARSE_MATRIX_TYPES = (csr_matrix, csc_matrix, csr_array, csc_array)


@dataclass(frozen=True)
class ResolvedRoles:
    """Concrete AnnData locations used by ProteoBench scoring."""

    proteins: str
    feature: str = "var_names"
    intensity: str = "X"

    def as_dict(self) -> JsonObject:
        """Return the ProteoBench role mapping stored in APB metadata."""
        return {
            "Proteins": f"var:{self.proteins}",
            "feature": self.feature,
            "Intensity": self.intensity,
            "Sample name": "obs:sample_name",
            "Condition": "obs:condition",
        }


@dataclass(frozen=True, slots=True)
class ExtractedProteoBenchLevel:
    """Calculation input plus the physical roles needed to persist its result."""

    calculation: ProteoBenchLevelInput
    roles: ResolvedRoles


def resolve_targets(obj: AnnData | MuData) -> list[AnnData]:
    """Return the AnnData levels held by one supported container."""
    if isinstance(obj, MuData):
        if not obj.mod:
            raise ValueError("MuData has no modality to score")
        return [_require_anndata(name, target) for name, target in obj.mod.items()]
    return [obj]


def resolve_roles(target: AnnData) -> tuple[ParseRule, ResolvedRoles]:
    """Resolve ProteoBench storage locations from the stored APB rule."""
    rule = require_stored_rule(target)
    proteins = rule.column_roles.protein_assignment
    if proteins is None:
        raise ValueError(
            "ProteoBench scoring requires column_roles.protein_assignment in the stored APB rule"
        )
    if proteins not in target.var.columns:
        raise ValueError(
            "stored APB rule maps column_roles.protein_assignment to missing "
            f"var column {proteins!r}"
        )
    return rule, ResolvedRoles(proteins=proteins)


def extract_quant_matrix(target: AnnData) -> QuantMatrix:
    """Extract a supported in-memory floating-point matrix from AnnData ``X``."""
    matrix = target.X
    if matrix is None:
        raise ValueError("ProteoBench scoring requires quantitative values in X")
    if isinstance(matrix, np.ndarray):
        if matrix.dtype == np.dtype(np.float32):
            return cast(NDArray[np.float32], matrix)
        if matrix.dtype == np.dtype(np.float64):
            return cast(NDArray[np.float64], matrix)
        raise TypeError("ProteoBench scoring requires float32 or float64 values in X")
    if isinstance(matrix, _SPARSE_MATRIX_TYPES):
        if matrix.dtype not in (np.dtype(np.float32), np.dtype(np.float64)):
            raise TypeError("ProteoBench scoring requires float32 or float64 values in X")
        return cast(SparseQuantMatrix, matrix)
    raise TypeError("ProteoBench scoring requires an in-memory X; load the object first")


def extract_observations(target: AnnData) -> pd.DataFrame:
    """Extract the in-memory observation table used for design alignment."""
    observations = target.obs
    if not isinstance(observations, pd.DataFrame):
        raise TypeError("ProteoBench scoring requires an in-memory obs DataFrame")
    return observations


def store_result(
    target: AnnData,
    result: ProteoBenchResult,
    roles: ResolvedRoles,
) -> None:
    """Persist one calculated result in the existing APB storage schema."""
    namespace: dict[str, JsonValue] = dict(_stored_proteobench(target))
    target.varm[_STORAGE_KEY] = result.intermediate.varm
    updates: dict[str, JsonValue] = {
        "schema_version": _STORAGE_SCHEMA_VERSION,
        "compatibility_version": PROTEOBENCH_COMPATIBILITY_VERSION,
        "source_revision": PROTEOBENCH_SOURCE_REVISION,
        "column_roles": roles.as_dict(),
        "protein_mapping": _serialize_protein_mapping(result.intermediate.protein_mapping),
        "scores": _serialize_scores(result.scores),
    }
    namespace.update(updates)
    update_namespace(target, {_STORAGE_KEY: namespace})


def read_level(target: AnnData) -> ExtractedProteoBenchLevel:
    """Extract the exact calculation values and storage roles for one level."""
    _require_available_storage(target)
    rule, roles = resolve_roles(target)
    return ExtractedProteoBenchLevel(
        calculation=ProteoBenchLevelInput(
            observations=extract_observations(target),
            matrix=extract_quant_matrix(target),
            feature_ids=target.var_names.copy(),
            reported_proteins=target.var[roles.proteins].copy(),
            level=rule.quantification_level,
        ),
        roles=roles,
    )


def _require_available_storage(target: AnnData) -> None:
    """Refuse to overwrite either part of an existing ProteoBench result."""
    if _STORAGE_KEY in target.varm:
        raise ValueError("varm['proteobench'] already exists; refusing to overwrite scores")
    if "scores" in _stored_proteobench(target):
        raise ValueError(
            "uns['anndata_proteomics']['proteobench']['scores'] already exists; "
            "refusing to overwrite scores"
        )


def _stored_proteobench(target: AnnData) -> Mapping[str, JsonValue]:
    """Read the already stored ProteoBench sub-namespace, if any."""
    stored = read_namespace(target).get(_STORAGE_KEY)
    return stored if isinstance(stored, Mapping) else {}


def _require_anndata(name: str, target: object) -> AnnData:
    """Reject a MuData modality that is not an AnnData."""
    if not isinstance(target, AnnData):
        raise TypeError(f"MuData modality {name!r} is not an AnnData")
    return target


def _serialize_protein_mapping(result: ProteinMappingProvenance) -> JsonObject:
    """Serialize typed protein-mapping provenance at the storage boundary."""
    document = to_json_compatible(result.model_dump(mode="json"))
    if not isinstance(document, dict):
        raise TypeError("protein-mapping serialization did not produce an object")
    return document


def _serialize_scores(result: ProteoBenchScores) -> JsonObject:
    """Serialize typed score results at the storage boundary."""
    document = to_json_compatible(json.loads(result.model_dump_json()))
    if not isinstance(document, dict):
        raise TypeError("ProteoBench score serialization did not produce an object")
    return document
