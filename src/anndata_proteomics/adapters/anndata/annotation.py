"""AnnData and MuData persistence adapter for sample annotation."""

from __future__ import annotations

import json

import pandas as pd
from anndata import AnnData
from mudata import MuData
from pydantic import BaseModel, ConfigDict

from anndata_proteomics.adapters.anndata.namespace import (
    MissingNamespaceText,
    read_namespace_text,
    update_namespace,
)
from anndata_proteomics.annotation.loader import AnnotationFileOrigin
from anndata_proteomics.annotation.sample import SampleAnnotationProvenance
from anndata_proteomics.serialization import JsonValue, to_json_compatible
from anndata_proteomics.workflows.sample_annotation import SampleAnnotationResult

_PROVENANCE_KEY = "obs_annotations_json"


class _ObsAnnotationProvenance(BaseModel):
    """Typed base for one serialized ``obs`` annotation record."""

    model_config = ConfigDict(frozen=True)

    match_on: str
    key_field: str
    obs_columns_added: tuple[str, ...]
    n_obs_matched: int


class _FileObsAnnotationProvenance(_ObsAnnotationProvenance):
    """Stored provenance for an annotation table read from a file."""

    source: str
    source_format: str


class _InMemoryObsAnnotationProvenance(_ObsAnnotationProvenance):
    """Stored provenance for a programmatically constructed annotation table."""

    source: None = None
    source_format: None = None


def read_observation_frames(target: AnnData | MuData) -> tuple[pd.DataFrame, ...]:
    """Read the primary obs frame followed by every MuData modality obs frame."""
    return tuple(_observation_frame(holder) for holder in _observation_holders(target))


def write_sample_annotation(
    target: AnnData | MuData,
    result: SampleAnnotationResult,
) -> None:
    """Persist annotated obs frames and provenance to AnnData or MuData."""
    holders = _observation_holders(target)
    if len(holders) != len(result.observations):
        raise ValueError("sample annotation result does not match the container observation axes")

    for holder, annotated in zip(holders, result.observations, strict=True):
        holder.obs = annotated.frame

    update_namespace(
        target,
        {_PROVENANCE_KEY: _append_annotation_provenance(target, result.provenance)},
    )


def _observation_holders(target: AnnData | MuData) -> tuple[AnnData | MuData, ...]:
    """Return the container itself followed by every MuData modality."""
    if isinstance(target, MuData):
        return (target, *target.mod.values())
    return (target,)


def _observation_frame(holder: AnnData | MuData) -> pd.DataFrame:
    """Require one in-memory observation frame from an AnnData-family holder."""
    observations = holder.obs
    if not isinstance(observations, pd.DataFrame):
        raise TypeError("sample annotation requires an in-memory obs DataFrame")
    return observations


def _append_annotation_provenance(
    target: AnnData | MuData,
    provenance: SampleAnnotationProvenance,
) -> str:
    """Append typed provenance at the AnnData serialization boundary."""
    stored = read_namespace_text(target, _PROVENANCE_KEY)
    existing_json = "[]" if isinstance(stored, MissingNamespaceText) else stored
    decoded = to_json_compatible(json.loads(existing_json))
    if not isinstance(decoded, list):
        raise ValueError("stored obs annotation provenance must be a JSON list")
    records: list[JsonValue] = list(decoded)
    records.append(_provenance_record(provenance).model_dump(mode="json"))
    return json.dumps(records)


def _provenance_record(provenance: SampleAnnotationProvenance) -> _ObsAnnotationProvenance:
    """Select the stored provenance variant matching the annotation origin."""
    if isinstance(provenance.origin, AnnotationFileOrigin):
        return _FileObsAnnotationProvenance(
            source=str(provenance.origin.path),
            source_format=provenance.origin.path.suffix.lower().lstrip("."),
            match_on=provenance.match_on,
            key_field=provenance.key_field,
            obs_columns_added=provenance.columns_added,
            n_obs_matched=provenance.matched_observation_count,
        )
    return _InMemoryObsAnnotationProvenance(
        match_on=provenance.match_on,
        key_field=provenance.key_field,
        obs_columns_added=provenance.columns_added,
        n_obs_matched=provenance.matched_observation_count,
    )
