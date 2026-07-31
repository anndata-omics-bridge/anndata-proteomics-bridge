"""Pure sample-annotation matching, alignment, diagnostics, and provenance."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from anndata_proteomics.annotation._sanitize import sanitize_columns
from anndata_proteomics.annotation.loader import (
    AnnotationOrigin,
    AnnotationTable,
)

_MAX_REPORTED = 5


@dataclass(frozen=True, slots=True)
class AnnotationMatch:
    """Annotation records selected for one observation axis."""

    annotation: pd.DataFrame
    observation_keys: pd.Index
    observations_in_annotation: NDArray[np.bool_]
    match_on: str
    key_field: str

    @property
    def observation_count(self) -> int:
        """Return the number of observations considered for matching."""
        return len(self.observation_keys)

    @property
    def matched_observation_count(self) -> int:
        """Return the number of observations matched to an annotation record."""
        return int(self.observations_in_annotation.sum())


@dataclass(frozen=True, slots=True)
class AnnotatedObservations:
    """A newly annotated observation frame and the columns it gained."""

    frame: pd.DataFrame
    columns_added: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class AnnotationDiagnostics:
    """Mismatch counts and examples produced by sample-annotation matching."""

    observation_count: int
    unmatched_observation_count: int
    unmatched_record_keys: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class SampleAnnotationProvenance:
    """Typed provenance for one sample-annotation operation."""

    origin: AnnotationOrigin
    match_on: str
    key_field: str
    columns_added: tuple[str, ...]
    matched_observation_count: int


def observation_keys(observations: pd.DataFrame, match_on: str) -> pd.Index:
    """Return string join keys from an observation frame."""
    if match_on == "index":
        return pd.Index(observations.index, dtype="object").astype(str)
    if match_on not in observations.columns:
        raise ValueError(
            f"match_on column {match_on!r} not found in obs columns: {list(observations.columns)}"
        )
    return pd.Index(observations[match_on].astype(str))


def match_sample_annotation(
    keys: pd.Index,
    annotation: AnnotationTable,
) -> AnnotationMatch:
    """Select the first declared annotation identifier matching observation keys."""
    key_fields = _annotation_key_fields(annotation)
    selected = _build_annotation_frame(annotation, key_fields[0], key_fields)
    matched_key_field = key_fields[0]
    if not keys.isin(selected.index).any():
        for key_field in key_fields[1:]:
            selected = _build_annotation_frame(annotation, key_field, key_fields)
            matched_key_field = key_field
            if keys.isin(selected.index).any():
                break

    in_annotation = keys.isin(selected.index)
    matched_count = int(in_annotation.sum())
    if matched_count == 0:
        raise ValueError(
            "no obs rows matched any annotation record on "
            f"match_on={annotation.match_on!r} (key_field={matched_key_field!r}). "
            f"first obs keys: {list(keys[:_MAX_REPORTED])}; "
            f"first record keys: {list(selected.index[:_MAX_REPORTED])}"
        )
    return AnnotationMatch(
        annotation=selected,
        observation_keys=keys,
        observations_in_annotation=in_annotation,
        match_on=annotation.match_on,
        key_field=matched_key_field,
    )


def annotate_observations(
    observations: pd.DataFrame,
    match: AnnotationMatch,
) -> AnnotatedObservations:
    """Return a copy of observations with aligned annotation columns added."""
    overlap = [column for column in match.annotation.columns if column in observations.columns]
    if overlap:
        raise ValueError(f"annotation columns already present in obs: {overlap}")

    keys = observation_keys(observations, match.match_on)
    aligned = match.annotation.reindex(keys)
    annotated = observations.copy()
    for column in match.annotation.columns:
        annotated[column] = aligned[column].to_numpy()
    return AnnotatedObservations(
        frame=annotated,
        columns_added=tuple(match.annotation.columns),
    )


def annotation_diagnostics(match: AnnotationMatch) -> AnnotationDiagnostics:
    """Return observation and record mismatches for a completed match."""
    unmatched_observation_count = int((~match.observations_in_annotation).sum())
    unmatched_record_keys: tuple[str, ...] = ()
    if not match.key_field.endswith("_aliases"):
        key_set = set(match.observation_keys)
        unmatched_record_keys = tuple(
            str(key) for key in match.annotation.index if key not in key_set
        )
    return AnnotationDiagnostics(
        observation_count=match.observation_count,
        unmatched_observation_count=unmatched_observation_count,
        unmatched_record_keys=unmatched_record_keys,
    )


def sample_annotation_provenance(
    origin: AnnotationOrigin,
    match: AnnotationMatch,
    annotated: AnnotatedObservations,
) -> SampleAnnotationProvenance:
    """Build typed provenance for an annotation operation."""
    return SampleAnnotationProvenance(
        origin=origin,
        match_on=match.match_on,
        key_field=match.key_field,
        columns_added=annotated.columns_added,
        matched_observation_count=match.matched_observation_count,
    )


def _annotation_key_fields(annotation: AnnotationTable) -> tuple[str, ...]:
    """Return the primary sample identifier followed by declared exact aliases."""
    fields = [annotation.key_field]
    alias_field = f"{annotation.key_field}_alias"
    if alias_field in annotation.samples:
        fields.append(alias_field)
    aliases_field = f"{annotation.key_field}_aliases"
    if aliases_field in annotation.samples:
        fields.append(aliases_field)
    return tuple(fields)


def _build_annotation_frame(
    annotation: AnnotationTable,
    join_field: str,
    identifier_fields: tuple[str, ...],
) -> pd.DataFrame:
    """Index annotation rows by one identifier and sanitize output columns."""
    frame = annotation.samples.copy()
    if join_field.endswith("_aliases"):
        frame = frame.explode(join_field)
    frame = frame.loc[frame[join_field].notna()].copy()
    frame[join_field] = frame[join_field].astype(str)
    if frame[join_field].duplicated().any():
        duplicates = sorted(frame[join_field][frame[join_field].duplicated()].unique())
        raise ValueError(f"duplicate {join_field!r} values in annotation table: {duplicates}")
    frame.index = pd.Index(frame[join_field], name=join_field)
    frame = frame.drop(columns=list(identifier_fields))
    frame.columns = sanitize_columns(list(frame.columns))
    return frame
