"""Backend-independent orchestration for sample annotation."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd
from loguru import logger

from anndata_proteomics.annotation.loader import AnnotationOrigin, AnnotationTable
from anndata_proteomics.annotation.sample import (
    AnnotatedObservations,
    AnnotationDiagnostics,
    SampleAnnotationProvenance,
    annotate_observations,
    annotation_diagnostics,
    match_sample_annotation,
    observation_keys,
    sample_annotation_provenance,
)

_MAX_REPORTED = 5


@dataclass(frozen=True, slots=True)
class SampleAnnotationResult:
    """Annotated observation frames, diagnostics, and persistence provenance."""

    observations: tuple[AnnotatedObservations, ...]
    diagnostics: AnnotationDiagnostics
    provenance: SampleAnnotationProvenance


def run_sample_annotation(
    observation_frames: tuple[pd.DataFrame, ...],
    annotation: AnnotationTable,
    origin: AnnotationOrigin,
) -> SampleAnnotationResult:
    """Annotate one primary observation frame and its propagated level frames."""
    if not observation_frames:
        raise ValueError("sample annotation requires at least one observation frame")

    primary_keys = observation_keys(observation_frames[0], annotation.match_on)
    match = match_sample_annotation(primary_keys, annotation)
    annotated = tuple(
        annotate_observations(observations, match) for observations in observation_frames
    )
    diagnostics = annotation_diagnostics(match)
    provenance = sample_annotation_provenance(origin, match, annotated[-1])
    _log_diagnostics(diagnostics)
    logger.info(
        "annotated obs: +{} column(s) {}, {}/{} rows matched",
        len(provenance.columns_added),
        list(provenance.columns_added),
        provenance.matched_observation_count,
        match.observation_count,
    )
    return SampleAnnotationResult(
        observations=annotated,
        diagnostics=diagnostics,
        provenance=provenance,
    )


def _log_diagnostics(diagnostics: AnnotationDiagnostics) -> None:
    """Log partial observation and annotation-record mismatches."""
    if diagnostics.unmatched_observation_count:
        logger.warning(
            "{}/{} obs rows had no matching annotation record",
            diagnostics.unmatched_observation_count,
            diagnostics.observation_count,
        )
    if diagnostics.unmatched_record_keys:
        shown = diagnostics.unmatched_record_keys[:_MAX_REPORTED]
        tail = " …" if len(diagnostics.unmatched_record_keys) > _MAX_REPORTED else ""
        logger.warning(
            "{} annotation record(s) matched no obs row: {}{}",
            len(diagnostics.unmatched_record_keys),
            list(shown),
            tail,
        )
