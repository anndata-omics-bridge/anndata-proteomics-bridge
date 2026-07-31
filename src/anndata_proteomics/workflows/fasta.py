"""Backend-independent orchestration for FASTA validation and MuLink edges."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from loguru import logger
from numpy.typing import NDArray
from scipy.sparse import csr_matrix

from anndata_proteomics.annotation.mulink import (
    combined_peptide_protein_matches,
    feature_positions,
    merge_owned_feature_mapping,
    peptide_feature_nodes,
    protein_feature_nodes,
    require_feature_positions,
    target_row_mask,
)
from anndata_proteomics.annotation.validate_fasta import (
    FastaValidationResult,
    FeatureMappingResult,
    PeptideFastaMatchingConfig,
    PeptideSequenceCollection,
    build_feature_mapping,
    combined_validation_summary,
    match_peptide_collections_to_fasta,
    unavailable_reported_protein_validation,
    validate_reported_proteins,
    validation_totals,
)
from anndata_proteomics.annotation.var_fasta import (
    ALL_FASTA_COLUMNS,
    FastaColumnSelection,
    ProteinFastaAnnotationConfig,
)
from anndata_proteomics.fasta.annotation import (
    DEFAULT_CLEAVAGE,
    CleavageRule,
    FastaAnnotationConfig,
    ResolvedCleavage,
    custom_cleavage,
    resolve_cleavage_name,
)
from anndata_proteomics.fasta.config import FastaConfig, ResolvedFastaConfig
from anndata_proteomics.fasta.parser import FastaSources
from anndata_proteomics.params.model import Parameters

_DEFAULT_MIN_LENGTH = 7
_DEFAULT_MAX_LENGTH = 30


@dataclass(frozen=True, slots=True)
class StoredCleavageOrTrypsin:
    """Use stored search enzyme, with a visible trypsin fallback."""


@dataclass(frozen=True, slots=True)
class NamedCleavage:
    """Use one explicit enzyme name."""

    enzyme: str


@dataclass(frozen=True, slots=True)
class CustomCleavage:
    """Use one explicit cleavage rule."""

    rule: CleavageRule


type CleavageSelection = StoredCleavageOrTrypsin | NamedCleavage | CustomCleavage

STORED_CLEAVAGE_OR_TRYPSIN = StoredCleavageOrTrypsin()


@dataclass(frozen=True, slots=True)
class StoredMinimumLengthOrDefault:
    """Use stored minimum peptide length, otherwise APB's default."""


@dataclass(frozen=True, slots=True)
class MinimumPeptideLength:
    """Use one explicit minimum peptide length."""

    value: int


type MinimumLengthSelection = StoredMinimumLengthOrDefault | MinimumPeptideLength

STORED_MINIMUM_LENGTH_OR_DEFAULT = StoredMinimumLengthOrDefault()


@dataclass(frozen=True, slots=True)
class StoredMaximumLengthOrDefault:
    """Use stored maximum peptide length, otherwise APB's default."""


@dataclass(frozen=True, slots=True)
class MaximumPeptideLength:
    """Use one explicit maximum peptide length."""

    value: int


type MaximumLengthSelection = StoredMaximumLengthOrDefault | MaximumPeptideLength

STORED_MAXIMUM_LENGTH_OR_DEFAULT = StoredMaximumLengthOrDefault()


@dataclass(frozen=True, slots=True)
class StoredProteinSearchParameters:
    """Search parameters physically present on the protein container."""

    parameters: Parameters


@dataclass(frozen=True, slots=True)
class MissingProteinSearchParameters:
    """The protein container has no stored search parameters."""


type ProteinSearchParameterState = StoredProteinSearchParameters | MissingProteinSearchParameters

MISSING_PROTEIN_SEARCH_PARAMETERS = MissingProteinSearchParameters()


@dataclass(frozen=True, slots=True)
class ProteinAnnotationInput:
    """Physical protein values plus explicit FASTA calculation selections."""

    protein_groups: pd.Series
    match_on: str
    search_parameters: ProteinSearchParameterState
    identifiers: FastaConfig = field(default_factory=FastaConfig)
    is_uniprot: bool = True
    cleavage: CleavageSelection = STORED_CLEAVAGE_OR_TRYPSIN
    minimum_length: MinimumLengthSelection = STORED_MINIMUM_LENGTH_OR_DEFAULT
    maximum_length: MaximumLengthSelection = STORED_MAXIMUM_LENGTH_OR_DEFAULT
    include_sequence: bool = False
    columns: FastaColumnSelection = ALL_FASTA_COLUMNS


@dataclass(frozen=True, slots=True)
class ProteinAnnotationProvenance:
    """Resolved values needed to describe persisted protein FASTA annotation."""

    match_on: str
    cleavage: ResolvedCleavage
    minimum_length: int
    maximum_length: int


@dataclass(frozen=True, slots=True)
class ProteinAnnotationCalculation:
    """Resolved protein FASTA calculation and persistence provenance."""

    protein_groups: pd.Series
    config: ProteinFastaAnnotationConfig
    provenance: ProteinAnnotationProvenance


@dataclass(frozen=True, slots=True)
class PeptideLevelInput:
    """One feature-level peptide series without reported protein assignments."""

    name: str
    sequences: pd.Series


@dataclass(frozen=True, slots=True)
class PeptideLevelWithReportedProteinsInput:
    """One peptide series and its aligned reported protein assignments."""

    name: str
    sequences: pd.Series
    reported_protein_field: str
    reported_proteins: pd.Series


type PeptideValidationInput = PeptideLevelInput | PeptideLevelWithReportedProteinsInput


@dataclass(slots=True)
class PeptideLevelValidation:
    """One FASTA validation whose reported-protein diagnostics are unavailable."""

    name: str
    matching: FastaValidationResult
    summary: pd.DataFrame


@dataclass(slots=True)
class PeptideLevelValidationWithReportedProteins:
    """One FASTA validation with reported-protein diagnostics."""

    name: str
    matching: FastaValidationResult
    summary: pd.DataFrame
    reported_protein_field: str


type PeptideValidation = PeptideLevelValidation | PeptideLevelValidationWithReportedProteins


@dataclass(slots=True)
class PeptideValidationWorkflowResult:
    """Validations produced for one or more levels by one shared FASTA scan."""

    levels: dict[str, PeptideValidation]
    fasta_config: ResolvedFastaConfig


@dataclass(frozen=True, slots=True)
class FeatureMappingAxes:
    """Backend-neutral feature axes and protein groups used to build MuLink edges."""

    global_feature_names: pd.Index
    peptide_feature_names: dict[str, pd.Index]
    protein_groups: pd.Series


@dataclass(slots=True)
class MuLinkFeatureMappingResult:
    """Calculated MuLink edges plus rows owned by this FASTA enrichment."""

    feature_mapping: FeatureMappingResult
    target_rows: NDArray[np.bool_]


@dataclass(slots=True)
class FeatureMappingState:
    """Existing shared MuLink state and APB's recorded contribution."""

    existing: csr_matrix
    owned: csr_matrix


@dataclass(frozen=True, slots=True)
class FeatureMappingStorageStats:
    """Summary of a calculated MuLink feature-mapping update."""

    n_fasta_edges: int
    n_unrepresented_fasta_proteins: int
    protein_match_on: str


@dataclass(slots=True)
class FeatureMappingUpdate:
    """Purely calculated shared mapping, ownership state, and provenance stats."""

    merged: csr_matrix
    owned: csr_matrix
    stats: FeatureMappingStorageStats


def resolve_protein_annotation_input(
    inputs: ProteinAnnotationInput,
) -> ProteinAnnotationCalculation:
    """Resolve digestion selections without consulting a storage backend."""
    cleavage = _resolve_cleavage(inputs.cleavage, inputs.search_parameters)
    minimum = _resolve_minimum_length(inputs.minimum_length, inputs.search_parameters)
    maximum = _resolve_maximum_length(inputs.maximum_length, inputs.search_parameters)
    config = ProteinFastaAnnotationConfig(
        fasta=FastaAnnotationConfig(
            identifiers=inputs.identifiers,
            is_uniprot=inputs.is_uniprot,
            cleavage=cleavage,
            min_length=minimum,
            max_length=maximum,
            include_sequence=inputs.include_sequence,
        ),
        columns=inputs.columns,
    )
    return ProteinAnnotationCalculation(
        protein_groups=inputs.protein_groups,
        config=config,
        provenance=ProteinAnnotationProvenance(
            match_on=inputs.match_on,
            cleavage=cleavage,
            minimum_length=minimum,
            maximum_length=maximum,
        ),
    )


def validate_peptide_levels(
    inputs: tuple[PeptideValidationInput, ...],
    fasta_sources: FastaSources,
    config: PeptideFastaMatchingConfig,
) -> PeptideValidationWorkflowResult:
    """Validate one or more peptide levels with one shared FASTA scan."""
    collections = tuple(
        PeptideSequenceCollection(name=level.name, sequences=level.sequences) for level in inputs
    )
    matching = match_peptide_collections_to_fasta(collections, fasta_sources, config)
    levels = {
        level.name: _complete_level_validation(level, matching[level.name]) for level in inputs
    }
    first = next(iter(levels.values()))
    _log_validation_totals(levels)
    return PeptideValidationWorkflowResult(
        levels=levels,
        fasta_config=first.matching.fasta_config,
    )


def build_mulink_feature_mapping(
    axes: FeatureMappingAxes,
    validations: PeptideValidationWorkflowResult,
) -> MuLinkFeatureMappingResult:
    """Build FASTA peptide-to-protein edges on one explicit global feature axis."""
    positions = feature_positions(axes.global_feature_names)
    peptide_names = _peptide_feature_names(axes.peptide_feature_names, validations)
    peptide_rows = require_feature_positions(peptide_names, positions)
    is_uniprot = _require_shared_is_uniprot(validations)
    matching = [validation.matching for validation in validations.levels.values()]
    peptide_nodes = peptide_feature_nodes(
        (result.normalized_sequences for result in matching),
        positions,
    )
    protein_nodes = protein_feature_nodes(
        axes.protein_groups,
        positions,
        is_uniprot=is_uniprot,
    )
    matches = combined_peptide_protein_matches(result.matches for result in matching)
    return MuLinkFeatureMappingResult(
        feature_mapping=build_feature_mapping(peptide_nodes, protein_nodes, matches),
        target_rows=target_row_mask(peptide_rows, positions),
    )


def calculate_feature_mapping_update(
    state: FeatureMappingState,
    calculation: MuLinkFeatureMappingResult,
    protein_match_on: str,
) -> FeatureMappingUpdate:
    """Replace APB-owned edges while preserving unrelated shared MuLink state."""
    expected_shape = calculation.feature_mapping.mapping.shape
    if state.existing.shape != expected_shape:
        raise ValueError(
            f"existing feature mapping has shape {state.existing.shape}, expected {expected_shape}"
        )
    if state.owned.shape != expected_shape:
        raise ValueError(
            f"owned feature mapping has shape {state.owned.shape}, expected {expected_shape}"
        )
    merge = merge_owned_feature_mapping(
        state.existing,
        state.owned,
        calculation.feature_mapping.mapping,
        calculation.target_rows,
    )
    return FeatureMappingUpdate(
        merged=merge.merged,
        owned=merge.owned,
        stats=FeatureMappingStorageStats(
            n_fasta_edges=calculation.feature_mapping.mapping.nnz,
            n_unrepresented_fasta_proteins=len(
                calculation.feature_mapping.unrepresented_accessions
            ),
            protein_match_on=protein_match_on,
        ),
    )


def _complete_level_validation(
    inputs: PeptideValidationInput,
    matching: FastaValidationResult,
) -> PeptideValidation:
    """Add reported-protein diagnostics to one peptide matching result."""
    if isinstance(inputs, PeptideLevelWithReportedProteinsInput):
        reported = validate_reported_proteins(inputs.reported_proteins, matching)
        return PeptideLevelValidationWithReportedProteins(
            name=inputs.name,
            matching=matching,
            summary=combined_validation_summary(matching, reported),
            reported_protein_field=inputs.reported_protein_field,
        )
    unavailable = unavailable_reported_protein_validation(matching.summary.index)
    return PeptideLevelValidation(
        name=inputs.name,
        matching=matching,
        summary=combined_validation_summary(matching, unavailable),
    )


def _peptide_feature_names(
    feature_names: dict[str, pd.Index],
    validations: PeptideValidationWorkflowResult,
) -> list[str]:
    """Return peptide feature names in validation-level order."""
    expected = set(validations.levels)
    observed = set(feature_names)
    if observed != expected:
        raise ValueError(
            "peptide feature axes and FASTA validations name different levels: "
            f"axes={sorted(observed)}, validations={sorted(expected)}"
        )
    return [
        str(feature_name) for level in validations.levels for feature_name in feature_names[level]
    ]


def _require_shared_is_uniprot(validations: PeptideValidationWorkflowResult) -> bool:
    """Require one accession interpretation across validation results."""
    values = {validation.matching.is_uniprot for validation in validations.levels.values()}
    if len(values) != 1:
        raise ValueError("FASTA validation results disagree on accession interpretation")
    return values.pop()


def _resolve_cleavage(
    selection: CleavageSelection,
    search_parameters: ProteinSearchParameterState,
) -> ResolvedCleavage:
    if isinstance(selection, NamedCleavage):
        return resolve_cleavage_name(selection.enzyme)
    if isinstance(selection, CustomCleavage):
        return custom_cleavage(selection.rule)
    if (
        isinstance(search_parameters, StoredProteinSearchParameters)
        and search_parameters.parameters.enzyme is not None
    ):
        return resolve_cleavage_name(search_parameters.parameters.enzyme)
    logger.warning(
        "no enzyme in search parameters and no cleavage override; "
        "using Trypsin for the peptide count"
    )
    return DEFAULT_CLEAVAGE


def _resolve_minimum_length(
    selection: MinimumLengthSelection,
    search_parameters: ProteinSearchParameterState,
) -> int:
    if isinstance(selection, MinimumPeptideLength):
        return selection.value
    if (
        isinstance(search_parameters, StoredProteinSearchParameters)
        and search_parameters.parameters.min_peptide_length is not None
    ):
        return int(search_parameters.parameters.min_peptide_length)
    return _DEFAULT_MIN_LENGTH


def _resolve_maximum_length(
    selection: MaximumLengthSelection,
    search_parameters: ProteinSearchParameterState,
) -> int:
    if isinstance(selection, MaximumPeptideLength):
        return selection.value
    if (
        isinstance(search_parameters, StoredProteinSearchParameters)
        and search_parameters.parameters.max_peptide_length is not None
    ):
        return int(search_parameters.parameters.max_peptide_length)
    return _DEFAULT_MAX_LENGTH


def _log_validation_totals(validations: dict[str, PeptideValidation]) -> None:
    """Log aggregate matching totals across validated levels."""
    totals = validation_totals(validation.matching for validation in validations.values())
    logger.info(
        "FASTA validation: {}/{} features matched across {} modality/modalities; "
        "{} unique peptide patterns, {} match sites",
        totals.n_matched_features,
        totals.n_features,
        totals.n_levels,
        totals.n_unique_patterns,
        totals.n_match_sites,
    )
