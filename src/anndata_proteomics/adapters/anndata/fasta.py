"""AnnData and MuData adapter for the pure FASTA calculations."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Literal

import pandas as pd
from anndata import AnnData
from loguru import logger
from mudata import MuData
from pydantic import BaseModel, ConfigDict
from scipy.sparse import csr_matrix

from anndata_proteomics.adapters.anndata.namespace import (
    MissingNamespaceText,
    has_namespace_key,
    read_namespace_text,
    update_namespace,
)
from anndata_proteomics.adapters.anndata.params import (
    has_search_parameters,
    require_search_parameters,
)
from anndata_proteomics.adapters.anndata.rules import (
    has_stored_column_role,
    require_stored_column_role,
)
from anndata_proteomics.annotation._sanitize import sanitize_columns
from anndata_proteomics.annotation.mulink import empty_feature_mapping
from anndata_proteomics.annotation.validate_fasta import (
    DEFAULT_PEPTIDE_FASTA_MATCHING_CONFIG,
    PeptideFastaMatchingConfig,
)
from anndata_proteomics.annotation.var_fasta import (
    ALL_FASTA_COLUMNS,
    FastaColumnSelection,
    ProteinFastaAnnotationResult,
)
from anndata_proteomics.fasta.annotation import (
    ResolvedCleavage,
)
from anndata_proteomics.fasta.config import FastaConfig, ResolvedFastaConfig
from anndata_proteomics.rules.schema import PEPTIDE_LEVELS
from anndata_proteomics.serialization import to_json_compatible
from anndata_proteomics.workflows.fasta import (
    MISSING_PROTEIN_SEARCH_PARAMETERS,
    STORED_CLEAVAGE_OR_TRYPSIN,
    STORED_MAXIMUM_LENGTH_OR_DEFAULT,
    STORED_MINIMUM_LENGTH_OR_DEFAULT,
    CleavageSelection,
    FeatureMappingAxes,
    FeatureMappingState,
    FeatureMappingStorageStats,
    FeatureMappingUpdate,
    MaximumLengthSelection,
    MinimumLengthSelection,
    PeptideLevelInput,
    PeptideLevelValidationWithReportedProteins,
    PeptideLevelWithReportedProteinsInput,
    PeptideValidation,
    PeptideValidationInput,
    ProteinAnnotationInput,
    ProteinAnnotationProvenance,
    StoredProteinSearchParameters,
)

_CONFIG_KEY = "fasta_config"
_LEVEL_KEY = "quantification_level"
_VAR_PROVENANCE_KEY = "var_annotations_json"
_SCHEMA_VERSION: Literal["0.2"] = "0.2"
_DEFAULT_SEQUENCE_FIELD = "ProForma_peptide"
_PROTEIN_ANNOTATION_KEY = "fasta"
_VALIDATION_KEY = "fasta_validation"
_FEATURE_MAPPING_KEY = "feature_mapping"
_OWNED_FEATURE_MAPPING_KEY = "_apb_fasta_feature_mapping_contribution"


class _VarAnnotationProvenance(BaseModel):
    """Typed base for one serialized ``var`` enrichment record."""

    model_config = ConfigDict(frozen=True)


class _ProteinAnnotationStorageProvenance(_VarAnnotationProvenance):
    """Stored provenance for protein-level FASTA annotation."""

    schema_version: Literal["0.2"] = _SCHEMA_VERSION
    source: Literal["fasta"] = "fasta"
    destination: str
    fasta_sources: tuple[str, ...]
    fasta_config: ResolvedFastaConfig
    match_on: str
    columns: tuple[str, ...]
    n_var_matched: int
    cleavage_enzyme: str
    min_peptide_length: int
    max_peptide_length: int


class _PeptideValidationStorageProvenance(_VarAnnotationProvenance):
    """Stored peptide-validation provenance without MuLink edges."""

    schema_version: Literal["0.2"] = _SCHEMA_VERSION
    source: Literal["fasta_validation"] = "fasta_validation"
    destination: str
    feature_mapping: Literal[None] = None
    feature_mapping_ownership: Literal[None] = None
    fasta_sources: tuple[str, ...]
    fasta_config: ResolvedFastaConfig
    requested_backend: str
    backend: str
    sequence_field: str
    leading_protein_field: str | None
    il_equivalent: bool
    n_features: int
    n_unique_sequences: int
    n_invalid_sequences: int
    n_matched_features: int
    n_unmatched_features: int
    n_feature_mapping_edges: Literal[0] = 0
    n_unrepresented_fasta_proteins: Literal[0] = 0
    protein_match_on: Literal[None] = None


class _PeptideValidationWithMappingStorageProvenance(_VarAnnotationProvenance):
    """Stored peptide-validation provenance with MuLink edges."""

    schema_version: Literal["0.2"] = _SCHEMA_VERSION
    source: Literal["fasta_validation"] = "fasta_validation"
    destination: str
    feature_mapping: str
    feature_mapping_ownership: str
    fasta_sources: tuple[str, ...]
    fasta_config: ResolvedFastaConfig
    requested_backend: str
    backend: str
    sequence_field: str
    leading_protein_field: str | None
    il_equivalent: bool
    n_features: int
    n_unique_sequences: int
    n_invalid_sequences: int
    n_matched_features: int
    n_unmatched_features: int
    n_feature_mapping_edges: int
    n_unrepresented_fasta_proteins: int
    protein_match_on: str


@dataclass(frozen=True, slots=True)
class StoredFastaAccessions:
    """Read the FASTA accession column declared by the stored APB rule."""


@dataclass(frozen=True, slots=True)
class FastaAccessionsColumn:
    """Read FASTA accessions from one explicit var column."""

    name: str


@dataclass(frozen=True, slots=True)
class FastaAccessionsIndex:
    """Read FASTA accessions from var names."""


type FastaAccessionsSelection = StoredFastaAccessions | FastaAccessionsColumn | FastaAccessionsIndex

STORED_FASTA_ACCESSIONS = StoredFastaAccessions()
FASTA_ACCESSIONS_INDEX = FastaAccessionsIndex()


@dataclass(frozen=True, slots=True)
class StoredReportedProteins:
    """Validate reported proteins when the stored APB rule declares their column."""


@dataclass(frozen=True, slots=True)
class ReportedProteinsColumn:
    """Validate reported proteins from one explicit var column."""

    name: str


type ReportedProteinsSelection = StoredReportedProteins | ReportedProteinsColumn

STORED_REPORTED_PROTEINS = StoredReportedProteins()


@dataclass(frozen=True, slots=True)
class AnnDataProteinFastaConfig:
    """Storage selections plus scientific settings for protein annotation."""

    match_on: FastaAccessionsSelection = STORED_FASTA_ACCESSIONS
    identifiers: FastaConfig = field(default_factory=FastaConfig)
    is_uniprot: bool = True
    cleavage: CleavageSelection = STORED_CLEAVAGE_OR_TRYPSIN
    minimum_length: MinimumLengthSelection = STORED_MINIMUM_LENGTH_OR_DEFAULT
    maximum_length: MaximumLengthSelection = STORED_MAXIMUM_LENGTH_OR_DEFAULT
    include_sequence: bool = False
    columns: FastaColumnSelection = ALL_FASTA_COLUMNS


DEFAULT_ANNDATA_PROTEIN_FASTA_CONFIG = AnnDataProteinFastaConfig()


@dataclass(frozen=True, slots=True)
class AnnDataPeptideFastaConfig:
    """Storage selections plus scientific settings for peptide validation."""

    sequence_field: str = _DEFAULT_SEQUENCE_FIELD
    matching: PeptideFastaMatchingConfig = DEFAULT_PEPTIDE_FASTA_MATCHING_CONFIG
    reported_proteins: ReportedProteinsSelection = STORED_REPORTED_PROTEINS


DEFAULT_ANNDATA_PEPTIDE_FASTA_CONFIG = AnnDataPeptideFastaConfig()


@dataclass(frozen=True, slots=True)
class MuDataPeptideFastaConfig:
    """Peptide validation settings plus protein-node accession selection."""

    validation: AnnDataPeptideFastaConfig = DEFAULT_ANNDATA_PEPTIDE_FASTA_CONFIG
    protein_match_on: FastaAccessionsSelection = STORED_FASTA_ACCESSIONS


DEFAULT_MUDATA_PEPTIDE_FASTA_CONFIG = MuDataPeptideFastaConfig()


@dataclass(frozen=True, slots=True)
class FeatureMappingExtraction:
    """Physical feature axes plus the selected protein storage field."""

    axes: FeatureMappingAxes
    protein_match_on: str


def read_protein_annotation_input(
    target: AnnData,
    config: AnnDataProteinFastaConfig = DEFAULT_ANNDATA_PROTEIN_FASTA_CONFIG,
) -> ProteinAnnotationInput:
    """Extract physical protein values and stored parameter state."""
    _require_protein_level(target)
    if _PROTEIN_ANNOTATION_KEY in target.varm:
        raise ValueError(
            f"varm[{_PROTEIN_ANNOTATION_KEY!r}] already present on the protein layer; "
            "delete it before re-annotating"
        )
    match_on, protein_groups = _read_fasta_accessions(target, config.match_on)
    search_parameters = (
        StoredProteinSearchParameters(require_search_parameters(target))
        if has_search_parameters(target)
        else MISSING_PROTEIN_SEARCH_PARAMETERS
    )
    return ProteinAnnotationInput(
        protein_groups=protein_groups,
        match_on=match_on,
        search_parameters=search_parameters,
        identifiers=config.identifiers,
        is_uniprot=config.is_uniprot,
        cleavage=config.cleavage,
        minimum_length=config.minimum_length,
        maximum_length=config.maximum_length,
        include_sequence=config.include_sequence,
        columns=config.columns,
    )


def require_protein_mudata_target(
    target: MuData,
) -> AnnData:
    """Return the concrete protein modality or raise precisely."""
    if "protein" not in target.mod:
        raise ValueError(f"MuData has no 'protein' modality; modalities: {list(target.mod)}")
    protein = target.mod["protein"]
    if not isinstance(protein, AnnData):
        raise TypeError("MuData 'protein' modality is not an AnnData")
    return protein


def write_protein_annotation(
    target: AnnData,
    result: ProteinFastaAnnotationResult,
    provenance: ProteinAnnotationProvenance,
) -> None:
    """Persist one calculated protein FASTA annotation."""
    target.varm[_PROTEIN_ANNOTATION_KEY] = result.frame
    _store_protein_annotation_provenance(
        target,
        result,
        provenance.match_on,
        provenance.cleavage,
        provenance.minimum_length,
        provenance.maximum_length,
    )
    logger.info(
        "stored protein annotation in varm[{!r}]: {} column(s) {}, {}/{} rows matched (enzyme={})",
        _PROTEIN_ANNOTATION_KEY,
        len(result.frame.columns),
        list(result.frame.columns),
        result.n_matched,
        len(result.frame),
        provenance.cleavage.enzyme,
    )


def read_peptide_anndata_input(
    target: AnnData,
    config: AnnDataPeptideFastaConfig = DEFAULT_ANNDATA_PEPTIDE_FASTA_CONFIG,
) -> PeptideValidationInput:
    """Extract one standalone peptide-derived level for FASTA validation."""
    return _read_peptide_input(
        _require_peptide_level(target),
        target,
        config,
    )


def read_peptide_mudata_modality_input(
    target: MuData,
    modality: str,
    config: MuDataPeptideFastaConfig = DEFAULT_MUDATA_PEPTIDE_FASTA_CONFIG,
) -> PeptideValidationInput:
    """Extract one selected peptide-derived MuData modality."""
    return _read_peptide_input(
        modality,
        require_peptide_mudata_target(target, modality),
        config.validation,
    )


def read_peptide_mudata_inputs(
    target: MuData,
    config: MuDataPeptideFastaConfig = DEFAULT_MUDATA_PEPTIDE_FASTA_CONFIG,
) -> tuple[PeptideValidationInput, ...]:
    """Extract every peptide-derived MuData modality for one shared scan."""
    extracted = tuple(
        _read_peptide_input(name, modality, config.validation)
        for name, modality in target.mod.items()
        if isinstance(modality, AnnData) and _is_peptide_level(modality)
    )
    if not extracted:
        raise ValueError("object has no peptide-derived modality to validate")
    return extracted


def require_peptide_mudata_target(target: MuData, modality: str) -> AnnData:
    """Return one peptide-derived modality or raise precisely."""
    if modality not in target.mod:
        raise ValueError(f"modality {modality!r} not in MuData; modalities: {list(target.mod)}")
    selected = target.mod[modality]
    if not isinstance(selected, AnnData):
        raise TypeError(f"modality {modality!r} is not an AnnData")
    if not _is_peptide_level(selected):
        raise ValueError(f"modality {modality!r} is not peptide-derived")
    return selected


def has_protein_mudata_target(target: MuData) -> bool:
    """Return whether MuData declares a protein modality."""
    return "protein" in target.mod


def read_feature_mapping_input(
    target: MuData,
    peptide_levels: tuple[str, ...],
    protein_selection: FastaAccessionsSelection,
) -> FeatureMappingExtraction:
    """Extract the global axis and protein groups needed for MuLink edges."""
    protein = require_protein_mudata_target(target)
    resolved_match_on, protein_groups = _read_fasta_accessions(protein, protein_selection)
    peptide_feature_names = {
        level: pd.Index(require_peptide_mudata_target(target, level).var_names).astype(str)
        for level in peptide_levels
    }
    return FeatureMappingExtraction(
        axes=FeatureMappingAxes(
            global_feature_names=pd.Index(target.var_names).astype(str),
            peptide_feature_names=peptide_feature_names,
            protein_groups=protein_groups,
        ),
        protein_match_on=resolved_match_on,
    )


def read_feature_mapping_state(
    target: MuData,
    expected_shape: tuple[int, int],
) -> FeatureMappingState:
    """Extract existing shared MuLink state without calculating an update."""
    return FeatureMappingState(
        existing=_stored_feature_mapping(target, _FEATURE_MAPPING_KEY, expected_shape),
        owned=_stored_feature_mapping(target, _OWNED_FEATURE_MAPPING_KEY, expected_shape),
    )


def write_feature_mapping(target: MuData, update: FeatureMappingUpdate) -> None:
    """Persist a previously calculated MuLink mapping and ownership state."""
    target.varp[_FEATURE_MAPPING_KEY] = update.merged
    target.varp[_OWNED_FEATURE_MAPPING_KEY] = update.owned


def write_peptide_validation(
    target: AnnData,
    result: PeptideValidation,
    sequence_field: str,
) -> None:
    """Persist one peptide FASTA validation without MuLink edges."""
    stored = result.summary.copy()
    stored.columns = sanitize_columns(list(stored.columns))
    target.varm[_VALIDATION_KEY] = stored
    _store_validation_provenance_without_feature_mapping(target, result, sequence_field)


def write_peptide_validation_with_feature_mapping(
    target: AnnData,
    result: PeptideValidation,
    sequence_field: str,
    feature_mapping: FeatureMappingStorageStats,
) -> None:
    """Persist one peptide FASTA validation with stored MuLink-edge provenance."""
    stored = result.summary.copy()
    stored.columns = sanitize_columns(list(stored.columns))
    target.varm[_VALIDATION_KEY] = stored
    _store_validation_provenance_with_feature_mapping(
        target,
        result,
        sequence_field,
        feature_mapping,
    )


def has_fasta_config(target: AnnData | MuData) -> bool:
    """Return whether a container stores a resolved FASTA configuration."""
    return has_namespace_key(target, _CONFIG_KEY)


def require_fasta_config(target: AnnData | MuData) -> ResolvedFastaConfig:
    """Read the stored resolved FASTA configuration or raise precisely."""
    payload = read_namespace_text(target, _CONFIG_KEY)
    if isinstance(payload, MissingNamespaceText):
        raise ValueError("container has no stored FASTA configuration")
    return ResolvedFastaConfig.model_validate_json(payload)


def write_anndata_fasta_config(target: AnnData, config: ResolvedFastaConfig) -> None:
    """Store a resolved FASTA configuration on AnnData."""
    _write_fasta_config(target, config)


def write_mudata_fasta_config(target: MuData, config: ResolvedFastaConfig) -> None:
    """Store a resolved FASTA configuration on MuData."""
    _write_fasta_config(target, config)


def _read_peptide_input(
    name: str,
    target: AnnData,
    config: AnnDataPeptideFastaConfig,
) -> PeptideValidationInput:
    if config.sequence_field not in target.var.columns:
        raise ValueError(
            f"sequence_field {config.sequence_field!r} not in var columns: "
            f"{list(target.var.columns)}"
        )
    sequences = target.var[config.sequence_field].copy()
    sequences.index = pd.Index(target.var_names).astype(str)
    if isinstance(config.reported_proteins, ReportedProteinsColumn):
        column = config.reported_proteins.name
        if column == "index":
            raise ValueError("leading_protein_field must name a var column, not 'index'")
        if column not in target.var.columns:
            raise ValueError(
                f"leading_protein_field {column!r} not in var columns: {list(target.var.columns)}"
            )
    else:
        if not has_stored_column_role(target, "fasta_accessions"):
            return PeptideLevelInput(name=name, sequences=sequences)
        column = require_stored_column_role(target, "fasta_accessions")
    values = target.var[column].copy()
    values.index = pd.Index(target.var_names).astype(str)
    return PeptideLevelWithReportedProteinsInput(
        name=name,
        sequences=sequences,
        reported_protein_field=column,
        reported_proteins=values,
    )


def _read_fasta_accessions(
    target: AnnData,
    selection: FastaAccessionsSelection,
) -> tuple[str, pd.Series]:
    if isinstance(selection, FastaAccessionsIndex):
        return "index", pd.Series(target.var_names, index=target.var_names, dtype="object")
    if isinstance(selection, FastaAccessionsColumn):
        column = selection.name
    else:
        if not has_stored_column_role(target, "fasta_accessions"):
            raise ValueError(
                "match_on was not provided and the stored APB rule does not declare "
                "column_roles.fasta_accessions"
            )
        column = require_stored_column_role(target, "fasta_accessions")
    if column not in target.var.columns:
        raise ValueError(f"match_on {column!r} not in var columns: {list(target.var.columns)}")
    values = target.var[column].copy()
    values.index = pd.Index(target.var_names).astype(str)
    return column, values


def _require_protein_level(target: AnnData) -> None:
    level = _require_namespace_level(target)
    if level != "protein":
        raise ValueError(
            "FASTA var-annotation applies to the protein layer only; got "
            f"quantification_level={level!r}"
        )


def _require_peptide_level(target: AnnData) -> str:
    level = _require_namespace_level(target)
    if level not in PEPTIDE_LEVELS:
        raise ValueError(
            "FASTA validation applies to peptide-derived layers "
            f"(ion/fragment/peptidoform/peptide), got {level!r}"
        )
    return level


def _is_peptide_level(target: AnnData) -> bool:
    return _has_namespace_level(target) and _require_namespace_level(target) in PEPTIDE_LEVELS


def _has_namespace_level(target: AnnData) -> bool:
    """Return whether the APB namespace contains a text quantification level."""
    return not isinstance(read_namespace_text(target, _LEVEL_KEY), MissingNamespaceText)


def _require_namespace_level(target: AnnData) -> str:
    """Return the stored quantification level or raise a precise boundary error."""
    level = read_namespace_text(target, _LEVEL_KEY)
    if isinstance(level, MissingNamespaceText):
        raise ValueError("container has no text quantification_level in the APB namespace")
    return level


def _write_fasta_config(
    target: AnnData | MuData,
    config: ResolvedFastaConfig,
) -> None:
    update_namespace(target, {_CONFIG_KEY: config.model_dump_json()})


def _store_protein_annotation_provenance(
    target: AnnData,
    result: ProteinFastaAnnotationResult,
    match_on: str,
    cleavage: ResolvedCleavage,
    min_length: int,
    max_length: int,
) -> None:
    entry = _ProteinAnnotationStorageProvenance(
        destination=f"varm[{_PROTEIN_ANNOTATION_KEY!r}]",
        fasta_sources=result.fasta_sources,
        fasta_config=result.fasta_config,
        match_on=match_on,
        columns=tuple(str(column) for column in result.frame.columns),
        n_var_matched=result.n_matched,
        cleavage_enzyme=cleavage.enzyme,
        min_peptide_length=min_length,
        max_peptide_length=max_length,
    )
    _append_var_provenance(target, entry)


def _store_validation_provenance_without_feature_mapping(
    target: AnnData,
    result: PeptideValidation,
    sequence_field: str,
) -> None:
    validation = result.matching
    leading_protein_field = (
        result.reported_protein_field
        if isinstance(result, PeptideLevelValidationWithReportedProteins)
        else None
    )
    entry = _PeptideValidationStorageProvenance(
        destination=f"varm[{_VALIDATION_KEY!r}]",
        fasta_sources=validation.fasta_sources,
        fasta_config=validation.fasta_config,
        requested_backend=validation.requested_backend,
        backend=validation.backend,
        sequence_field=sequence_field,
        leading_protein_field=leading_protein_field,
        il_equivalent=validation.il_equivalent,
        n_features=validation.n_features,
        n_unique_sequences=validation.n_unique_sequences,
        n_invalid_sequences=validation.n_invalid_sequences,
        n_matched_features=validation.n_matched_features,
        n_unmatched_features=validation.n_unmatched_features,
    )
    _append_var_provenance(target, entry)


def _store_validation_provenance_with_feature_mapping(
    target: AnnData,
    result: PeptideValidation,
    sequence_field: str,
    feature_mapping: FeatureMappingStorageStats,
) -> None:
    validation = result.matching
    leading_protein_field = (
        result.reported_protein_field
        if isinstance(result, PeptideLevelValidationWithReportedProteins)
        else None
    )
    entry = _PeptideValidationWithMappingStorageProvenance(
        destination=f"varm[{_VALIDATION_KEY!r}]",
        feature_mapping=f"varp[{_FEATURE_MAPPING_KEY!r}]",
        feature_mapping_ownership=f"varp[{_OWNED_FEATURE_MAPPING_KEY!r}]",
        fasta_sources=validation.fasta_sources,
        fasta_config=validation.fasta_config,
        requested_backend=validation.requested_backend,
        backend=validation.backend,
        sequence_field=sequence_field,
        leading_protein_field=leading_protein_field,
        il_equivalent=validation.il_equivalent,
        n_features=validation.n_features,
        n_unique_sequences=validation.n_unique_sequences,
        n_invalid_sequences=validation.n_invalid_sequences,
        n_matched_features=validation.n_matched_features,
        n_unmatched_features=validation.n_unmatched_features,
        n_feature_mapping_edges=feature_mapping.n_fasta_edges,
        n_unrepresented_fasta_proteins=(feature_mapping.n_unrepresented_fasta_proteins),
        protein_match_on=feature_mapping.protein_match_on,
    )
    _append_var_provenance(target, entry)


def _append_var_provenance(
    target: AnnData,
    entry: _VarAnnotationProvenance,
) -> None:
    stored = read_namespace_text(target, _VAR_PROVENANCE_KEY)
    existing_payload = "[]" if isinstance(stored, MissingNamespaceText) else stored
    existing = to_json_compatible(json.loads(existing_payload))
    if not isinstance(existing, list) or any(not isinstance(item, dict) for item in existing):
        raise TypeError("stored var annotation provenance must be a JSON array of objects")
    serialized = to_json_compatible(entry.model_dump(mode="json"))
    if not isinstance(serialized, dict):
        raise TypeError("var annotation provenance did not serialize to an object")
    existing.append(serialized)
    update_namespace(target, {_VAR_PROVENANCE_KEY: json.dumps(existing)})


def _stored_feature_mapping(
    target: MuData,
    key: str,
    expected_shape: tuple[int, int],
) -> csr_matrix:
    existing = target.varp.get(key)
    mapping = empty_feature_mapping(expected_shape) if existing is None else csr_matrix(existing)
    if mapping.shape != expected_shape:
        raise ValueError(
            f"existing varp[{key!r}] has shape {mapping.shape}, expected {expected_shape}"
        )
    return mapping
