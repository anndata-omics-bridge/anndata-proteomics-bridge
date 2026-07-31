"""Pure peptide-to-FASTA validation and feature-edge calculations."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from prozor.annotate import annotate_peptides_streaming
from scipy.sparse import csr_matrix

from anndata_proteomics.annotation.var_fasta import leading_accession
from anndata_proteomics.fasta.annotation import (
    describe_sources,
    materialize_sources,
    parse_header_id,
    uniprot_proteinname,
)
from anndata_proteomics.fasta.config import (
    FastaConfig,
    FastaConfigAccumulator,
    ResolvedFastaConfig,
    matches_any,
)
from anndata_proteomics.fasta.parser import FastaSource, FastaSources, iter_fasta

_MATCH_COLUMNS = [
    "sequence",
    "protein_id",
    "proteinname",
    "start",
    "end",
    "length",
    "is_decoy",
    "is_contaminant",
]
_MAX_REPORTED = 5


@dataclass(frozen=True, slots=True)
class PeptideFastaMatchingConfig:
    """Complete scientific configuration for peptide-to-FASTA matching."""

    backend: str = "auto"
    identifiers: FastaConfig = field(default_factory=FastaConfig)
    il_equivalent: bool = False
    is_uniprot: bool = True


DEFAULT_PEPTIDE_FASTA_MATCHING_CONFIG = PeptideFastaMatchingConfig()


@dataclass(frozen=True, slots=True)
class PeptideSequenceCollection:
    """One named feature-level peptide sequence series."""

    name: str
    sequences: pd.Series


@dataclass(slots=True)
class FastaValidationResult:
    """FASTA matches and per-feature peptide-presence facts for one level."""

    summary: pd.DataFrame
    matches: pd.DataFrame
    normalized_sequences: pd.Series
    fasta_proteins: frozenset[str]
    fasta_sources: tuple[str, ...]
    n_features: int
    n_unique_sequences: int
    n_invalid_sequences: int
    n_matched_features: int
    n_unmatched_features: int
    requested_backend: str
    backend: str
    il_equivalent: bool
    is_uniprot: bool
    fasta_config: ResolvedFastaConfig
    unmatched_sequences: tuple[str, ...]

    @property
    def fraction_unmatched(self) -> float:
        """Return the fraction of features whose peptide was not found."""
        return self.n_unmatched_features / self.n_features if self.n_features else 0.0

    def sample_unmatched(self, k: int = _MAX_REPORTED) -> list[str]:
        """Return up to *k* distinct normalized unmatched sequences."""
        return list(self.unmatched_sequences[:k])


@dataclass(slots=True)
class ProteinAssignmentValidation:
    """Per-feature checks against reported leading proteins."""

    summary: pd.DataFrame


@dataclass(frozen=True, slots=True)
class ValidationTotals:
    """Aggregate matching counts across one or more validated levels."""

    n_levels: int
    n_features: int
    n_matched_features: int
    n_unique_patterns: int
    n_match_sites: int


def combined_validation_summary(
    matching: FastaValidationResult,
    proteins: ProteinAssignmentValidation,
) -> pd.DataFrame:
    """Join per-feature FASTA facts with reported-protein diagnostics."""
    return pd.concat([matching.summary, proteins.summary], axis="columns")


def validation_totals(results: Iterable[FastaValidationResult]) -> ValidationTotals:
    """Aggregate matching counts and distinct match sites across levels."""
    collected = list(results)
    unique_patterns = {
        sequence
        for result in collected
        for sequence in result.normalized_sequences
        if isinstance(sequence, str)
    }
    match_sites = pd.concat(
        [result.matches for result in collected],
        ignore_index=True,
    ).drop_duplicates()
    return ValidationTotals(
        n_levels=len(collected),
        n_features=sum(result.n_features for result in collected),
        n_matched_features=sum(result.n_matched_features for result in collected),
        n_unique_patterns=len(unique_patterns),
        n_match_sites=len(match_sites),
    )


@dataclass(frozen=True, slots=True)
class PeptideFeatureNode:
    """One normalized peptide feature at its global feature-axis position."""

    position: int
    sequence: str


@dataclass(frozen=True, slots=True)
class PeptideFeatureNodes:
    """Peptide feature nodes participating in one MuLink edge calculation."""

    nodes: tuple[PeptideFeatureNode, ...]
    total_nodes: int


@dataclass(frozen=True, slots=True)
class ProteinFeatureNodes:
    """Global protein-feature positions indexed by normalized accession."""

    positions_by_accession: dict[str, tuple[int, ...]]
    total_nodes: int


@dataclass(frozen=True, slots=True)
class PeptideProteinMatches:
    """FASTA protein accessions reached by each normalized peptide sequence."""

    accessions_by_sequence: dict[str, frozenset[str]]
    all_accessions: frozenset[str]


@dataclass(slots=True)
class FeatureMappingResult:
    """Calculated peptide-to-protein adjacency and representation diagnostics."""

    mapping: csr_matrix[np.int8]
    represented_accessions: frozenset[str]
    unrepresented_accessions: frozenset[str]


@dataclass(slots=True)
class _FastaScan:
    matches: pd.DataFrame
    fasta_proteins: frozenset[str]
    fasta_config: ResolvedFastaConfig
    requested_backend: str
    resolved_backend: str
    fasta_sources: tuple[str, ...]


def match_peptides_to_fasta(
    sequences: pd.Series,
    fasta_sources: FastaSources,
    config: PeptideFastaMatchingConfig = DEFAULT_PEPTIDE_FASTA_MATCHING_CONFIG,
) -> FastaValidationResult:
    """Match one feature-level peptide series to a FASTA database."""
    collection = PeptideSequenceCollection(name="features", sequences=sequences)
    return match_peptide_collections_to_fasta((collection,), fasta_sources, config)[collection.name]


def match_peptide_collections_to_fasta(
    collections: tuple[PeptideSequenceCollection, ...],
    fasta_sources: FastaSources,
    config: PeptideFastaMatchingConfig = DEFAULT_PEPTIDE_FASTA_MATCHING_CONFIG,
) -> dict[str, FastaValidationResult]:
    """Match multiple sequence collections with one shared FASTA scan."""
    if not collections:
        raise ValueError("at least one peptide sequence collection is required")
    names = [collection.name for collection in collections]
    if len(set(names)) != len(names):
        raise ValueError(f"peptide sequence collection names must be unique: {names}")

    normalized = {
        collection.name: normalize_peptide_sequences(
            collection.sequences,
            il_equivalent=config.il_equivalent,
        )
        for collection in collections
    }
    patterns = sorted(
        {
            sequence
            for values in normalized.values()
            for sequence in values
            if isinstance(sequence, str)
        }
    )
    scan = _scan_fasta(patterns, fasta_sources, config)
    per_pattern = _per_pattern_stats(scan.matches)
    return {
        name: _build_validation_result(values, scan, per_pattern, config)
        for name, values in normalized.items()
    }


def normalize_peptide_sequences(
    sequences: pd.Series,
    *,
    il_equivalent: bool,
) -> pd.Series:
    """Normalize valid peptide strings and mark all other values as missing."""
    normalized = pd.Series(
        pd.NA,
        index=pd.Index(sequences.index).astype(str),
        dtype="object",
    )
    for position, value in enumerate(sequences):
        if not isinstance(value, str) or not value.strip():
            continue
        sequence = value.strip().upper()
        normalized.iloc[position] = sequence.replace("I", "L") if il_equivalent else sequence
    return normalized


def validate_reported_proteins(
    reported_proteins: pd.Series,
    matches: FastaValidationResult,
) -> ProteinAssignmentValidation:
    """Validate reported leading proteins against peptide FASTA matches."""
    reported = reported_proteins.copy()
    reported.index = pd.Index(reported.index).astype(str)
    if not reported.index.equals(matches.normalized_sequences.index):
        raise ValueError("reported proteins and matched peptide features are not aligned")

    per_pattern = _per_pattern_stats(matches.matches)
    leading_in_fasta = pd.Series(
        pd.NA,
        index=matches.normalized_sequences.index,
        dtype="boolean",
    )
    peptide_in_leading = leading_in_fasta.copy()
    for position, (sequence, raw_protein) in enumerate(
        zip(
            matches.normalized_sequences,
            reported,
            strict=True,
        )
    ):
        if pd.isna(raw_protein) or not str(raw_protein).strip():
            continue
        leading = leading_accession(str(raw_protein), is_uniprot=matches.is_uniprot)
        proteins = per_pattern.get(sequence, (0, ()))[1] if isinstance(sequence, str) else ()
        leading_in_fasta.iloc[position] = leading in matches.fasta_proteins
        if isinstance(sequence, str):
            peptide_in_leading.iloc[position] = leading in proteins
    return ProteinAssignmentValidation(
        summary=_protein_assignment_frame(
            matches.normalized_sequences.index,
            leading_in_fasta,
            peptide_in_leading,
        )
    )


def unavailable_reported_protein_validation(
    feature_names: pd.Index,
) -> ProteinAssignmentValidation:
    """Return explicit unavailable values when no reported protein field exists."""
    missing = pd.Series(
        pd.NA,
        index=pd.Index(feature_names).astype(str),
        dtype="boolean",
    )
    return ProteinAssignmentValidation(
        summary=_protein_assignment_frame(feature_names, missing, missing)
    )


def peptide_protein_matches(matches: pd.DataFrame) -> PeptideProteinMatches:
    """Index matched FASTA protein accessions by peptide sequence."""
    by_sequence = {
        str(sequence): frozenset(group["proteinname"].astype(str))
        for sequence, group in matches.groupby("sequence", sort=False)
    }
    return PeptideProteinMatches(
        accessions_by_sequence=by_sequence,
        all_accessions=frozenset(matches["proteinname"].astype(str)),
    )


def build_feature_mapping(
    peptide_nodes: PeptideFeatureNodes,
    protein_nodes: ProteinFeatureNodes,
    matches: PeptideProteinMatches,
) -> FeatureMappingResult:
    """Calculate directed peptide-feature to protein-feature edges."""
    if peptide_nodes.total_nodes != protein_nodes.total_nodes:
        raise ValueError(
            "peptide and protein nodes describe different global feature-axis sizes: "
            f"{peptide_nodes.total_nodes} != {protein_nodes.total_nodes}"
        )
    rows: list[int] = []
    columns: list[int] = []
    represented: set[str] = set()
    for peptide in peptide_nodes.nodes:
        for accession in matches.accessions_by_sequence.get(peptide.sequence, frozenset()):
            protein_positions = protein_nodes.positions_by_accession.get(accession, ())
            rows.extend(peptide.position for _position in protein_positions)
            columns.extend(protein_positions)
            if protein_positions:
                represented.add(accession)
    mapping = csr_matrix(
        (
            np.ones(len(rows), dtype=np.int8),
            (
                np.asarray(rows, dtype=np.int64),
                np.asarray(columns, dtype=np.int64),
            ),
        ),
        shape=(peptide_nodes.total_nodes, peptide_nodes.total_nodes),
    )
    mapping.sum_duplicates()
    if mapping.nnz:
        mapping.data[:] = 1
    represented_accessions = frozenset(represented)
    return FeatureMappingResult(
        mapping=mapping,
        represented_accessions=represented_accessions,
        unrepresented_accessions=matches.all_accessions - represented_accessions,
    )


def _scan_fasta(
    patterns: list[str],
    fasta_sources: FastaSources,
    config: PeptideFastaMatchingConfig,
) -> _FastaScan:
    sources = materialize_sources(fasta_sources)
    source_descriptions = tuple(describe_sources(sources))
    accumulator = FastaConfigAccumulator(config.identifiers)
    fasta_proteins: set[str] = set()
    records = _protein_records(
        sources,
        accumulator,
        fasta_proteins,
        il_equivalent=config.il_equivalent,
        is_uniprot=config.is_uniprot,
    )
    if patterns:
        annotations = annotate_peptides_streaming(patterns, records, backend=config.backend)
    else:
        for _record in records:
            pass
        annotations = annotate_peptides_streaming([], (), backend=config.backend)

    effective_config = accumulator.resolve()
    rows = [
        (
            annotation.peptide,
            annotation.protein_id,
            uniprot_proteinname(annotation.protein_id)
            if config.is_uniprot
            else annotation.protein_id,
            annotation.start,
            annotation.end,
            annotation.length,
            matches_any(annotation.protein_id, effective_config.decoy.patterns),
            matches_any(annotation.protein_id, effective_config.contaminant.patterns),
        )
        for annotation in annotations
    ]
    frame = pd.DataFrame(rows, columns=_MATCH_COLUMNS)
    for column in ("start", "end", "length"):
        frame[column] = frame[column].astype("int64")
    for column in ("is_decoy", "is_contaminant"):
        frame[column] = frame[column].astype("bool")
    return _FastaScan(
        matches=frame,
        fasta_proteins=frozenset(fasta_proteins),
        fasta_config=effective_config,
        requested_backend=config.backend,
        resolved_backend=annotations.resolved_backend,
        fasta_sources=source_descriptions,
    )


def _protein_records(
    fasta_sources: list[FastaSource],
    accumulator: FastaConfigAccumulator,
    fasta_proteins: set[str],
    *,
    il_equivalent: bool,
    is_uniprot: bool,
) -> Iterable[tuple[str, str]]:
    for source in fasta_sources:
        for record in iter_fasta(source):
            parsed = parse_header_id(record.header)
            accumulator.observe(parsed.identifier)
            proteinname = (
                uniprot_proteinname(parsed.identifier) if is_uniprot else parsed.identifier
            )
            fasta_proteins.add(proteinname)
            sequence = record.sequence.upper()
            yield parsed.identifier, sequence.replace("I", "L") if il_equivalent else sequence


def _build_validation_result(
    normalized: pd.Series,
    scan: _FastaScan,
    per_pattern: dict[str, tuple[int, tuple[str, ...]]],
    config: PeptideFastaMatchingConfig,
) -> FastaValidationResult:
    summary = _build_match_summary(normalized, per_pattern)
    unique_sequences = {value for value in normalized if isinstance(value, str)}
    unmatched = tuple(
        sorted(sequence for sequence in unique_sequences if sequence not in per_pattern)
    )
    n_matched = int(summary["peptide_in_fasta"].sum())
    return FastaValidationResult(
        summary=summary,
        matches=scan.matches[scan.matches["sequence"].isin(unique_sequences)].copy(),
        normalized_sequences=normalized,
        fasta_proteins=scan.fasta_proteins,
        fasta_sources=scan.fasta_sources,
        n_features=len(summary),
        n_unique_sequences=len(unique_sequences),
        n_invalid_sequences=int(normalized.isna().sum()),
        n_matched_features=n_matched,
        n_unmatched_features=len(summary) - n_matched,
        requested_backend=scan.requested_backend,
        backend=scan.resolved_backend,
        il_equivalent=config.il_equivalent,
        is_uniprot=config.is_uniprot,
        fasta_config=scan.fasta_config,
        unmatched_sequences=unmatched,
    )


def _per_pattern_stats(
    matches: pd.DataFrame,
) -> dict[str, tuple[int, tuple[str, ...]]]:
    """Map sequence to total match sites and distinct matching accessions."""
    stats: dict[str, tuple[int, tuple[str, ...]]] = {}
    for sequence, group in matches.groupby("sequence", sort=False):
        proteins = tuple(sorted(set(group["proteinname"].astype(str))))
        stats[str(sequence)] = (len(group), proteins)
    return stats


def _build_match_summary(
    normalized: pd.Series,
    per_pattern: dict[str, tuple[int, tuple[str, ...]]],
) -> pd.DataFrame:
    peptide_in_fasta: list[bool] = []
    site_counts: list[int] = []
    protein_counts: list[int] = []
    protein_ids: list[str] = []
    for sequence in normalized:
        site_count, proteins = (
            per_pattern.get(sequence, (0, ())) if isinstance(sequence, str) else (0, ())
        )
        peptide_in_fasta.append(site_count > 0)
        site_counts.append(site_count)
        protein_counts.append(len(proteins))
        protein_ids.append(";".join(proteins))
    return pd.DataFrame(
        {
            "peptide_in_fasta": peptide_in_fasta,
            "fasta_match_site_count": site_counts,
            "fasta_matching_protein_count": protein_counts,
            "fasta_matching_protein_ids": protein_ids,
        },
        index=normalized.index,
    )


def _protein_assignment_frame(
    feature_names: pd.Index,
    leading_in_fasta: pd.Series,
    peptide_in_leading: pd.Series,
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "leading_protein_in_fasta": leading_in_fasta.array,
            "peptide_in_leading_protein": peptide_in_leading.array,
        },
        index=pd.Index(feature_names).astype(str),
    )
