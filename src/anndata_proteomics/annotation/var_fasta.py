"""Pure protein-level FASTA annotation calculations.

This module knows how protein-group values join to FASTA records and how the
resulting annotation frame is aligned.  It does not know where protein groups
are stored or how the result is persisted.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd
from loguru import logger

from anndata_proteomics.annotation._sanitize import sanitize_columns
from anndata_proteomics.fasta.annotation import (
    DEFAULT_FASTA_ANNOTATION_CONFIG,
    FastaAnnotationConfig,
    describe_sources,
    fasta_to_dataframe_with_config,
    materialize_sources,
    uniprot_proteinname,
)
from anndata_proteomics.fasta.config import ResolvedFastaConfig
from anndata_proteomics.fasta.parser import FastaSources

_MAX_REPORTED = 5
_JOIN_KEY = "proteinname"


@dataclass(frozen=True, slots=True)
class AllFastaColumns:
    """Store every column produced by FASTA annotation."""


@dataclass(frozen=True, slots=True)
class SelectedFastaColumns:
    """Store the named FASTA annotation columns that are available."""

    names: tuple[str, ...]


type FastaColumnSelection = AllFastaColumns | SelectedFastaColumns

ALL_FASTA_COLUMNS = AllFastaColumns()


@dataclass(frozen=True, slots=True)
class ProteinFastaAnnotationConfig:
    """Complete configuration for protein-to-FASTA annotation."""

    fasta: FastaAnnotationConfig = field(default_factory=lambda: DEFAULT_FASTA_ANNOTATION_CONFIG)
    columns: FastaColumnSelection = ALL_FASTA_COLUMNS


DEFAULT_PROTEIN_FASTA_ANNOTATION_CONFIG = ProteinFastaAnnotationConfig()


@dataclass(slots=True)
class ProteinFastaAnnotationResult:
    """A var-aligned protein annotation and its diagnostics."""

    frame: pd.DataFrame
    fasta_config: ResolvedFastaConfig
    fasta_sources: tuple[str, ...]
    join_keys: pd.Index
    n_matched: int
    unmatched_protein_groups: tuple[str, ...]
    unmatched_fasta_records: tuple[str, ...]


def annotate_proteins_from_fasta(
    protein_groups: pd.Series,
    fasta_sources: FastaSources,
    config: ProteinFastaAnnotationConfig = DEFAULT_PROTEIN_FASTA_ANNOTATION_CONFIG,
) -> ProteinFastaAnnotationResult:
    """Build a FASTA annotation aligned to one protein-feature series."""
    sources = materialize_sources(fasta_sources)
    source_descriptions = tuple(describe_sources(sources))
    annotation, resolved_config = fasta_to_dataframe_with_config(
        sources,
        config.fasta,
    )
    annotation = _index_by_join_key(annotation)
    keys = protein_group_join_keys(
        protein_groups,
        is_uniprot=config.fasta.is_uniprot,
    )
    in_table = keys.isin(annotation.index)
    n_matched = int(in_table.sum())
    if n_matched == 0:
        raise ValueError(
            "no var rows matched any FASTA record "
            "(leading accession of the protein group). "
            f"first var keys: {list(keys[:_MAX_REPORTED])}; "
            f"first FASTA proteinnames: {list(annotation.index[:_MAX_REPORTED])}"
        )

    frame = _align_annotation(
        protein_groups.index,
        keys,
        annotation,
        config.columns,
    )
    unmatched_groups = tuple(str(key) for key in keys[~in_table])
    key_set = set(keys)
    unmatched_records = tuple(str(key) for key in annotation.index if key not in key_set)
    _log_mismatches(len(keys), unmatched_groups, unmatched_records)
    return ProteinFastaAnnotationResult(
        frame=frame,
        fasta_config=resolved_config,
        fasta_sources=source_descriptions,
        join_keys=keys,
        n_matched=n_matched,
        unmatched_protein_groups=unmatched_groups,
        unmatched_fasta_records=unmatched_records,
    )


def protein_group_join_keys(
    protein_groups: pd.Series,
    *,
    is_uniprot: bool,
) -> pd.Index:
    """Return leading-accession join keys for protein-feature values."""
    return pd.Index(
        [leading_accession(str(value), is_uniprot=is_uniprot) for value in protein_groups],
        dtype="object",
    )


def leading_accession(group_value: str, *, is_uniprot: bool) -> str:
    """Return the first normalized accession represented by a protein group."""
    token = group_value.removeprefix("prt:").strip().split(";")[0].strip()
    return uniprot_proteinname(token) if is_uniprot else token


def protein_group_accessions(group_value: str, *, is_uniprot: bool) -> tuple[str, ...]:
    """Return all normalized accessions represented by a protein group."""
    raw = group_value.removeprefix("prt:").strip()
    tokens = (token.strip() for token in raw.split(";"))
    return tuple(uniprot_proteinname(token) if is_uniprot else token for token in tokens if token)


def _index_by_join_key(annotation: pd.DataFrame) -> pd.DataFrame:
    """Index by accession using target, curated, then input-order precedence."""
    if _JOIN_KEY not in annotation.columns:
        raise ValueError(f"FASTA frame is missing the {_JOIN_KEY!r} join column")
    indexed = annotation.copy()
    indexed["_input_order"] = range(len(indexed))
    indexed["_database_priority"] = indexed["fasta.id"].map(_database_priority)
    indexed = indexed.sort_values(
        ["is_decoy", "_database_priority", "_input_order"],
        kind="stable",
    )
    duplicated = indexed[_JOIN_KEY].duplicated()
    if duplicated.any():
        duplicates = sorted(indexed[_JOIN_KEY][duplicated].unique())[:_MAX_REPORTED]
        logger.warning(
            "{} duplicate {!r} value(s) in FASTA; using target > sp| > tr| > input order: {}",
            int(duplicated.sum()),
            _JOIN_KEY,
            duplicates,
        )
        indexed = indexed.loc[~duplicated]
    return (
        indexed.sort_values("_input_order", kind="stable")
        .drop(columns=["_input_order", "_database_priority"])
        .set_index(_JOIN_KEY)
    )


def _database_priority(fasta_id: str) -> int:
    """Return curated-UniProt, unreviewed-UniProt, then other priority."""
    if "sp|" in fasta_id:
        return 0
    if "tr|" in fasta_id:
        return 1
    return 2


def _align_annotation(
    feature_names: pd.Index,
    keys: pd.Index,
    annotation: pd.DataFrame,
    columns: FastaColumnSelection,
) -> pd.DataFrame:
    """Align selected FASTA annotation columns to the protein feature axis."""
    available = list(annotation.columns)
    if isinstance(columns, SelectedFastaColumns):
        selected = [name for name in columns.names if name in available]
    else:
        selected = available
    aligned = annotation.reindex(keys)[selected]
    for column in ("is_decoy", "is_contaminant"):
        if column in aligned.columns:
            aligned[column] = aligned[column].astype("boolean")
    for column in aligned.columns:
        dtype = aligned[column].dtype
        is_text = pd.api.types.is_object_dtype(dtype) or isinstance(dtype, pd.StringDtype)
        if is_text and aligned[column].isna().any():
            aligned[column] = pd.Categorical(aligned[column])
    aligned.columns = sanitize_columns(selected)
    aligned.index = pd.Index(feature_names).astype(str)
    return aligned


def _log_mismatches(
    n_protein_groups: int,
    unmatched_groups: tuple[str, ...],
    unmatched_records: tuple[str, ...],
) -> None:
    if unmatched_groups:
        logger.warning(
            "{}/{} var rows had no matching FASTA record",
            len(unmatched_groups),
            n_protein_groups,
        )
    if unmatched_records:
        shown = list(unmatched_records[:_MAX_REPORTED])
        tail = " …" if len(unmatched_records) > _MAX_REPORTED else ""
        logger.info(
            "{} FASTA record(s) matched no var row: {}{}",
            len(unmatched_records),
            shown,
            tail,
        )
