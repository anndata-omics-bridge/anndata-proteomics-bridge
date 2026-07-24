"""Validate peptide identifications against FASTA with Aho--Corasick.

Validation is annotation only.  It never removes quantified features, including
features reported for decoys or contaminants.  Per-feature facts are stored in
``varm['fasta_validation']``.  For MuData objects, relationships that can be
represented by existing feature nodes are added as directed peptide-feature to
protein-feature edges in the MuLink-compatible
``varp['feature_mapping']`` adjacency matrix.
"""

from __future__ import annotations

import json
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from loguru import logger
from prozor.annotate import annotate_peptides_streaming
from scipy.sparse import csr_matrix

from anndata_proteomics.annotation._sanitize import sanitize_columns
from anndata_proteomics.annotation.var_fasta import (
    leading_accession,
    protein_group_accessions,
    resolve_match_on,
)
from anndata_proteomics.fasta.anndata_io import write_fasta_config
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
from anndata_proteomics.fasta.parser import FastaSource, iter_fasta

_SCHEMA_VERSION = "0.2"
_DEFAULT_SEQUENCE_FIELD = "ProForma_peptide"
_VARM_KEY = "fasta_validation"
_FEATURE_MAPPING_KEY = "feature_mapping"
_OWNED_FEATURE_MAPPING_KEY = "_apb_fasta_feature_mapping_contribution"
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
_PEPTIDE_LEVELS = frozenset({"ion", "fragment", "peptidoform", "peptide"})
_LEADING_PROTEIN_FIELDS = (
    "Protein_Group",
    "Leading_Razor_Protein",
    "PG_ProteinAccessions",
    "PG_ProteinGroups",
    "Protein_ID",
    "Accession",
    "protein_group",
    "Proteins",
    "Protein",
)
_MAX_REPORTED = 5


@dataclass(slots=True)
class FastaValidationResult:
    """Validation result for one peptide-derived feature modality."""

    summary: pd.DataFrame
    matches: pd.DataFrame
    n_features: int
    n_unique_sequences: int
    n_invalid_sequences: int
    n_matched_features: int
    n_unmatched_features: int
    requested_backend: str
    backend: str
    sequence_field: str
    leading_protein_field: str | None
    il_equivalent: bool
    fasta_config: ResolvedFastaConfig
    unmatched_sequences: list[str]

    @property
    def fraction_unmatched(self) -> float:
        """Fraction of feature rows whose peptide was not found in the FASTA."""
        return self.n_unmatched_features / self.n_features if self.n_features else 0.0

    def sample_unmatched(self, k: int = _MAX_REPORTED) -> list[str]:
        """Return up to *k* distinct normalized unmatched sequences."""
        return self.unmatched_sequences[:k]


@dataclass(frozen=True, slots=True)
class MuLinkStorageStats:
    """Summary of the MuLink-compatible feature mapping update."""

    n_fasta_edges: int = 0
    n_unrepresented_fasta_proteins: int = 0
    protein_match_on: str | None = None


@dataclass(slots=True)
class _TargetInput:
    name: str
    target: Any
    normalized_sequences: pd.Series
    leading_proteins: pd.Series
    leading_protein_field: str | None


def validate_peptides_against_fasta(  # noqa: PLR0913 - stable public API
    obj: Any,
    fasta_sources: FastaSource | Iterable[FastaSource],
    *,
    sequence_field: str = _DEFAULT_SEQUENCE_FIELD,
    backend: str = "auto",
    fasta_config: FastaConfig | ResolvedFastaConfig | None = None,
    decoy_pattern: str | None = None,
    contaminant_pattern: str | None = None,
    leading_protein_field: str | None = None,
    protein_match_on: str | None = None,
    il_equivalent: bool = False,
    is_uniprot: bool = True,
    modality: str | None = None,
    store: bool = True,
) -> FastaValidationResult:
    """Validate one peptide-derived AnnData or selected MuData modality.

    The unmodified sequence comes from ``sequence_field``.  ``modality`` is only
    needed when a MuData has several peptide-derived modalities; the unified
    :func:`validate_peptide_modalities_against_fasta` operation validates all of
    them in one FASTA scan.
    """
    name, target = _resolve_feature_target(obj, modality)
    results = _validate_targets(
        obj,
        {name: target},
        fasta_sources,
        sequence_field=sequence_field,
        backend=backend,
        fasta_config=fasta_config,
        decoy_pattern=decoy_pattern,
        contaminant_pattern=contaminant_pattern,
        leading_protein_field=leading_protein_field,
        protein_match_on=protein_match_on,
        il_equivalent=il_equivalent,
        is_uniprot=is_uniprot,
        store=store,
    )
    return results[name]


def validate_peptide_modalities_against_fasta(  # noqa: PLR0913 - stable public API
    obj: Any,
    fasta_sources: FastaSource | Iterable[FastaSource],
    *,
    sequence_field: str = _DEFAULT_SEQUENCE_FIELD,
    backend: str = "auto",
    fasta_config: FastaConfig | ResolvedFastaConfig | None = None,
    decoy_pattern: str | None = None,
    contaminant_pattern: str | None = None,
    leading_protein_field: str | None = None,
    protein_match_on: str | None = None,
    il_equivalent: bool = False,
    is_uniprot: bool = True,
    store: bool = True,
) -> dict[str, FastaValidationResult]:
    """Validate every peptide-derived modality using one automaton and FASTA scan."""
    targets = _resolve_all_feature_targets(obj)
    return _validate_targets(
        obj,
        targets,
        fasta_sources,
        sequence_field=sequence_field,
        backend=backend,
        fasta_config=fasta_config,
        decoy_pattern=decoy_pattern,
        contaminant_pattern=contaminant_pattern,
        leading_protein_field=leading_protein_field,
        protein_match_on=protein_match_on,
        il_equivalent=il_equivalent,
        is_uniprot=is_uniprot,
        store=store,
    )


def _validate_targets(  # noqa: PLR0913 - explicit validation contract
    owner: Any,
    targets: dict[str, Any],
    fasta_sources: FastaSource | Iterable[FastaSource],
    *,
    sequence_field: str,
    backend: str,
    fasta_config: FastaConfig | ResolvedFastaConfig | None,
    decoy_pattern: str | None,
    contaminant_pattern: str | None,
    leading_protein_field: str | None,
    protein_match_on: str | None,
    il_equivalent: bool,
    is_uniprot: bool,
    store: bool,
) -> dict[str, FastaValidationResult]:
    if fasta_config is not None and (decoy_pattern is not None or contaminant_pattern is not None):
        raise ValueError("pass either fasta_config or decoy_pattern/contaminant_pattern, not both")
    if not targets:
        raise ValueError("object has no peptide-derived modality to validate")

    sources = materialize_sources(fasta_sources)
    source_descriptions = describe_sources(sources)
    prepared = {
        name: _prepare_target(
            name,
            target,
            sequence_field=sequence_field,
            leading_protein_field=leading_protein_field,
            il_equivalent=il_equivalent,
            is_uniprot=is_uniprot,
        )
        for name, target in targets.items()
    }
    patterns = sorted(
        {
            sequence
            for item in prepared.values()
            for sequence in item.normalized_sequences
            if sequence is not None
        }
    )
    input_config: FastaConfig | None = None
    resolved_config: ResolvedFastaConfig | None = None
    if isinstance(fasta_config, ResolvedFastaConfig):
        resolved_config = fasta_config
    elif isinstance(fasta_config, FastaConfig):
        input_config = fasta_config
    else:
        input_config = FastaConfig(
            decoy_patterns=_single_pattern(decoy_pattern),
            contaminant_patterns=_single_pattern(contaminant_pattern),
        )
    matches, fasta_proteins, resolved_config, resolved_backend = _scan_fasta(
        patterns,
        sources,
        backend=backend,
        input_config=input_config,
        resolved_config=resolved_config,
        il_equivalent=il_equivalent,
        is_uniprot=is_uniprot,
    )
    per_pattern = _per_pattern_stats(matches)
    results: dict[str, FastaValidationResult] = {}
    for name, item in prepared.items():
        summary = _build_summary(item, per_pattern, fasta_proteins)
        unique_sequences = {
            sequence for sequence in item.normalized_sequences if sequence is not None
        }
        unmatched_sequences = sorted(
            sequence for sequence in unique_sequences if sequence not in per_pattern
        )
        n_matched = int(summary["peptide_in_fasta"].sum())
        results[name] = FastaValidationResult(
            summary=summary,
            matches=matches[matches["sequence"].isin(unique_sequences)].copy(),
            n_features=len(summary),
            n_unique_sequences=len(unique_sequences),
            n_invalid_sequences=int(item.normalized_sequences.isna().sum()),
            n_matched_features=n_matched,
            n_unmatched_features=len(summary) - n_matched,
            requested_backend=backend,
            backend=resolved_backend,
            sequence_field=sequence_field,
            leading_protein_field=item.leading_protein_field,
            il_equivalent=il_equivalent,
            fasta_config=resolved_config,
            unmatched_sequences=unmatched_sequences,
        )

    mulink_stats = MuLinkStorageStats()
    if store:
        write_fasta_config(owner, resolved_config)
        if hasattr(owner, "mod") and "protein" in owner.mod:
            mulink_stats = _store_mulink_feature_mapping(
                owner,
                prepared,
                matches,
                protein_match_on=protein_match_on,
                is_uniprot=is_uniprot,
            )
        for name, item in prepared.items():
            _store(
                item.target,
                summary=results[name].summary,
                fasta_sources=source_descriptions,
                result=results[name],
                mulink_stats=mulink_stats,
            )

    total_features = sum(result.n_features for result in results.values())
    total_matched = sum(result.n_matched_features for result in results.values())
    logger.info(
        "FASTA validation: {}/{} features matched across {} modality/modalities; "
        "{} unique peptide patterns, {} match sites",
        total_matched,
        total_features,
        len(results),
        len(patterns),
        len(matches),
    )
    return results


def _prepare_target(
    name: str,
    target: Any,
    *,
    sequence_field: str,
    leading_protein_field: str | None,
    il_equivalent: bool,
    is_uniprot: bool,
) -> _TargetInput:
    normalized = _feature_sequences(target, sequence_field, il_equivalent=il_equivalent)
    resolved_leading_field = _resolve_leading_protein_field(target, leading_protein_field)
    leading = _feature_leading_proteins(
        target,
        resolved_leading_field,
        is_uniprot=is_uniprot,
    )
    return _TargetInput(
        name=name,
        target=target,
        normalized_sequences=normalized,
        leading_proteins=leading,
        leading_protein_field=resolved_leading_field,
    )


def _resolve_feature_target(obj: Any, modality: str | None) -> tuple[str, Any]:
    if hasattr(obj, "mod"):
        if modality is not None:
            if modality not in obj.mod:
                raise ValueError(
                    f"modality {modality!r} not in MuData; modalities: {list(obj.mod)}"
                )
            target = obj.mod[modality]
            if _level(target) not in _PEPTIDE_LEVELS:
                raise ValueError(f"modality {modality!r} is not peptide-derived")
            return modality, target
        candidates = _resolve_all_feature_targets(obj)
        if len(candidates) == 1:
            return next(iter(candidates.items()))
        if not candidates:
            raise ValueError("MuData has no peptide-derived modality")
        raise ValueError(
            f"MuData has multiple peptide-derived modalities {list(candidates)}; pass modality="
        )
    level = _level(obj)
    if level not in _PEPTIDE_LEVELS:
        raise ValueError(
            "FASTA validation applies to peptide-derived layers "
            f"(ion/fragment/peptidoform/peptide), got {level!r}"
        )
    return level or "features", obj


def _resolve_all_feature_targets(obj: Any) -> dict[str, Any]:
    if hasattr(obj, "mod"):
        return {
            name: target for name, target in obj.mod.items() if _level(target) in _PEPTIDE_LEVELS
        }
    level = _level(obj)
    return {level or "features": obj} if level in _PEPTIDE_LEVELS else {}


def _level(obj: Any) -> str | None:
    return (obj.uns.get("anndata_proteomics") or {}).get("quantification_level")


def _feature_sequences(target: Any, sequence_field: str, *, il_equivalent: bool) -> pd.Series:
    if sequence_field not in target.var.columns:
        raise ValueError(
            f"sequence_field {sequence_field!r} not in var columns: {list(target.var.columns)}"
        )
    values = [
        _normalize_sequence(value, il_equivalent=il_equivalent)
        for value in target.var[sequence_field]
    ]
    return pd.Series(values, index=pd.Index(target.var_names).astype(str), dtype="object")


def _normalize_sequence(value: object, *, il_equivalent: bool) -> str | None:
    if not isinstance(value, str):
        return None
    sequence = value.strip().upper()
    if not sequence:
        return None
    return sequence.replace("I", "L") if il_equivalent else sequence


def _resolve_leading_protein_field(target: Any, requested: str | None) -> str | None:
    if requested is not None:
        if requested not in target.var.columns:
            raise ValueError(
                f"leading_protein_field {requested!r} not in var columns: "
                f"{list(target.var.columns)}"
            )
        return requested
    return next(
        (field for field in _LEADING_PROTEIN_FIELDS if field in target.var.columns),
        None,
    )


def _feature_leading_proteins(
    target: Any,
    field: str | None,
    *,
    is_uniprot: bool,
) -> pd.Series:
    index = pd.Index(target.var_names).astype(str)
    if field is None:
        return pd.Series([None] * len(index), index=index, dtype="object")
    proteins: list[str | None] = []
    for value in target.var[field]:
        if pd.isna(value) or not str(value).strip():
            proteins.append(None)
        else:
            proteins.append(leading_accession(str(value), is_uniprot=is_uniprot))
    return pd.Series(proteins, index=index, dtype="object")


def _scan_fasta(  # noqa: C901 - single-pass FASTA state machine
    patterns: list[str],
    fasta_sources: list[FastaSource],
    *,
    backend: str,
    input_config: FastaConfig | None,
    resolved_config: ResolvedFastaConfig | None,
    il_equivalent: bool,
    is_uniprot: bool,
) -> tuple[pd.DataFrame, set[str], ResolvedFastaConfig, str]:
    accumulator = FastaConfigAccumulator(input_config) if input_config is not None else None
    fasta_proteins: set[str] = set()

    def protein_records() -> Iterable[tuple[str, str]]:
        for source in fasta_sources:
            for record in iter_fasta(source):
                fasta_id, _ = parse_header_id(record.header)
                if accumulator is not None:
                    accumulator.observe(fasta_id)
                proteinname = uniprot_proteinname(fasta_id) if is_uniprot else fasta_id
                fasta_proteins.add(proteinname)
                sequence = record.sequence.upper()
                if il_equivalent:
                    sequence = sequence.replace("I", "L")
                yield fasta_id, sequence

    records = protein_records()
    if patterns:
        annotations = annotate_peptides_streaming(patterns, records, backend=backend)
    else:
        for _ in records:
            pass
        annotations = annotate_peptides_streaming([], (), backend=backend)

    if resolved_config is None:
        if accumulator is None:  # pragma: no cover - guarded by _validate_targets
            raise AssertionError("an input or resolved FASTA configuration is required")
        effective_config = accumulator.resolve()
    else:
        effective_config = resolved_config
    rows = [
        (
            annotation.peptide,
            annotation.protein_id,
            uniprot_proteinname(annotation.protein_id) if is_uniprot else annotation.protein_id,
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
    return frame, fasta_proteins, effective_config, annotations.resolved_backend


def _per_pattern_stats(
    matches: pd.DataFrame,
) -> dict[str, tuple[int, tuple[str, ...]]]:
    """Map sequence to total match sites and distinct matching accessions."""
    if matches.empty:
        return {}
    stats: dict[str, tuple[int, tuple[str, ...]]] = {}
    for sequence, group in matches.groupby("sequence", sort=False):
        proteins = tuple(sorted(set(group["proteinname"].astype(str))))
        stats[str(sequence)] = (len(group), proteins)
    return stats


def _build_summary(
    item: _TargetInput,
    per_pattern: dict[str, tuple[int, tuple[str, ...]]],
    fasta_proteins: set[str],
) -> pd.DataFrame:
    peptide_in_fasta: list[bool] = []
    site_counts: list[int] = []
    protein_counts: list[int] = []
    protein_ids: list[str] = []
    leading_in_fasta: list[object] = []
    peptide_in_leading: list[object] = []

    for sequence, leading in zip(
        item.normalized_sequences,
        item.leading_proteins,
        strict=True,
    ):
        site_count, proteins = (
            per_pattern.get(sequence, (0, ())) if sequence is not None else (0, ())
        )
        peptide_in_fasta.append(site_count > 0)
        site_counts.append(site_count)
        protein_counts.append(len(proteins))
        protein_ids.append(";".join(proteins))
        if leading is None:
            leading_in_fasta.append(pd.NA)
            peptide_in_leading.append(pd.NA)
        else:
            leading_in_fasta.append(leading in fasta_proteins)
            peptide_in_leading.append(pd.NA if sequence is None else leading in proteins)

    return pd.DataFrame(
        {
            "peptide_in_fasta": peptide_in_fasta,
            "fasta_match_site_count": site_counts,
            "fasta_matching_protein_count": protein_counts,
            "fasta_matching_protein_ids": protein_ids,
            "leading_protein_in_fasta": pd.Series(
                leading_in_fasta,
                dtype="boolean",
            ).array,
            "peptide_in_leading_protein": pd.Series(
                peptide_in_leading,
                dtype="boolean",
            ).array,
        },
        index=pd.Index(item.target.var_names).astype(str),
    )


def _store_mulink_feature_mapping(  # noqa: C901, PLR0912 - modality topology validation
    mdata: Any,
    targets: dict[str, _TargetInput],
    matches: pd.DataFrame,
    *,
    protein_match_on: str | None,
    is_uniprot: bool,
) -> MuLinkStorageStats:
    if not mdata.var_names.is_unique:
        raise ValueError("MuLink feature_mapping requires globally unique MuData var_names")
    protein = mdata.mod["protein"]
    resolved_match_on = resolve_match_on(protein, protein_match_on)
    if resolved_match_on == "index":
        raw_groups = pd.Series(protein.var_names, index=protein.var_names)
    else:
        raw_groups = protein.var[resolved_match_on]

    accession_to_nodes: dict[str, list[str]] = {}
    for node, raw_group in zip(protein.var_names, raw_groups, strict=True):
        for accession in protein_group_accessions(
            str(raw_group),
            is_uniprot=is_uniprot,
        ):
            accession_to_nodes.setdefault(accession, []).append(str(node))

    matched_by_sequence = {
        str(sequence): set(group["proteinname"].astype(str))
        for sequence, group in matches.groupby("sequence", sort=False)
    }
    global_positions = {str(name): position for position, name in enumerate(mdata.var_names)}
    missing_feature_names = [
        str(feature_name)
        for item in targets.values()
        for feature_name in item.target.var_names
        if str(feature_name) not in global_positions
    ]
    if missing_feature_names:
        raise ValueError(
            "peptide-derived feature names are absent from the MuData global var axis: "
            f"{missing_feature_names[:_MAX_REPORTED]}"
        )
    target_row_mask = np.zeros(mdata.n_vars, dtype=bool)
    for item in targets.values():
        target_row_mask[
            [global_positions[str(feature_name)] for feature_name in item.target.var_names]
        ] = True
    rows: list[int] = []
    columns: list[int] = []
    represented_accessions: set[str] = set()
    all_matched_accessions = set(matches["proteinname"].astype(str))
    for item in targets.values():
        for feature_name, sequence in item.normalized_sequences.items():
            if sequence is None:
                continue
            for accession in matched_by_sequence.get(sequence, set()):
                for protein_node in accession_to_nodes.get(accession, []):
                    rows.append(global_positions[str(feature_name)])
                    columns.append(global_positions[protein_node])
                    represented_accessions.add(accession)

    new_mapping = csr_matrix(
        (
            np.ones(len(rows), dtype=np.int8),
            (
                np.asarray(rows, dtype=np.int64),
                np.asarray(columns, dtype=np.int64),
            ),
        ),
        shape=(mdata.n_vars, mdata.n_vars),
        # MuData's update path temporarily writes -1 while reindexing varp;
        # use a signed type so h5mu round-trips remain valid.
    )
    new_mapping.sum_duplicates()
    if new_mapping.nnz:
        new_mapping.data[:] = 1

    existing = csr_matrix(
        mdata.varp.get(
            _FEATURE_MAPPING_KEY,
            _empty_int8_csr(new_mapping.shape),
        )
    )
    if existing.shape != new_mapping.shape:
        raise ValueError(
            f"existing varp[{_FEATURE_MAPPING_KEY!r}] has shape {existing.shape}, "
            f"expected {new_mapping.shape}"
        )

    # ``feature_mapping`` is shared MuLink state.  Keep APB's additive
    # contribution separately so a later validation can replace only the
    # edges generated by the previous FASTA, without touching relationships
    # owned by another producer (including an existing edge at the same
    # coordinate).
    old_owned = csr_matrix(
        mdata.varp.get(
            _OWNED_FEATURE_MAPPING_KEY,
            _empty_int8_csr(new_mapping.shape),
        )
    )
    if old_owned.shape != new_mapping.shape:
        raise ValueError(
            f"existing varp[{_OWNED_FEATURE_MAPPING_KEY!r}] has shape "
            f"{old_owned.shape}, expected {new_mapping.shape}"
        )
    merged, owned = _replace_owned_feature_mapping(
        existing,
        old_owned,
        new_mapping,
        target_row_mask,
    )
    mdata.varp[_FEATURE_MAPPING_KEY] = merged
    mdata.varp[_OWNED_FEATURE_MAPPING_KEY] = owned

    return MuLinkStorageStats(
        n_fasta_edges=new_mapping.nnz,
        n_unrepresented_fasta_proteins=len(all_matched_accessions - represented_accessions),
        protein_match_on=resolved_match_on,
    )


def _empty_int8_csr(shape: tuple[int, int]) -> csr_matrix[np.int8]:
    """Create a typed empty signed matrix for MuData feature mappings."""
    empty_indices = np.empty(0, dtype=np.int64)
    return csr_matrix(
        (
            np.empty(0, dtype=np.int8),
            (empty_indices, empty_indices),
        ),
        shape=shape,
    )


def _replace_owned_feature_mapping(
    existing: csr_matrix,
    old_owned: csr_matrix,
    new_mapping: csr_matrix,
    target_row_mask: np.ndarray,
) -> tuple[csr_matrix, csr_matrix]:
    """Replace APB's contribution on selected rows of shared MuLink state."""
    old_owned_coo = old_owned.tocoo()
    targeted = target_row_mask[old_owned_coo.row]
    retained_owned = csr_matrix(
        (
            old_owned_coo.data[~targeted],
            (old_owned_coo.row[~targeted], old_owned_coo.col[~targeted]),
        ),
        shape=old_owned.shape,
        dtype=old_owned.dtype,
    )

    # An owned entry is removed only while its value is still the one APB
    # wrote.  If another producer changed that coordinate, preserve the new
    # value and relinquish ownership.  Filtering COO entries avoids casting or
    # subtracting the entire matrix, so even uint64 weights remain exact.
    target_owned_rows = old_owned_coo.row[targeted].astype(np.int64, copy=False)
    target_owned_cols = old_owned_coo.col[targeted].astype(np.int64, copy=False)
    target_owned_data = old_owned_coo.data[targeted]
    if target_owned_data.size:
        current = np.asarray(existing[target_owned_rows, target_owned_cols]).ravel()
        unchanged = current == target_owned_data
        remove_keys = (
            target_owned_rows[unchanged] * existing.shape[1] + target_owned_cols[unchanged]
        )
        existing_coo = existing.tocoo()
        existing_rows = existing_coo.row.astype(np.int64, copy=False)
        existing_cols = existing_coo.col.astype(np.int64, copy=False)
        existing_keys = existing_rows * existing.shape[1] + existing_cols
        keep = ~np.isin(existing_keys, remove_keys, assume_unique=True)
        base = csr_matrix(
            (
                existing_coo.data[keep],
                (existing_coo.row[keep], existing_coo.col[keep]),
            ),
            shape=existing.shape,
            dtype=existing.dtype,
        )
    else:
        base = existing.copy()

    occupied = base.astype(bool).astype("int8")
    new_owned = new_mapping - new_mapping.multiply(occupied)
    new_owned.eliminate_zeros()
    merged = base + new_owned.astype(base.dtype)
    merged.eliminate_zeros()
    return csr_matrix(merged), csr_matrix(retained_owned + new_owned)


def _store(
    target: Any,
    *,
    summary: pd.DataFrame,
    fasta_sources: list[str],
    result: FastaValidationResult,
    mulink_stats: MuLinkStorageStats,
) -> None:
    """Attach the per-feature validation frame and an auditable provenance entry."""
    stored = summary.copy()
    stored.columns = sanitize_columns(list(stored.columns))
    target.varm[_VARM_KEY] = stored

    entry = {
        "schema_version": _SCHEMA_VERSION,
        "source": "fasta_validation",
        "destination": f"varm[{_VARM_KEY!r}]",
        "feature_mapping": (
            f"varp[{_FEATURE_MAPPING_KEY!r}]" if mulink_stats.protein_match_on is not None else None
        ),
        "feature_mapping_ownership": (
            f"varp[{_OWNED_FEATURE_MAPPING_KEY!r}]"
            if mulink_stats.protein_match_on is not None
            else None
        ),
        "fasta_sources": fasta_sources,
        "fasta_config": result.fasta_config.model_dump(mode="json"),
        "requested_backend": result.requested_backend,
        "backend": result.backend,
        "sequence_field": result.sequence_field,
        "leading_protein_field": result.leading_protein_field,
        "il_equivalent": result.il_equivalent,
        "n_features": result.n_features,
        "n_unique_sequences": result.n_unique_sequences,
        "n_invalid_sequences": result.n_invalid_sequences,
        "n_matched_features": result.n_matched_features,
        "n_unmatched_features": result.n_unmatched_features,
        "n_feature_mapping_edges": mulink_stats.n_fasta_edges,
        "n_unrepresented_fasta_proteins": (mulink_stats.n_unrepresented_fasta_proteins),
        "protein_match_on": mulink_stats.protein_match_on,
    }
    namespace = dict(target.uns.get("anndata_proteomics", {}))
    existing = json.loads(namespace.get("var_annotations_json", "[]"))
    existing.append(entry)
    namespace["var_annotations_json"] = json.dumps(existing)
    target.uns["anndata_proteomics"] = namespace


def _single_pattern(pattern: str | None) -> tuple[str, ...] | None:
    if pattern is None:
        return None
    return (pattern,) if pattern else ()
