"""Pure MuLink feature-node assembly and shared-adjacency merging.

Node assembly and the owned-contribution merge are calculations over ordinary pandas,
NumPy, and SciPy values. They live here rather than in a workflow or a storage adapter so
that orchestration only orders them and any backend can supply their inputs.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy.sparse import csr_matrix

from anndata_proteomics.annotation.validate_fasta import (
    PeptideFeatureNode,
    PeptideFeatureNodes,
    PeptideProteinMatches,
    ProteinFeatureNodes,
    peptide_protein_matches,
)
from anndata_proteomics.annotation.var_fasta import protein_group_accessions

_MAX_REPORTED = 5


@dataclass(frozen=True, slots=True)
class FeaturePositions:
    """Every feature's row/column position on one global feature axis."""

    by_name: dict[str, int]
    total: int

    def require(self, name: str) -> int:
        """Return one feature's position or raise precisely."""
        if name not in self.by_name:
            raise ValueError(f"feature name {name!r} is absent from the global feature axis")
        return self.by_name[name]


@dataclass(frozen=True, slots=True)
class OwnedFeatureMappingMerge:
    """Shared adjacency and APB's own contribution after one replacement."""

    merged: csr_matrix
    owned: csr_matrix


def feature_positions(global_feature_names: pd.Index) -> FeaturePositions:
    """Index one global feature axis by name, requiring unique names."""
    if not global_feature_names.is_unique:
        raise ValueError("MuLink feature_mapping requires globally unique feature names")
    by_name = {
        str(feature_name): position for position, feature_name in enumerate(global_feature_names)
    }
    return FeaturePositions(by_name=by_name, total=len(global_feature_names))


def require_feature_positions(names: Sequence[str], positions: FeaturePositions) -> list[int]:
    """Resolve every name to its axis position, reporting all absent names at once."""
    missing = [name for name in names if name not in positions.by_name]
    if missing:
        raise ValueError(
            "peptide-derived feature names are absent from the global feature axis: "
            f"{missing[:_MAX_REPORTED]}"
        )
    return [positions.by_name[name] for name in names]


def peptide_feature_nodes(
    normalized_sequences: Iterable[pd.Series],
    positions: FeaturePositions,
) -> PeptideFeatureNodes:
    """Build peptide domain nodes from normalized per-level sequence series."""
    nodes: list[PeptideFeatureNode] = []
    for level_sequences in normalized_sequences:
        for feature_name, sequence in level_sequences.items():
            if isinstance(sequence, str):
                nodes.append(
                    PeptideFeatureNode(
                        position=positions.require(str(feature_name)),
                        sequence=sequence,
                    )
                )
    return PeptideFeatureNodes(nodes=tuple(nodes), total_nodes=positions.total)


def protein_feature_nodes(
    protein_groups: pd.Series,
    positions: FeaturePositions,
    *,
    is_uniprot: bool,
) -> ProteinFeatureNodes:
    """Build protein domain nodes from explicit protein-group values."""
    positions_by_accession: dict[str, list[int]] = {}
    for feature_name, group in protein_groups.items():
        position = positions.require(str(feature_name))
        for accession in protein_group_accessions(str(group), is_uniprot=is_uniprot):
            positions_by_accession.setdefault(accession, []).append(position)
    return ProteinFeatureNodes(
        positions_by_accession={
            accession: tuple(node_positions)
            for accession, node_positions in positions_by_accession.items()
        },
        total_nodes=positions.total,
    )


def combined_peptide_protein_matches(
    match_frames: Iterable[pd.DataFrame],
) -> PeptideProteinMatches:
    """Combine and deduplicate match sites across validated levels."""
    combined = pd.concat(list(match_frames), ignore_index=True).drop_duplicates()
    return peptide_protein_matches(combined)


def target_row_mask(row_positions: Iterable[int], positions: FeaturePositions) -> NDArray[np.bool_]:
    """Mark the axis rows this enrichment owns."""
    mask = np.zeros(positions.total, dtype=np.bool_)
    mask[list(row_positions)] = True
    return mask


def empty_feature_mapping(shape: tuple[int, int]) -> csr_matrix[np.int8]:
    """Return the all-absent adjacency used when no mapping is stored yet."""
    empty_indices = np.empty(0, dtype=np.int64)
    return csr_matrix(
        (
            np.empty(0, dtype=np.int8),
            (empty_indices, empty_indices),
        ),
        shape=shape,
    )


def merge_owned_feature_mapping(
    existing: csr_matrix,
    old_owned: csr_matrix,
    new_mapping: csr_matrix,
    target_rows: NDArray[np.bool_],
) -> OwnedFeatureMappingMerge:
    """Replace APB's contribution on selected rows of shared MuLink state.

    Rows outside ``target_rows`` keep both their shared edges and APB's recorded
    ownership. On targeted rows APB's previous contribution is withdrawn from the shared
    matrix only where the stored value still matches what APB wrote, so edges another
    producer has since changed are never removed.
    """
    old_owned_coo = old_owned.tocoo()
    old_owned_coo.sum_duplicates()
    targeted = target_rows[old_owned_coo.row]
    retained_owned = csr_matrix(
        (
            old_owned_coo.data[~targeted],
            (old_owned_coo.row[~targeted], old_owned_coo.col[~targeted]),
        ),
        shape=old_owned.shape,
        dtype=old_owned.dtype,
    )

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
        existing_coo.sum_duplicates()
        existing_rows = existing_coo.row.astype(np.int64, copy=False)
        existing_cols = existing_coo.col.astype(np.int64, copy=False)
        existing_keys = existing_rows * existing.shape[1] + existing_cols
        keep = ~np.isin(existing_keys, remove_keys)
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
    return OwnedFeatureMappingMerge(
        merged=csr_matrix(merged),
        owned=csr_matrix(retained_owned + new_owned),
    )
