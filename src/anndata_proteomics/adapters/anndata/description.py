"""AnnData composition of the backend-independent description workflow.

Both composition roots — the ``apb`` CLI and APB Studio — call these functions rather than
wiring the HDF5 readers themselves.
"""

from __future__ import annotations

from pathlib import Path

from anndata_proteomics.adapters.anndata import summary_hdf5
from anndata_proteomics.description import JsonObject
from anndata_proteomics.workflows import summary as summary_workflow

HDF5_DESCRIPTION_READERS = summary_workflow.StoredDescriptionReaders(
    read_anndata=summary_hdf5.read_anndata_description_source,
    read_mudata=summary_hdf5.read_mudata_description_source,
    read_mudata_modality=summary_hdf5.read_mudata_modality_description_source,
)


def describe_path(path: Path | str) -> JsonObject:
    """Describe one stored HDF5 container without reading quantitative values."""
    return summary_workflow.describe_path(path, HDF5_DESCRIPTION_READERS)


def describe_modality_path(path: Path | str, modality: str) -> JsonObject:
    """Describe one named modality in a stored HDF5 MuData container."""
    return summary_workflow.describe_modality_path(path, modality, HDF5_DESCRIPTION_READERS)
