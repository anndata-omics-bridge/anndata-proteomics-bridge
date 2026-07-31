"""Backend-independent orchestration for lightweight APB container descriptions.

The suffix decides which extraction a stored artifact needs; the extraction itself is
injected, so this module names no storage backend. The AnnData composition lives in
``adapters/anndata/description.py``.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from anndata_proteomics.description import (
    AnnDataDescriptionSource,
    JsonObject,
    MuDataDescriptionSource,
    calculate_anndata_description,
    calculate_mudata_description,
)

ANNDATA_SUFFIX = ".h5ad"
MUDATA_SUFFIX = ".h5mu"


@dataclass(frozen=True, slots=True)
class StoredDescriptionReaders:
    """One backend's file-oriented description extraction functions."""

    read_anndata: Callable[[Path], AnnDataDescriptionSource]
    read_mudata: Callable[[Path], MuDataDescriptionSource]
    read_mudata_modality: Callable[[Path, str], AnnDataDescriptionSource]


def describe_path(path: Path | str, readers: StoredDescriptionReaders) -> JsonObject:
    """Describe one complete stored container without reading quantitative values."""
    result_path = Path(path).expanduser()
    if result_path.suffix == ANNDATA_SUFFIX:
        return calculate_anndata_description(readers.read_anndata(result_path))
    if result_path.suffix == MUDATA_SUFFIX:
        return calculate_mudata_description(readers.read_mudata(result_path))
    raise ValueError(f"unsupported converted result type: {result_path}")


def describe_modality_path(
    path: Path | str,
    modality: str,
    readers: StoredDescriptionReaders,
) -> JsonObject:
    """Describe one named modality in a stored multi-level container."""
    result_path = Path(path).expanduser()
    if result_path.suffix == ANNDATA_SUFFIX:
        raise ValueError("--modality applies only to MuData (.h5mu) inputs")
    if result_path.suffix != MUDATA_SUFFIX:
        raise ValueError(f"unsupported converted result type: {result_path}")
    return calculate_anndata_description(readers.read_mudata_modality(result_path, modality))
