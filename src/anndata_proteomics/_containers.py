"""Structural types for the AnnData-shaped objects APB reads metadata from.

Most APB operations take a real :class:`anndata.AnnData` or :class:`mudata.MuData` and
should say so. These protocols exist for the two readers that deliberately accept less:
``readers.summary`` describes an HDF5 container from shape and ``uns`` alone, without
loading ``X`` or any layer dataset, so it feeds its own lightweight metadata records
through the same accessors a real container would use.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, Protocol


class UnsHolder(Protocol):
    """Any object exposing an AnnData-style ``uns`` mapping."""

    @property
    def uns(self) -> Mapping[str, Any]:
        """Unstructured metadata, including APB's ``anndata_proteomics`` namespace."""
        ...
