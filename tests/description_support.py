"""Test composition for live-container description extraction and calculation."""

from __future__ import annotations

from anndata import AnnData
from mudata import MuData

from anndata_proteomics.adapters.anndata.summary import (
    read_anndata_description_source,
    read_mudata_description_source,
)
from anndata_proteomics.description import (
    calculate_anndata_description,
    calculate_mudata_description,
)
from anndata_proteomics.serialization import JsonObject


def describe_anndata(target: AnnData) -> JsonObject:
    """Compose live AnnData extraction with pure description calculation."""
    return calculate_anndata_description(read_anndata_description_source(target))


def describe_mudata(target: MuData) -> JsonObject:
    """Compose live MuData extraction with pure description calculation."""
    return calculate_mudata_description(read_mudata_description_source(target))
