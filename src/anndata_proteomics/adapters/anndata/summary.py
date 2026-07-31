"""Extract lightweight description inputs from live AnnData and MuData objects."""

from __future__ import annotations

from anndata import AnnData
from mudata import MuData

from anndata_proteomics.adapters.anndata.namespace import read_namespace
from anndata_proteomics.adapters.anndata.summary_metadata import parse_description_metadata
from anndata_proteomics.description import (
    AnnDataDescriptionSource,
    MuDataDescriptionSource,
)


def read_anndata_description_source(target: AnnData) -> AnnDataDescriptionSource:
    """Extract only shape, layer names, and APB-owned metadata from AnnData."""
    return AnnDataDescriptionSource(
        n_runs=int(target.n_obs),
        n_features=int(target.n_vars),
        layers=tuple(str(name) for name in target.layers.keys() if name is not None),
        metadata=parse_description_metadata(read_namespace(target)),
    )


def read_mudata_description_source(target: MuData) -> MuDataDescriptionSource:
    """Extract only root and modality metadata from MuData."""
    modalities = {
        str(name): read_anndata_description_source(_require_anndata(str(name), modality))
        for name, modality in target.mod.items()
    }
    return MuDataDescriptionSource(
        n_runs=int(target.n_obs),
        n_features=int(target.n_vars),
        modalities=modalities,
        metadata=parse_description_metadata(read_namespace(target)),
    )


def _require_anndata(name: str, target: AnnData | MuData) -> AnnData:
    if not isinstance(target, AnnData):
        raise TypeError(f"MuData modality {name!r} is not an AnnData")
    return target
