"""Extract typed description inputs directly from AnnData-family HDF5 files."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
from anndata.io import read_elem

from anndata_proteomics.adapters.anndata.namespace import NAMESPACE, parse_namespace
from anndata_proteomics.adapters.anndata.summary_metadata import parse_description_metadata
from anndata_proteomics.description import (
    AnnDataDescriptionSource,
    MuDataDescriptionSource,
)
from anndata_proteomics.serialization import JsonObject


@dataclass(frozen=True, slots=True)
class _MissingShape:
    """An HDF5 element does not encode a two-dimensional matrix shape."""


_MISSING_SHAPE = _MissingShape()


def read_anndata_description_source(path: Path) -> AnnDataDescriptionSource:
    """Read shape, layer names, and APB metadata from one .h5ad file."""
    with h5py.File(path, "r") as store:
        return _read_anndata_group(store)


def read_mudata_description_source(path: Path) -> MuDataDescriptionSource:
    """Read root and modality description inputs from one .h5mu file."""
    with h5py.File(path, "r") as store:
        modalities = _require_modalities(store)
        order = _modality_order(modalities)
        sources = {name: _read_anndata_group(_require_group(modalities, name)) for name in order}
        return MuDataDescriptionSource(
            n_runs=_frame_length(store, "obs"),
            n_features=_frame_length(store, "var"),
            modalities=sources,
            metadata=parse_description_metadata(_read_apb_namespace(store)),
        )


def read_mudata_modality_description_source(
    path: Path,
    modality: str,
) -> AnnDataDescriptionSource:
    """Read description inputs for one named modality in a .h5mu file."""
    with h5py.File(path, "r") as store:
        modalities = _require_modalities(store)
        order = _modality_order(modalities)
        if modality not in modalities:
            raise ValueError(f"modality {modality!r} not in MuData; modalities: {order}")
        return _read_anndata_group(_require_group(modalities, modality))


def _read_anndata_group(store: h5py.Group | h5py.File) -> AnnDataDescriptionSource:
    """Read only shape, layer names, and APB-owned ``uns`` metadata."""
    n_runs, n_features = _matrix_or_axis_shape(store)
    layers = store.get("layers")
    layer_names = tuple(str(name) for name in layers) if isinstance(layers, h5py.Group) else ()
    return AnnDataDescriptionSource(
        n_runs=n_runs,
        n_features=n_features,
        layers=layer_names,
        metadata=parse_description_metadata(_read_apb_namespace(store)),
    )


def _matrix_or_axis_shape(store: h5py.Group | h5py.File) -> tuple[int, int]:
    element = store.get("X")
    if isinstance(element, h5py.Dataset):
        if len(element.shape) != 2:
            raise ValueError(f"invalid AnnData X shape: {element.shape}")
        return int(element.shape[0]), int(element.shape[1])
    if isinstance(element, h5py.Group):
        shape = _group_matrix_shape(element)
        if not isinstance(shape, _MissingShape):
            return shape
    return _frame_length(store, "obs"), _frame_length(store, "var")


def _group_matrix_shape(element: h5py.Group) -> tuple[int, int] | _MissingShape:
    """Decode a sparse matrix group's two-dimensional shape attribute."""
    value = element.attrs.get("shape")
    if isinstance(value, np.ndarray):
        if value.ndim != 1 or value.size != 2:
            return _MISSING_SHAPE
        first = value[0]
        second = value[1]
    elif isinstance(value, (list, tuple)) and len(value) == 2:
        first = value[0]
        second = value[1]
    else:
        return _MISSING_SHAPE
    if not isinstance(first, (int, np.integer)):
        raise ValueError(f"invalid matrix shape value: {first!r}")
    if not isinstance(second, (int, np.integer)):
        raise ValueError(f"invalid matrix shape value: {second!r}")
    return int(first), int(second)


def _frame_length(store: h5py.Group | h5py.File, axis: str) -> int:
    """Read an encoded dataframe's row count from its index dataset shape."""
    frame = store.get(axis)
    if not isinstance(frame, h5py.Group):
        raise ValueError(f"invalid container: missing {axis!r} dataframe")
    raw_index = frame.attrs.get("_index")
    if isinstance(raw_index, bytes):
        raw_index = raw_index.decode("utf-8")
    if not isinstance(raw_index, str):
        raise ValueError(f"invalid container: {axis!r} dataframe has no string index key")
    element = frame.get(raw_index)
    if not isinstance(element, (h5py.Group, h5py.Dataset)):
        raise ValueError(f"invalid container: cannot determine {axis!r} dataframe length")
    return _encoded_vector_length(element, axis)


def _encoded_vector_length(element: h5py.Group | h5py.Dataset, axis: str) -> int:
    """Return a one-dimensional encoded element's length from HDF5 metadata."""
    if isinstance(element, h5py.Dataset):
        if not element.shape:
            raise ValueError(f"invalid container: cannot determine {axis!r} dataframe length")
        return int(element.shape[0])
    for key in ("values", "codes", "mask"):
        child = element.get(key)
        if isinstance(child, (h5py.Group, h5py.Dataset)):
            return _encoded_vector_length(child, axis)
    raise ValueError(f"invalid container: cannot determine {axis!r} dataframe length")


def _read_apb_namespace(store: h5py.Group | h5py.File) -> JsonObject:
    """Decode only APB's ``uns`` namespace, never another HDF5 slot."""
    uns = store.get("uns")
    if not isinstance(uns, h5py.Group):
        return {}
    namespace = uns.get(NAMESPACE)
    if namespace is None:
        return {}
    if not isinstance(namespace, (h5py.Group, h5py.Dataset)):
        raise ValueError(f"invalid uns[{NAMESPACE!r}] HDF5 element: {type(namespace).__name__}")
    decoded = read_elem(namespace)
    return parse_namespace(decoded)


def _require_modalities(store: h5py.File) -> h5py.Group:
    modalities = store.get("mod")
    if not isinstance(modalities, h5py.Group):
        raise ValueError("invalid MuData container: missing 'mod' group")
    return modalities


def _modality_order(modalities: h5py.Group) -> list[str]:
    """Return stored MuData modality order, with group order as a fallback."""
    raw = modalities.attrs.get("mod-order")
    if raw is None:
        return [str(name) for name in modalities]
    if isinstance(raw, np.ndarray):
        if raw.ndim != 1:
            raise ValueError("invalid MuData container: 'mod-order' must be a sequence")
        values = raw.tolist()
    elif isinstance(raw, (list, tuple)):
        values = raw
    else:
        raise ValueError("invalid MuData container: 'mod-order' must be a sequence")
    names: list[str] = []
    for value in values:
        if not isinstance(value, (str, bytes)):
            raise ValueError(f"invalid MuData modality name: {value!r}")
        names.append(_modality_name(value))
    return names


def _modality_name(value: str | bytes) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return value


def _require_group(parent: h5py.Group, key: str) -> h5py.Group:
    value = parent.get(key)
    if not isinstance(value, h5py.Group):
        raise ValueError(f"invalid MuData container: modality {key!r} is not a group")
    return value
