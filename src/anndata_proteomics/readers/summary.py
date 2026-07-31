"""Lightweight presentation views over converted AnnData and MuData containers."""

from __future__ import annotations

import json
from collections.abc import Mapping
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import h5py
import numpy as np
from anndata.io import read_elem

from anndata_proteomics._matrix_types import named_layers
from anndata_proteomics.params.anndata_io import read_search_parameters
from anndata_proteomics.rules.anndata_io import read_stored_rule
from anndata_proteomics.rules.schema import ColumnGroup

_NAMESPACE = "anndata_proteomics"
_VIEW_SCHEMA_VERSION = "1"
_CONVERSION_FIELDS = (
    "schema_version",
    "quantification_level",
    "software_name",
    "rule_selection_method",
    "search_parameters_path",
    "search_parameters_version_status",
)


@dataclass(frozen=True)
class _AnnDataMetadata:
    """The AnnData fields needed by :func:`describe`, without quantitative values."""

    n_obs: int
    n_vars: int
    layers: dict[str, None]
    uns: dict[str, Any]


@dataclass(frozen=True)
class _MuDataMetadata:
    """The MuData fields needed by :func:`describe`, without quantitative values."""

    n_obs: int
    n_vars: int
    mod: dict[str, _AnnDataMetadata]
    uns: dict[str, Any]


def describe(obj: Any) -> dict[str, Any]:
    """Return a JSON-compatible view without scanning ``X`` or quantitative layers."""
    if _is_mudata(obj):
        result = {
            "schema_version": _VIEW_SCHEMA_VERSION,
            "container_type": "mudata",
            "quantification": {
                "n_runs": int(obj.n_obs),
                "n_features": int(obj.n_vars),
                "modalities": list(obj.mod),
            },
            "modalities": {name: describe(modality) for name, modality in obj.mod.items()},
        }
        _add_metadata_views(result, obj)
        return result

    metadata = obj.uns.get(_NAMESPACE) or {}
    params = read_search_parameters(obj)
    result = {
        "schema_version": _VIEW_SCHEMA_VERSION,
        "container_type": "anndata",
        "quantification": {
            "n_runs": int(obj.n_obs),
            "n_features": int(obj.n_vars),
            "level": metadata.get("quantification_level"),
            "software_name": metadata.get("software_name"),
            "software_version": params.software_version if params is not None else None,
            "layers": [str(name) for name in named_layers(obj)],
        },
    }
    column_mapping = _column_mapping(obj)
    if column_mapping is not None:
        result["column_mapping"] = column_mapping
    _add_metadata_views(result, obj)
    return result


def describe_path(path: Path | str, modality: str | None = None) -> dict[str, Any]:
    """Describe one HDF5 container without loading ``X`` or layer datasets."""
    result_path = Path(path).expanduser()
    if result_path.suffix == ".h5ad":
        if modality is not None:
            raise ValueError("--modality applies only to MuData (.h5mu) inputs")
        with h5py.File(result_path, "r") as store:
            return describe(_read_anndata_metadata(store))

    if result_path.suffix == ".h5mu":
        with h5py.File(result_path, "r") as store:
            obj = _read_mudata_metadata(store)
            if modality is None:
                return describe(obj)
            if modality not in obj.mod:
                raise ValueError(
                    f"modality {modality!r} not in MuData; modalities: {list(obj.mod)}"
                )
            return describe(obj.mod[modality])

    raise ValueError(f"unsupported converted result type: {result_path}")


def _read_anndata_metadata(store: h5py.Group | h5py.File) -> _AnnDataMetadata:
    """Read only shape, layer names, and APB-owned ``uns`` metadata."""
    shape = _matrix_shape(store.get("X"))
    if shape is None:
        shape = (_frame_length(store, "obs"), _frame_length(store, "var"))
    layers = store.get("layers")
    layer_names = list(layers) if isinstance(layers, h5py.Group) else []
    return _AnnDataMetadata(
        n_obs=shape[0],
        n_vars=shape[1],
        layers={str(name): None for name in layer_names},
        uns=_read_apb_namespace(store),
    )


def _read_mudata_metadata(store: h5py.File) -> _MuDataMetadata:
    """Read the root and modality metadata needed for a MuData view."""
    modalities = store.get("mod")
    if not isinstance(modalities, h5py.Group):
        raise ValueError("invalid MuData container: missing 'mod' group")
    order = _modality_order(modalities)
    mod = {name: _read_anndata_metadata(_require_group(modalities, name)) for name in order}
    return _MuDataMetadata(
        n_obs=_frame_length(store, "obs"),
        n_vars=_frame_length(store, "var"),
        mod=mod,
        uns=_read_apb_namespace(store),
    )


def _matrix_shape(element: h5py.Group | h5py.Dataset | Any | None) -> tuple[int, int] | None:
    """Return an encoded dense/sparse matrix shape without reading its values."""
    if isinstance(element, h5py.Dataset):
        if len(element.shape) != 2:
            raise ValueError(f"invalid AnnData X shape: {element.shape}")
        return int(element.shape[0]), int(element.shape[1])
    if isinstance(element, h5py.Group):
        shape = element.attrs.get("shape")
        if shape is not None and len(shape) == 2:
            return int(shape[0]), int(shape[1])
    return None


def _frame_length(store: h5py.Group | h5py.File, axis: str) -> int:
    """Read an encoded dataframe's row count from index dataset shapes."""
    frame = store.get(axis)
    if not isinstance(frame, h5py.Group):
        raise ValueError(f"invalid container: missing {axis!r} dataframe")
    raw_index = frame.attrs.get("_index")
    if isinstance(raw_index, bytes):
        raw_index = raw_index.decode("utf-8")
    if not isinstance(raw_index, str):
        raise ValueError(f"invalid container: {axis!r} dataframe has no string index key")
    index = frame.get(raw_index)
    length = _encoded_vector_length(index)
    if length is None:
        raise ValueError(f"invalid container: cannot determine {axis!r} dataframe length")
    return length


def _encoded_vector_length(element: h5py.Group | h5py.Dataset | Any | None) -> int | None:
    """Return a one-dimensional encoded element's length from HDF5 metadata."""
    if isinstance(element, h5py.Dataset):
        return int(element.shape[0]) if element.shape else None
    if isinstance(element, h5py.Group):
        for key in ("values", "codes", "mask"):
            child = element.get(key)
            length = _encoded_vector_length(child)
            if length is not None:
                return length
    return None


def _read_apb_namespace(store: h5py.Group | h5py.File) -> dict[str, Any]:
    """Decode only APB's ``uns`` namespace, never another HDF5 slot."""
    uns = store.get("uns")
    if not isinstance(uns, h5py.Group):
        return {}
    namespace = uns.get(_NAMESPACE)
    if namespace is None:
        return {}
    if not isinstance(namespace, (h5py.Group, h5py.Dataset)):
        raise ValueError(f"invalid uns[{_NAMESPACE!r}] HDF5 element: {type(namespace).__name__}")
    decoded = read_elem(namespace)
    if not isinstance(decoded, Mapping):
        raise ValueError(
            f"invalid uns[{_NAMESPACE!r}] value: expected mapping, got {type(decoded).__name__}"
        )
    return {_NAMESPACE: dict(decoded)}


def _modality_order(modalities: h5py.Group) -> list[str]:
    """Return stored MuData modality order, with group order as a fallback."""
    raw = modalities.attrs.get("mod-order")
    if raw is None:
        return list(modalities)
    values = raw.tolist() if isinstance(raw, np.ndarray) else list(raw)
    return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values]


def _require_group(parent: h5py.Group, key: str) -> h5py.Group:
    value = parent.get(key)
    if not isinstance(value, h5py.Group):
        raise ValueError(f"invalid MuData container: modality {key!r} is not a group")
    return value


def _add_metadata_views(result: dict[str, Any], obj: Any) -> None:
    namespace = obj.uns.get(_NAMESPACE)
    if not isinstance(namespace, Mapping):
        return

    conversion = {
        field: _to_json_compatible(namespace[field])
        for field in _CONVERSION_FIELDS
        if field in namespace
    }
    if conversion:
        result["conversion"] = conversion

    params = read_search_parameters(obj)
    if params is not None:
        result["search_parameters"] = params.model_dump(mode="json")

    annotations = _annotation_provenance(namespace)
    if annotations:
        result["annotations"] = annotations

    qc = _decode_json_value(namespace.get("qc"))
    if qc is not None:
        result["qc"] = qc

    proteobench = namespace.get("proteobench")
    if isinstance(proteobench, Mapping) and "scores" in proteobench:
        result["proteobench"] = _to_json_compatible(proteobench)


def _annotation_provenance(namespace: Mapping[str, Any]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    obs = _decode_json_value(namespace.get("obs_annotations_json"))
    if obs:
        result["obs"] = obs
    var = _decode_json_value(namespace.get("var_annotations_json"))
    if var:
        result["var"] = var
    fasta_config = _decode_json_value(namespace.get("fasta_config"))
    if fasta_config is not None:
        result["fasta_config"] = fasta_config
    return result


def _decode_json_value(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    if isinstance(value, str):
        return _to_json_compatible(json.loads(value))
    return _to_json_compatible(value)


def _column_mapping(obj: Any) -> dict[str, Any] | None:
    """Describe where the effective rule placed vendor and computed columns."""
    rule = read_stored_rule(obj)
    if rule is None:
        return None

    layers_by_name = {layer.name: layer for layer in rule.layers}
    materialized_layers = set(named_layers(obj)) | {rule.axis.x_layer}
    layers = {
        name: {
            "source": layer.source,
            "source_kind": "column" if rule.input_shape == "long" else "pattern",
        }
        for name, layer in layers_by_name.items()
        if name in materialized_layers
    }
    x_layer = layers_by_name[rule.axis.x_layer]
    source_kind = "column" if rule.input_shape == "long" else "pattern"
    return {
        "X": {
            "layer": x_layer.name,
            "source": x_layer.source,
            "source_kind": source_kind,
        },
        "layers": layers,
        "obs": _column_group_mapping(rule.columns.obs),
        "var": _column_group_mapping(rule.columns.var),
    }


def _column_group_mapping(group: ColumnGroup) -> dict[str, str]:
    """Map output column names to their vendor source or compute operation."""
    mapping = dict(group.select)
    mapping.update({column.name: f"computed:{column.how}" for column in group.compute})
    return mapping


def _to_json_compatible(value: Any) -> Any:
    """Copy an HDF5-decoded value into ordinary JSON-compatible values."""
    if isinstance(value, Mapping):
        return {str(key): _to_json_compatible(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return [_to_json_compatible(item) for item in value.tolist()]
    if isinstance(value, (list, tuple)):
        return [_to_json_compatible(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return deepcopy(value)


def _is_mudata(obj: Any) -> bool:
    return hasattr(obj, "mod")
