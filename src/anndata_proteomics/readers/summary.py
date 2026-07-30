"""Stage-owned descriptive summaries for converted AnnData and MuData objects."""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np

from anndata_proteomics._matrix_types import is_sparse_matrix, named_layers
from anndata_proteomics.params.anndata_io import read_search_parameters
from anndata_proteomics.rules.schema import ColumnGroup, ParseRule

_NAMESPACE = "anndata_proteomics"
_SUMMARY_KEY = "descriptive_summary"
_SCHEMA_VERSION = "5"
_FASTA_ANNOTATION_KEY = "fasta"
_FASTA_VALIDATION_KEY = "fasta_validation"
_FASTA_MATCHED_KEY = "peptide_in_fasta"
_FASTA_PROTEIN_COUNT = "fasta_matching_protein_count"


def store_quantification_summary(obj: Any) -> None:
    """Compute and store the conversion-owned summary component.

    AnnData quantification metrics are computed from its final layers. For MuData,
    modalities that do not already carry that component are summarized and the
    container-level shape and modality index are stored on the MuData object.
    """
    if _is_mudata(obj):
        for modality in obj.mod.values():
            payload = _read_payload(modality)
            if "quantification" not in payload:
                _store_quantification_summary_anndata(modality)
        payload = _read_payload(obj)
        payload.update(
            {
                "schema_version": _SCHEMA_VERSION,
                "container_type": "mudata",
                "quantification": {
                    "n_runs": int(obj.n_obs),
                    "n_features": int(obj.n_vars),
                    "modalities": list(obj.mod),
                },
            }
        )
        _write_payload(obj, payload)
        return
    _store_quantification_summary_anndata(obj)


def store_annotation_summary(
    obj: Any,
    *,
    fields: Sequence[str],
    n_annotated_runs: int,
) -> None:
    """Merge a compact sample-annotation component without recomputing conversion metrics."""
    targets = (obj, *obj.mod.values()) if _is_mudata(obj) else (obj,)
    for target in targets:
        payload = _read_payload(target)
        payload.setdefault("schema_version", _SCHEMA_VERSION)
        payload.setdefault(
            "container_type",
            "mudata" if _is_mudata(target) else "anndata",
        )
        payload["annotation"] = {
            "annotated_run_count": int(n_annotated_runs),
            "fields": list(fields),
            "group_counts": {
                field: int(target.obs[field].nunique(dropna=True)) for field in fields
            },
        }
        _write_payload(target, payload)


def store_fasta_summary(obj: Any) -> None:
    """Merge a small level-specific FASTA component into each applicable modality."""
    targets = obj.mod.values() if _is_mudata(obj) else (obj,)
    for target in targets:
        fasta_summary = _fasta_summary(target)
        if fasta_summary is None:
            continue
        payload = _read_payload(target)
        payload.setdefault("schema_version", _SCHEMA_VERSION)
        payload.setdefault("container_type", "anndata")
        payload["fasta"] = fasta_summary
        _write_payload(target, payload)


def _fasta_summary(target: Any) -> dict[str, int] | None:
    validation = target.varm.get(_FASTA_VALIDATION_KEY)
    if validation is not None and _FASTA_PROTEIN_COUNT in validation.columns:
        protein_counts = np.asarray(validation[_FASTA_PROTEIN_COUNT])
        matched = (
            np.asarray(validation[_FASTA_MATCHED_KEY])
            if _FASTA_MATCHED_KEY in validation.columns
            else protein_counts > 0
        )
        return {
            "feature_count": int(target.n_vars),
            "matched_feature_count": int(np.count_nonzero(matched)),
            "proteotypic_feature_count": int(np.count_nonzero(protein_counts == 1)),
        }

    annotation = target.varm.get(_FASTA_ANNOTATION_KEY)
    if annotation is None:
        return None
    annotated = np.asarray(annotation.notna().any(axis=1))
    return {
        "feature_count": int(target.n_vars),
        "annotated_feature_count": int(np.count_nonzero(annotated)),
    }


def describe(obj: Any) -> dict[str, Any]:
    """Return the descriptive summary for an in-memory AnnData or MuData.

    Newly converted objects use their stored stage-owned components. For legacy
    AnnData objects without a stored quantification component, that component is
    computed from the layers in memory without mutating the object.
    """
    if _is_mudata(obj):
        payload = _read_payload(obj)
        payload.setdefault("schema_version", _SCHEMA_VERSION)
        payload.setdefault("container_type", "mudata")
        payload.setdefault(
            "quantification",
            {
                "n_runs": int(obj.n_obs),
                "n_features": int(obj.n_vars),
                "modalities": list(obj.mod),
            },
        )
        payload["modalities"] = {name: describe(modality) for name, modality in obj.mod.items()}
        return payload

    payload = _read_payload(obj)
    payload.setdefault("schema_version", _SCHEMA_VERSION)
    payload.setdefault("container_type", "anndata")
    if "quantification" not in payload:
        payload["quantification"] = _quantification_summary(obj)
    if "column_mapping" not in payload:
        column_mapping = _column_mapping(obj)
        if column_mapping is not None:
            payload["column_mapping"] = column_mapping
    proteobench = obj.uns.get("proteobench")
    if isinstance(proteobench, Mapping) and "scores" in proteobench:
        payload["proteobench"] = _to_json_compatible(proteobench)
    return payload


def describe_path(path: Path | str, modality: str | None = None) -> dict[str, Any]:
    """Read and describe one converted container or one MuData modality.

    Args:
        path: Path to an APB ``.h5ad`` or ``.h5mu`` container.
        modality: Optional modality name. Valid only for MuData inputs.

    Returns:
        A JSON-compatible descriptive-summary dictionary.

    Raises:
        ValueError: The suffix is unsupported or a modality target is invalid.
    """
    result_path = Path(path).expanduser()
    if result_path.suffix == ".h5ad":
        if modality is not None:
            raise ValueError("--modality applies only to MuData (.h5mu) inputs")
        import anndata as ad

        obj = ad.read_h5ad(result_path, backed="r")
        try:
            return describe(obj)
        finally:
            obj.file.close()

    if result_path.suffix == ".h5mu":
        import mudata

        with mudata.set_options(pull_on_update=False):
            obj = mudata.read_h5mu(result_path, backed="r")
            try:
                if modality is None:
                    return describe(obj)
                if modality not in obj.mod:
                    raise ValueError(
                        f"modality {modality!r} not in MuData; modalities: {list(obj.mod)}"
                    )
                return describe(obj.mod[modality])
            finally:
                _close_mudata(obj)

    raise ValueError(f"unsupported converted result type: {result_path}")


def _store_quantification_summary_anndata(obj: Any) -> None:
    payload = _read_payload(obj)
    payload.update(
        {
            "schema_version": _SCHEMA_VERSION,
            "container_type": "anndata",
            "quantification": _quantification_summary(obj),
        }
    )
    column_mapping = _column_mapping(obj)
    if column_mapping is None:
        payload.pop("column_mapping", None)
    else:
        payload["column_mapping"] = column_mapping
    _write_payload(obj, payload)


def _quantification_summary(obj: Any) -> dict[str, Any]:
    metadata = obj.uns.get(_NAMESPACE) or {}
    params = read_search_parameters(obj)
    return {
        "n_runs": int(obj.n_obs),
        "n_features": int(obj.n_vars),
        "level": metadata.get("quantification_level"),
        "software_name": metadata.get("software_name"),
        "software_version": params.software_version if params is not None else None,
        "layers": {
            str(name): _layer_summary(layer, n_obs=obj.n_obs)
            for name, layer in named_layers(obj).items()
        },
    }


def _column_mapping(obj: Any) -> dict[str, Any] | None:
    """Describe where the effective rule placed vendor and computed columns."""
    rule = _stored_rule(obj)
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


def _stored_rule(obj: Any) -> ParseRule | None:
    """Return a valid stored effective rule, or ``None`` for legacy metadata."""
    namespace = obj.uns.get(_NAMESPACE)
    if not isinstance(namespace, Mapping):
        return None
    raw = namespace.get("rule_json")
    try:
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8")
        if isinstance(raw, str):
            return ParseRule.model_validate_json(raw)
        if isinstance(raw, Mapping):
            return ParseRule.model_validate(dict(raw))
    except (TypeError, UnicodeDecodeError, ValueError):
        return None
    return None


def _column_group_mapping(group: ColumnGroup) -> dict[str, str]:
    """Map output column names to their vendor source or compute operation."""
    mapping = dict(group.select)
    mapping.update({column.name: f"computed:{column.how}" for column in group.compute})
    return mapping


def _layer_summary(layer: Any, *, n_obs: int) -> dict[str, Any]:
    values = _matrix_values(layer)
    finite = np.isfinite(values)
    present_per_feature = np.count_nonzero(finite, axis=0)
    missing_per_feature = n_obs - present_per_feature
    histogram = np.bincount(missing_per_feature, minlength=n_obs + 1)
    observed = values[finite]
    return {
        "missingness_histogram": {
            str(missing_runs): int(feature_count)
            for missing_runs, feature_count in enumerate(histogram)
        },
        "summary": _numeric_summary(observed),
    }


def _numeric_summary(values: np.ndarray) -> dict[str, float | None]:
    """Return the six statistics produced by R's summary() for numeric values."""
    if not values.size:
        return {
            "min": None,
            "first_quartile": None,
            "median": None,
            "mean": None,
            "third_quartile": None,
            "max": None,
        }
    first_quartile, median, third_quartile = np.quantile(
        values,
        (0.25, 0.5, 0.75),
        method="linear",
    )
    return {
        "min": float(np.min(values)),
        "first_quartile": float(first_quartile),
        "median": float(median),
        "mean": float(np.mean(values)),
        "third_quartile": float(third_quartile),
        "max": float(np.max(values)),
    }


def _matrix_values(layer: Any) -> np.ndarray:
    if hasattr(layer, "to_memory"):
        layer = layer.to_memory()
    if is_sparse_matrix(layer):
        return layer.toarray()
    return np.asarray(layer)


def _to_json_compatible(value: Any) -> Any:
    """Copy an HDF5-decoded mapping into ordinary JSON-compatible values."""
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


def _read_payload(obj: Any) -> dict[str, Any]:
    namespace = obj.uns.get(_NAMESPACE) or {}
    raw = namespace.get(_SUMMARY_KEY)
    if raw is None:
        return {}
    if isinstance(raw, bytes):
        raw = raw.decode("utf-8")
    if isinstance(raw, str):
        payload = json.loads(raw)
    elif isinstance(raw, Mapping):
        payload = deepcopy(dict(raw))
    else:
        raise ValueError(f"invalid stored descriptive summary: {type(raw).__name__}")
    version = payload.get("schema_version")
    if version == "1":
        payload = _upgrade_v1_payload(payload)
        version = payload["schema_version"]
    if version == "2":
        payload = _upgrade_v2_payload(payload)
        version = payload["schema_version"]
    if version == "3":
        payload = _upgrade_v3_payload(payload)
        version = payload["schema_version"]
    if version == "4":
        payload = _upgrade_v4_payload(payload)
        version = payload["schema_version"]
    if version != _SCHEMA_VERSION:
        raise ValueError(f"unsupported descriptive-summary schema version: {version!r}")
    return payload


def _upgrade_v1_payload(payload: dict[str, Any]) -> dict[str, Any]:
    """Convert present-run histogram lists into missing-run count mappings."""
    quantification = payload.get("quantification") or {}
    n_runs = quantification.get("n_runs")
    if isinstance(n_runs, int):
        for layer in (quantification.get("layers") or {}).values():
            histogram = layer.get("missingness_histogram")
            if isinstance(histogram, list):
                layer["missingness_histogram"] = {
                    str(missing_runs): int(histogram[n_runs - missing_runs])
                    for missing_runs in range(n_runs + 1)
                }
    payload["schema_version"] = "2"
    return payload


def _upgrade_v2_payload(payload: dict[str, Any]) -> dict[str, Any]:
    """Rename generic numeric statistics and reserve unavailable legacy hinges."""
    quantification = payload.get("quantification") or {}
    for layer in (quantification.get("layers") or {}).values():
        intensity = layer.pop("intensity", None)
        if isinstance(intensity, Mapping):
            layer["fivenum"] = {
                "min": intensity.get("min"),
                "lower_hinge": None,
                "median": intensity.get("median"),
                "upper_hinge": None,
                "max": intensity.get("max"),
            }
    payload["schema_version"] = "3"
    return payload


def _upgrade_v3_payload(payload: dict[str, Any]) -> dict[str, Any]:
    """Rename five-number statistics and reserve unavailable summary values."""
    quantification = payload.get("quantification") or {}
    for layer in (quantification.get("layers") or {}).values():
        fivenum = layer.pop("fivenum", None)
        if isinstance(fivenum, Mapping):
            layer["summary"] = {
                "min": fivenum.get("min"),
                "first_quartile": None,
                "median": fivenum.get("median"),
                "mean": None,
                "third_quartile": None,
                "max": fivenum.get("max"),
            }
    payload["schema_version"] = "4"
    return payload


def _upgrade_v4_payload(payload: dict[str, Any]) -> dict[str, Any]:
    """Advance legacy summaries; column mapping is unavailable without the object."""
    payload["schema_version"] = _SCHEMA_VERSION
    return payload


def _write_payload(obj: Any, payload: dict[str, Any]) -> None:
    obj.uns.setdefault(_NAMESPACE, {})
    obj.uns[_NAMESPACE][_SUMMARY_KEY] = json.dumps(
        payload,
        allow_nan=False,
        sort_keys=True,
    )


def _is_mudata(obj: Any) -> bool:
    return hasattr(obj, "mod")


def _close_mudata(obj: Any) -> None:
    file_manager = getattr(obj, "file", None)
    if file_manager is not None:
        file_manager.close()
        return
    for modality in obj.mod.values():
        modality.file.close()
