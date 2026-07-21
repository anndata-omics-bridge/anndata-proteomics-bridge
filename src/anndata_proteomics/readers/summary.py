"""Stage-owned descriptive summaries for converted AnnData and MuData objects."""

from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse

from anndata_proteomics.params.anndata_io import read_search_parameters

_NAMESPACE = "anndata_proteomics"
_SUMMARY_KEY = "descriptive_summary"
_SCHEMA_VERSION = "1"
_FASTA_VALIDATION_KEY = "fasta_validation"
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


def store_fasta_summary(obj: Any) -> None:
    """Merge the FASTA-owned summary component when validation results exist."""
    targets = obj.mod.values() if _is_mudata(obj) else (obj,)
    for target in targets:
        validation = target.varm.get(_FASTA_VALIDATION_KEY)
        if validation is None or _FASTA_PROTEIN_COUNT not in validation.columns:
            continue
        protein_counts = np.asarray(validation[_FASTA_PROTEIN_COUNT])
        payload = _read_payload(target)
        payload.setdefault("schema_version", _SCHEMA_VERSION)
        payload.setdefault("container_type", "anndata")
        payload["fasta"] = {"proteotypic_feature_count": int(np.count_nonzero(protein_counts == 1))}
        _write_payload(target, payload)


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
            str(name): _layer_summary(layer, n_obs=obj.n_obs) for name, layer in obj.layers.items()
        },
    }


def _layer_summary(layer: Any, *, n_obs: int) -> dict[str, Any]:
    values = _matrix_values(layer)
    finite = np.isfinite(values)
    present_per_feature = np.count_nonzero(finite, axis=0)
    histogram = np.bincount(present_per_feature, minlength=n_obs + 1)
    observed = values[finite]
    intensity = (
        {"min": None, "median": None, "max": None}
        if not observed.size
        else {
            "min": float(np.min(observed)),
            "median": float(np.median(observed)),
            "max": float(np.max(observed)),
        }
    )
    return {
        "missingness_histogram": histogram.astype(int).tolist(),
        "intensity": intensity,
    }


def _matrix_values(layer: Any) -> np.ndarray:
    if hasattr(layer, "to_memory"):
        layer = layer.to_memory()
    if sparse.issparse(layer):
        return layer.toarray()
    return np.asarray(layer)


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
        payload = dict(raw)
    else:
        raise ValueError(f"invalid stored descriptive summary: {type(raw).__name__}")
    version = payload.get("schema_version")
    if version != _SCHEMA_VERSION:
        raise ValueError(f"unsupported descriptive-summary schema version: {version!r}")
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
