"""Read and write :class:`Parameters` to ``AnnData.uns``.

Storage layer: a JSON string under
``uns['anndata_proteomics']['search_parameters']`` (h5py can't serialize
heterogeneous nested dicts, so JSON is the wire format).

API layer: ``read_search_parameters`` returns a validated
:class:`Parameters` instance; ``write_search_parameters`` round-trips one
back. This mirrors how :class:`ParseRule` is the typed handle for the
TOML side.
"""

from __future__ import annotations

import json

import anndata as ad

from anndata_proteomics.params.model import Parameters

_UNS_KEY = "anndata_proteomics"
_PARAMS_KEY = "search_parameters"
_PARAMS_PATH_KEY = "search_parameters_path"


def read_search_parameters(adata: ad.AnnData) -> Parameters | None:
    """Return the stored :class:`Parameters`, or ``None`` if none were stored.

    Raises :class:`pydantic.ValidationError` if the stored payload no longer
    fits the current :class:`Parameters` schema.
    """
    container = adata.uns.get(_UNS_KEY)
    if not container:
        return None
    raw = container.get(_PARAMS_KEY)
    if not raw:
        return None
    return Parameters.model_validate(json.loads(raw))


def write_search_parameters(
    adata: ad.AnnData,
    params: Parameters,
    *,
    source_path: str | None = None,
) -> None:
    """Serialize ``params`` into ``adata.uns`` (and the original source path)."""
    adata.uns.setdefault(_UNS_KEY, {})
    adata.uns[_UNS_KEY][_PARAMS_KEY] = json.dumps(params.model_dump(mode="json"))
    if source_path is not None:
        adata.uns[_UNS_KEY][_PARAMS_PATH_KEY] = str(source_path)


def get_search_parameters_path(adata: ad.AnnData) -> str | None:
    """Return the original source-path provenance, if recorded."""
    container = adata.uns.get(_UNS_KEY)
    if not container:
        return None
    value = container.get(_PARAMS_PATH_KEY)
    return str(value) if value is not None else None
