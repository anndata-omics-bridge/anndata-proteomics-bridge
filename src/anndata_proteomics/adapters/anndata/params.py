"""AnnData-family persistence for validated search parameters."""

from __future__ import annotations

import json
from dataclasses import dataclass

from anndata import AnnData
from mudata import MuData

from anndata_proteomics.adapters.anndata.namespace import (
    MissingNamespaceText,
    has_namespace_key,
    read_namespace_text,
    update_namespace,
)
from anndata_proteomics.params.model import Parameters

_PARAMETERS_KEY = "search_parameters"


@dataclass(frozen=True, slots=True)
class MissingSearchParameters:
    """No search parameters are stored on this container."""


MISSING_SEARCH_PARAMETERS = MissingSearchParameters()


def has_search_parameters(target: AnnData | MuData) -> bool:
    """Return whether the APB namespace contains a search-parameter payload."""
    return has_namespace_key(target, _PARAMETERS_KEY)


def read_search_parameters(target: AnnData | MuData) -> Parameters | MissingSearchParameters:
    """Read the stored search parameters, or report their absence, with one validation."""
    payload = read_namespace_text(target, _PARAMETERS_KEY)
    if isinstance(payload, MissingNamespaceText):
        return MISSING_SEARCH_PARAMETERS
    return Parameters.model_validate(json.loads(payload))


def require_search_parameters(target: AnnData | MuData) -> Parameters:
    """Read and validate the stored search parameters.

    Raises:
        ValueError: If no search parameters were stored.
        TypeError: If the stored payload is not JSON text.
        json.JSONDecodeError: If the stored payload is malformed JSON.
        pydantic.ValidationError: If the decoded payload violates the schema.
    """
    parameters = read_search_parameters(target)
    if isinstance(parameters, MissingSearchParameters):
        raise ValueError("converted object has no stored search parameters")
    return parameters


def write_search_parameters(target: AnnData | MuData, parameters: Parameters) -> None:
    """Serialize validated search parameters into the APB namespace."""
    update_namespace(
        target,
        {_PARAMETERS_KEY: json.dumps(parameters.model_dump(mode="json"))},
    )
