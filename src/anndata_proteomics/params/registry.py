"""Software-name → parameter-parser dispatch."""

from __future__ import annotations

import re
from collections.abc import Callable
from pathlib import Path

from anndata_proteomics.params.model import Parameters
from anndata_proteomics.params.parsers.alphadia import extract_params as _alphadia_extract
from anndata_proteomics.params.parsers.alphapept import extract_params as _alphapept_extract
from anndata_proteomics.params.parsers.diann import extract_params as _diann_extract
from anndata_proteomics.params.parsers.fragpipe import extract_params as _fragpipe_extract
from anndata_proteomics.params.parsers.maxquant import extract_params as _maxquant_extract
from anndata_proteomics.params.parsers.metamorpheus import extract_params as _metamorpheus_extract
from anndata_proteomics.params.parsers.msaid import extract_params as _msaid_extract
from anndata_proteomics.params.parsers.peaks import extract_params as _peaks_extract
from anndata_proteomics.params.parsers.sage import extract_params as _sage_extract
from anndata_proteomics.params.parsers.spectronaut import extract_params as _spectronaut_extract
from anndata_proteomics.params.parsers.wombat import extract_params as _wombat_extract

ParseFn = Callable[..., Parameters]


_REGISTRY: dict[str, ParseFn] = {
    "alphadia": _alphadia_extract,
    "alphapept": _alphapept_extract,
    "dia-nn": _diann_extract,
    "diann": _diann_extract,
    "fragpipe": _fragpipe_extract,
    "maxquant": _maxquant_extract,
    "metamorpheus": _metamorpheus_extract,
    "msaid": _msaid_extract,
    "peaks": _peaks_extract,
    "sage": _sage_extract,
    "spectronaut": _spectronaut_extract,
    "wombat": _wombat_extract,
}


def get_parser(software: str) -> ParseFn:
    """Look up a parser by software name (case-insensitive)."""
    key = software.lower()
    if key not in _REGISTRY:
        raise KeyError(
            f"no parameter parser registered for {software!r}; available: {sorted(_REGISTRY)}"
        )
    return _REGISTRY[key]


def parse_params(path: str | Path, software: str) -> Parameters:
    """Convenience: look up a parser and run it on ``path``."""
    return get_parser(software)(path)


def parser_slug(software_name: str) -> str | None:
    """Resolve the primary parameter parser named by a catalog software label.

    Compound labels place the workflow owner first, for example
    ``FragPipe (DIA-NN quant)``. The earliest registered software token is
    therefore the parameter-file owner, while result-table recognition remains
    the responsibility of the parsing-rule registry.
    """
    normalized_name = _normalize_software_name(software_name)
    canonical: dict[str, str] = {}
    for candidate in _REGISTRY:
        token = _normalize_software_name(candidate)
        current = canonical.get(token)
        if current is None or candidate == token:
            canonical[token] = candidate
    matches = [
        (normalized_name.find(token), -len(token), candidate)
        for token, candidate in canonical.items()
        if token and token in normalized_name
    ]
    return min(matches)[2] if matches else None


def available_software() -> list[str]:
    return sorted(_REGISTRY)


def _normalize_software_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.lower())
