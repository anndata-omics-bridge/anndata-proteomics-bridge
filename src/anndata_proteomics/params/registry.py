"""Software-name → parameter-parser dispatch."""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Union

from anndata_proteomics.params.alphapept import extract_params as _alphapept_extract
from anndata_proteomics.params.model import Parameters
from anndata_proteomics.params.sage import extract_params as _sage_extract
from anndata_proteomics.params.wombat import extract_params as _wombat_extract


ParseFn = Callable[[Union[str, Path]], Parameters]


_REGISTRY: dict[str, ParseFn] = {
    "alphapept": _alphapept_extract,
    "sage": _sage_extract,
    "wombat": _wombat_extract,
}


def get_parser(software: str) -> ParseFn:
    """Look up a parser by software name (case-insensitive)."""
    key = software.lower()
    if key not in _REGISTRY:
        raise KeyError(
            f"no parameter parser registered for {software!r}; "
            f"available: {sorted(_REGISTRY)}"
        )
    return _REGISTRY[key]


def parse_params(path: Union[str, Path], software: str) -> Parameters:
    """Convenience: look up a parser and run it on ``path``."""
    return get_parser(software)(path)


def available_software() -> list[str]:
    return sorted(_REGISTRY)
