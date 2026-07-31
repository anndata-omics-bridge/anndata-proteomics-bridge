"""Concrete AnnData extraction for matrix metadata."""

from __future__ import annotations

from anndata import AnnData


def layer_names(target: AnnData) -> tuple[str, ...]:
    """Return the names of explicitly stored quantitative layers."""
    return tuple(str(name) for name in target.layers.keys() if name is not None)
