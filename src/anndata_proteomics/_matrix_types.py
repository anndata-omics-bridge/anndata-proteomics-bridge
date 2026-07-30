"""Structural types for sparse matrix operations used across APB."""

from __future__ import annotations

from typing import Any, Protocol, TypeIs

import numpy as np
from numpy.typing import NDArray
from scipy import sparse


class SparseMatrix(Protocol):
    """Operations shared by SciPy sparse matrices and sparse arrays."""

    data: NDArray[np.generic]

    def toarray(self) -> NDArray[np.generic]:
        """Materialize the sparse values as a dense NumPy array."""
        ...


def is_sparse_matrix(value: object) -> TypeIs[SparseMatrix]:
    """Narrow values recognized by SciPy to their shared sparse operations."""
    return sparse.issparse(value)


def named_layers(adata: Any) -> dict[str, Any]:
    """Return only the explicitly named layers of an AnnData.

    From anndata 0.13, ``adata.layers`` also yields ``X`` under a ``None`` key. APB
    always writes ``X`` as a copy of ``axis.x_layer``, so that entry is a duplicate
    under a name no rule declares; iterating it raw leaks a ``"None"`` layer into
    summaries and logs.
    """
    return {name: layer for name, layer in adata.layers.items() if name is not None}
