"""Structural types for sparse matrix operations used across APB."""

from __future__ import annotations

from typing import Any, Protocol, TypeIs

import numpy as np
from numpy.typing import NDArray
from scipy import sparse


class SparseMatrix(Protocol):
    """Operations shared by SciPy sparse matrices and sparse arrays."""

    # Element type is intentionally unconstrained: SciPy's containers are generic in it,
    # and a mutable attribute is invariant, so pinning it here would reject csr_matrix[float64].
    data: NDArray[Any]

    @property
    def shape(self) -> tuple[int, ...]:
        """Matrix dimensions, ``(n_obs, n_vars)`` for a quantitative layer."""
        ...

    def toarray(self) -> NDArray[np.generic]:
        """Materialize the sparse values as a dense NumPy array."""
        ...


type QuantMatrix = NDArray[np.floating[Any]] | SparseMatrix
"""A quantitative layer as stored: dense NumPy values or a SciPy sparse container.

Slicing is deliberately absent from :class:`SparseMatrix`: SciPy's ``__getitem__`` is
declared over its own shape type variables and no structural signature accepts every
container, so the two helpers that slice a layer by row take an unconstrained matrix.
"""


def is_sparse_matrix(value: object) -> TypeIs[SparseMatrix]:
    """Narrow values recognized by SciPy to their shared sparse operations."""
    return sparse.issparse(value)


def named_layers(adata: Any) -> dict[str, Any]:
    """Return only the explicitly named layers of an AnnData.

    From anndata 0.13, ``adata.layers`` also yields ``X`` under a ``None`` key. APB
    always writes ``X`` as a copy of ``axis.x_layer``, so that entry is a duplicate
    under a name no rule declares; iterating it raw leaks a ``"None"`` layer into
    summaries and logs.

    Structurally typed: ``readers.summary`` passes its own lightweight metadata record
    here rather than a real AnnData, and ``AnnData.layers`` is a plain attribute that no
    read-only protocol can describe without excluding one of the two callers.
    """
    return {name: layer for name, layer in adata.layers.items() if name is not None}
