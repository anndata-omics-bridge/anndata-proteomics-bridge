"""Structural types for sparse matrix operations used across APB."""

from __future__ import annotations

from typing import Protocol, TypeIs

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
