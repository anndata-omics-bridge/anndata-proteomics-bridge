"""Exact in-memory matrix types accepted by backend-neutral calculations."""

from __future__ import annotations

from typing import TypeIs

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import csc_array, csc_matrix, csr_array, csr_matrix

type DenseQuantMatrix = NDArray[np.float32] | NDArray[np.float64]
type SparseQuantMatrix = (
    csr_matrix[np.float32]
    | csr_matrix[np.float64]
    | csc_matrix[np.float32]
    | csc_matrix[np.float64]
    | csr_array[np.float32]
    | csr_array[np.float64]
    | csc_array[np.float32]
    | csc_array[np.float64]
)
type QuantMatrix = DenseQuantMatrix | SparseQuantMatrix

_SPARSE_TYPES = (csr_matrix, csc_matrix, csr_array, csc_array)


def is_sparse_matrix(value: QuantMatrix) -> TypeIs[SparseQuantMatrix]:
    """Narrow a quantitative matrix to the supported SciPy containers."""
    return isinstance(value, _SPARSE_TYPES)
