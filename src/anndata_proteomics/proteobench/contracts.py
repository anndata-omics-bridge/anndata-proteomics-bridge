"""Typed in-memory contracts shared by ProteoBench calculations."""

from __future__ import annotations

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
type FloatArray = NDArray[np.float32] | NDArray[np.float64]
type FloatDType = type[np.float32] | type[np.float64]
