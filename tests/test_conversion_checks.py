"""Layer-occupancy conversion contract checks."""

from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse

from anndata_proteomics.converters.checks import (
    LayerContractError,
    check_layer_occupancy,
    layer_occupancies,
)


def _full(value: float = 1.0) -> np.ndarray:
    return np.full((4, 5), value, dtype=np.float64)


def _empty() -> np.ndarray:
    return np.full((4, 5), np.nan, dtype=np.float64)


def test_occupancy_counts_only_finite_cells() -> None:
    matrix = _empty()
    matrix[0, 0] = 1.0
    matrix[1, 1] = np.inf

    (occupancy,) = layer_occupancies({"Quantity": matrix})

    assert (occupancy.present, occupancy.total) == (1, 20)
    assert occupancy.ratio == pytest.approx(0.05)
    assert "Quantity: 1/20 (5.00%)" == occupancy.describe()


def test_sparse_layer_is_measured_without_densifying() -> None:
    dense = np.zeros((4, 5), dtype=np.float64)
    dense[0, 0] = 3.0
    dense[2, 3] = 7.0

    (occupancy,) = layer_occupancies({"Quantity": sparse.csr_matrix(dense)})

    assert (occupancy.present, occupancy.total) == (2, 20)


def test_empty_x_layer_beside_a_populated_sibling_is_an_error() -> None:
    # The parse-failure fingerprint: the source column was read but nothing survived
    # coercion, while a sibling column from the same file is full.
    with pytest.raises(LayerContractError, match="PG_Quantity.*effectively empty"):
        check_layer_occupancy(
            {"PG_Quantity": _empty(), "PG_RunEvidenceCount": _full()},
            x_layer="PG_Quantity",
        )


def test_other_empty_layers_warn_and_only_fail_under_strict() -> None:
    layers = {"PG_Quantity": _full(), "PG_Cscore": _empty()}

    occupancies = check_layer_occupancy(layers, x_layer="PG_Quantity")
    assert [item.name for item in occupancies] == ["PG_Quantity", "PG_Cscore"]

    with pytest.raises(LayerContractError, match="PG_Cscore"):
        check_layer_occupancy(layers, x_layer="PG_Quantity", strict=True)


def test_uniformly_sparse_layers_are_accepted() -> None:
    # Single-cell acquisitions are legitimately sparse in every layer at once; only a
    # sparse layer *beside a populated sibling* indicates lost data.
    sparse_layer = _empty()
    sparse_layer[0, 0] = 1.0

    occupancies = check_layer_occupancy(
        {"Quantity": sparse_layer, "Score": sparse_layer.copy()},
        x_layer="Quantity",
        strict=True,
    )

    assert all(item.present == 1 for item in occupancies)


def test_a_single_populated_layer_is_not_compared_against_itself() -> None:
    check_layer_occupancy({"Quantity": _full()}, x_layer="Quantity", strict=True)
