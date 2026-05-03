"""Tests for converters/assemble.py."""

from __future__ import annotations

import numpy as np
import pandas as pd

from anndata_proteomics.converters._pieces import ConversionPieces
from anndata_proteomics.converters.assemble import to_anndata
from anndata_proteomics.rules.loader import load_packaged_rule


def test_to_anndata_shape_and_uns() -> None:
    rule = load_packaged_rule("diann", "ion")
    obs = pd.DataFrame({"Run": ["S1", "S2"]}, index=["S1", "S2"])
    var = pd.DataFrame({"Modified_Sequence": ["P1", "P2", "P3"]}, index=["P1_2", "P2_2", "P3_3"])
    X = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    layers = {
        rule.axis.x_layer: X,
        "Q_Value": np.array([[0.01, 0.02, 0.03], [0.04, 0.05, 0.06]]),
    }
    pieces = ConversionPieces(X=X, obs=obs, var=var, layers=layers)
    adata = to_anndata(pieces, rule)

    assert adata.shape == (2, 3)
    assert "anndata_proteomics" in adata.uns
    assert adata.uns["anndata_proteomics"]["software_name"] == "DIA-NN"
    assert adata.uns["anndata_proteomics"]["input_shape"] == "long"
    assert "Q_Value" in adata.layers
