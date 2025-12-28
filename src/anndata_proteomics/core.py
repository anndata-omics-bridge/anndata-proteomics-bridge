"""
Core Converter class with strategy pattern.

Usage:
    converter = Converter(MaxQuantStrategy())
    adata = converter.convert("evidence.txt", "annotation.csv")

    # Or use the builder for auto-detection:
    converter = ConverterBuilder.from_file("evidence.txt")
    adata = converter.convert("evidence.txt", "annotation.csv")
"""

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData

from .annotation import load_annotation, match_samples


class Converter:
    """
    Proteomics to AnnData converter.

    Initialized with a strategy that handles format-specific loading.
    Uses strategy's get_obs(), get_var(), get_layers() methods for data extraction.
    """

    def __init__(self, strategy):
        """
        Initialize converter with a strategy.

        Parameters
        ----------
        strategy : object
            Strategy with obs_id, var_id attributes and get_obs(), get_var(),
            get_layers() methods.
        """
        self.strategy = strategy

    def convert(
        self,
        data_path: str | Path,
        annotation_path: str | Path,
        annotation_id_col: str = "sample",
        factor_cols: Optional[list[str]] = None,
        label_cols: Optional[list[str]] = None,
        x_layer: Optional[str] = None,
        log2_transform: bool = True,
    ) -> AnnData:
        """
        Convert proteomics data to AnnData format.

        Parameters
        ----------
        data_path : str or Path
            Path to proteomics data file.
        annotation_path : str or Path
            Path to sample annotation CSV/TSV file.
        annotation_id_col : str
            Column in annotation file containing sample IDs.
        factor_cols : list of str, optional
            Annotation columns to use as factors. Auto-detects if None.
        label_cols : list of str, optional
            Annotation columns to use as labels.
        x_layer : str, optional
            Which layer column to use for X matrix. Default is first layer column.
        log2_transform : bool
            Whether to log2-transform intensities in X.

        Returns
        -------
        AnnData
            AnnData object with samples as obs and precursors as var.
        """
        data_path = Path(data_path)
        annotation_path = Path(annotation_path)

        # Load data using strategy
        df = self.strategy.load(data_path)

        # Get obs_id and var_id from strategy
        obs_id = self.strategy.obs_id
        var_id = self.strategy.var_id

        # Get var metadata from strategy
        var_df = self.strategy.get_var(df)

        # Get layers data from strategy
        layers_df = self.strategy.get_layers(df)

        # Identify layer columns (everything except obs_id and var_id)
        layer_cols = [c for c in layers_df.columns if c not in {obs_id, var_id}]

        if not layer_cols:
            raise ValueError(f"No layer columns found. Available: {list(layers_df.columns)}")

        # Pivot each layer column to wide format
        layers = {}
        for col in layer_cols:
            pivot = layers_df.pivot_table(
                index=var_id,
                columns=obs_id,
                values=col,
                aggfunc="first"
            )
            layers[col] = pivot

        # Get sample order from first layer
        first_layer = layers[layer_cols[0]]
        sample_order = list(first_layer.columns)
        var_order = list(first_layer.index)

        # Load annotation and match to samples
        annotation_df = load_annotation(str(annotation_path), id_col=annotation_id_col)
        obs_df = match_samples(sample_order, annotation_df, annotation_id_col=annotation_id_col)

        # Ensure var matches layer order
        var_df = var_df.loc[var_order]

        # Auto-detect factor columns
        if factor_cols is None:
            factor_cols = [
                col for col in obs_df.columns
                if obs_df[col].dtype == "object" or obs_df[col].dtype.name == "category"
            ]

        label_cols = label_cols or []

        # Build numpy arrays for layers
        layers_arrays = {}
        for col in layer_cols:
            # Reorder columns to match obs_df order
            layer_pivot = layers[col][obs_df.index]
            layer_data = layer_pivot.values.astype(np.float64)
            # Transpose: layers should be (n_obs, n_var)
            layers_arrays[col] = layer_data.T

        # X is first layer (or user-specified)
        x_layer_name = x_layer or layer_cols[0]
        if x_layer_name not in layers_arrays:
            raise ValueError(f"Layer '{x_layer_name}' not found. Available: {list(layers_arrays.keys())}")

        X_raw = layers_arrays[x_layer_name]

        # Log2 transform for X
        if log2_transform:
            X = np.log2(X_raw + 1)
        else:
            X = X_raw.copy()

        # Create AnnData
        adata = AnnData(
            X=X,
            obs=obs_df,
            var=var_df,
        )

        # Add layers
        for col, data in layers_arrays.items():
            adata.layers[col] = data

        # Add metadata
        adata.uns["exploreDE"] = {
            "column_roles": {
                "var": {
                    "description": ["description"] if "description" in var_df.columns else [],
                    "label": [c for c in ["proforma", "stripped_sequence"] if c in var_df.columns],
                },
                "obs": {
                    "factor": factor_cols,
                    "label": label_cols or ["sample"],
                },
            },
            "quantification_level": "ion",
            "software": getattr(self.strategy, "name", "Unknown"),
            "layers": {col: f"Layer: {col}" for col in layer_cols},
            "primary_quantification": x_layer_name,
        }

        return adata

    def __repr__(self) -> str:
        return f"Converter(strategy={self.strategy!r})"
