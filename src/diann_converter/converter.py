"""
Core converter: Proteomics data to AnnData format.

Converts DIA-NN, Spectronaut, and other proteomics software outputs to AnnData h5ad files.
"""

import numpy as np
import pandas as pd
from anndata import AnnData
from typing import Optional, List, Dict, Union

from .reader import load_file
from .annotation import match_samples, load_annotation
from .parse_settings import load_config


def proteomics_to_anndata(
    data_path: str,
    annotation_path: str,
    software: str = None,
    annotation_id_col: str = "sample",
    factor_cols: Optional[List[str]] = None,
    label_cols: Optional[List[str]] = None,
    run_col: str = "run",
    qvalue_threshold: Optional[float] = 0.01,
    log2_transform: bool = True,
) -> AnnData:
    """
    Convert proteomics data to AnnData format.

    Parameters
    ----------
    data_path : str
        Path to proteomics data file (DIA-NN report.tsv, Spectronaut export, etc.).
    annotation_path : str
        Path to sample annotation CSV/TSV file.
    software : str, optional
        Software name ("diann", "spectronaut"). If None, auto-detects.
    annotation_id_col : str, optional
        Column name in annotation file containing sample IDs. Defaults to "sample".
    factor_cols : list of str, optional
        Annotation columns to use as factors (experimental conditions).
        If None, auto-detects columns with categorical data.
    label_cols : list of str, optional
        Annotation columns to use as labels (sample identifiers).
        If None, uses index.
    run_col : str, optional
        Column name for run/sample identifier in data file. Defaults to "run".
    qvalue_threshold : float, optional
        Q-value threshold for filtering precursors. Defaults to 0.01 (1% FDR).
        If None, no filtering.
    log2_transform : bool, optional
        Whether to log2-transform intensities in X. Defaults to True.

    Returns
    -------
    AnnData
        AnnData object with:
        - X: Log2-transformed intensities (precursors × samples)
        - layers['intensities']: Raw intensities
        - obs: Sample annotations
        - var: Precursor metadata (sequence, charge, proteins, genes)
        - uns['exploreDE']: Metadata for exploreDE compatibility

    Examples
    --------
    >>> adata = proteomics_to_anndata(
    ...     "report.tsv",
    ...     "annotation.csv",
    ...     software="diann",
    ...     factor_cols=["condition", "batch"]
    ... )
    >>> adata.write_h5ad("proteomics.h5ad")
    """
    # Load configuration
    config = load_config(software) if software else None

    # Load data
    df = load_file(data_path, software=software)

    # Load annotation
    annotation_df = load_annotation(annotation_path, id_col=annotation_id_col)

    # Apply column mapping
    if config:
        df = config.rename_columns(df)
        df = config.convert_modifications(df)
        df = config.filter_decoys(df)
        df = config.mark_contaminants(df)

    # Filter by Q-value
    if qvalue_threshold is not None:
        qvalue_cols = [col for col in ["qvalue", "global_qvalue", "protein_qvalue"] if col in df.columns]
        if qvalue_cols:
            qval_col = qvalue_cols[0]  # Use first available
            n_before = len(df)
            df = df[df[qval_col] <= qvalue_threshold]
            n_after = len(df)
            print(f"Filtered by {qval_col} <= {qvalue_threshold}: {n_before} -> {n_after} precursors")

    # Create precursor ID
    if "proforma" in df.columns and "charge" in df.columns:
        df["precursor_id"] = df["proforma"] + "/" + df["charge"].astype(str)
    elif "sequence" in df.columns and "charge" in df.columns:
        df["precursor_id"] = df["sequence"] + "/" + df["charge"].astype(str)
    else:
        raise ValueError("Cannot create precursor ID: missing sequence/proforma and/or charge columns")

    # Identify run column
    if run_col not in df.columns:
        # Try alternatives
        for alt in ["run", "file_name", "raw_file", "Run", "File.Name"]:
            if alt in df.columns:
                run_col = alt
                break
        else:
            raise ValueError(f"Run column '{run_col}' not found. Available columns: {list(df.columns)}")

    # Match samples
    unique_runs = df[run_col].unique().tolist()
    matched_annotation = match_samples(unique_runs, annotation_df, annotation_id_col=annotation_id_col)

    # Pivot to wide format: precursors (rows) × samples (columns)
    intensity_col = "intensity" if "intensity" in df.columns else "quantity"
    if intensity_col not in df.columns:
        raise ValueError(f"Intensity column not found. Available columns: {list(df.columns)}")

    # Create pivot table
    pivot_df = df.pivot_table(
        index="precursor_id",
        columns=run_col,
        values=intensity_col,
        aggfunc="first"  # If duplicates, take first
    )

    # Reorder columns to match annotation
    pivot_df = pivot_df[matched_annotation.index]

    # Create var (precursor metadata)
    var_df = df.groupby("precursor_id").first()[
        [col for col in ["proforma", "sequence", "charge", "proteins", "genes"] if col in df.columns]
    ]

    # Add description column (required by exploreDE)
    if "proteins" in var_df.columns:
        var_df["description"] = var_df["proteins"]
    else:
        var_df["description"] = var_df.index

    # Ensure var matches pivot index
    var_df = var_df.loc[pivot_df.index]

    # Create obs (sample metadata)
    obs_df = matched_annotation.copy()

    # Auto-detect factor and label columns if not specified
    if factor_cols is None:
        # Use categorical columns as factors
        factor_cols = [col for col in obs_df.columns if obs_df[col].dtype == "object" or obs_df[col].dtype.name == "category"]

    if label_cols is None:
        # Use index as label
        label_cols = []

    # Create intensity matrix (precursors × samples)
    # Keep NaN for missing values - DIA-NN only reports detected precursors
    X_raw = pivot_df.values.astype(np.float64)

    # Log2 transform for X
    if log2_transform:
        X = np.log2(X_raw + 1)  # Add pseudocount to avoid log(0)
    else:
        X = X_raw.copy()

    # Create AnnData object (transpose to match AnnData convention: samples × features)
    adata = AnnData(
        X=X.T,  # Transpose: samples × precursors
        obs=obs_df,
        var=var_df,
    )

    # Add layers
    adata.layers["intensities"] = X_raw.T

    # Add exploreDE metadata
    adata.uns["exploreDE"] = {
        "column_roles": {
            "var": {
                "description": ["description"],
                "label": [col for col in ["proforma", "sequence", "proteins", "genes"] if col in var_df.columns],
            },
            "obs": {
                "factor": factor_cols,
                "label": label_cols if label_cols else ["sample"],
            },
        }
    }

    return adata


def convert_file(
    data_path: str,
    annotation_path: str,
    output_path: str,
    **kwargs
) -> None:
    """
    Convert proteomics file to AnnData h5ad format.

    Convenience function that calls proteomics_to_anndata() and saves the result.

    Parameters
    ----------
    data_path : str
        Path to proteomics data file.
    annotation_path : str
        Path to sample annotation file.
    output_path : str
        Path for output h5ad file.
    **kwargs
        Additional arguments passed to proteomics_to_anndata().

    Examples
    --------
    >>> convert_file(
    ...     "report.tsv",
    ...     "annotation.csv",
    ...     "output.h5ad",
    ...     software="diann",
    ...     factor_cols=["condition"]
    ... )
    """
    adata = proteomics_to_anndata(data_path, annotation_path, **kwargs)
    adata.write_h5ad(output_path)
    print(f"Saved AnnData object to {output_path}")
    print(f"Shape: {adata.shape} (samples × precursors)")
    print(f"Factors: {adata.uns['exploreDE']['column_roles']['obs']['factor']}")
