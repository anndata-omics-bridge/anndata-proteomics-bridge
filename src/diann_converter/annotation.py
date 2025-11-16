"""
Sample annotation handling and matching.

Matches sample identifiers from proteomics files to user-provided annotation.
"""

import pandas as pd
from typing import Optional, List
import re


def clean_sample_name(sample_name: str) -> str:
    """
    Clean sample name by removing common file extensions and paths.

    Parameters
    ----------
    sample_name : str
        Sample name to clean.

    Returns
    -------
    str
        Cleaned sample name.

    Examples
    --------
    >>> clean_sample_name("/path/to/sample.mzML")
    'sample'
    >>> clean_sample_name("sample.raw")
    'sample'
    >>> clean_sample_name("sample.d")
    'sample'
    """
    # Remove path
    name = sample_name.split("/")[-1]
    name = name.split("\\")[-1]

    # Remove common extensions
    extensions = [
        ".mzML", ".mzML.gz", ".mzXML", ".raw", ".d",
        ".wiff", ".wiff2", ".RAW", ".MZML", ".MZXML"
    ]
    for ext in extensions:
        if name.endswith(ext):
            name = name[:-len(ext)]

    return name


def match_samples(
    data_samples: List[str],
    annotation_df: pd.DataFrame,
    annotation_id_col: str = "sample",
    fuzzy: bool = True
) -> pd.DataFrame:
    """
    Match data samples to annotation file.

    Parameters
    ----------
    data_samples : list of str
        Sample names from proteomics data.
    annotation_df : pd.DataFrame
        Annotation dataframe with sample metadata.
    annotation_id_col : str, optional
        Column name in annotation_df containing sample identifiers. Defaults to "sample".
    fuzzy : bool, optional
        If True, tries fuzzy matching (cleaning filenames). Defaults to True.

    Returns
    -------
    pd.DataFrame
        Annotation dataframe indexed by matched sample names.

    Raises
    ------
    ValueError
        If samples cannot be matched or annotation_id_col not found.

    Examples
    --------
    >>> data_samples = ["sample1.mzML", "sample2.mzML"]
    >>> annotation = pd.DataFrame({
    ...     "sample": ["sample1", "sample2"],
    ...     "condition": ["A", "B"]
    ... })
    >>> matched = match_samples(data_samples, annotation)
    """
    if annotation_id_col not in annotation_df.columns:
        raise ValueError(
            f"Annotation ID column '{annotation_id_col}' not found. "
            f"Available columns: {list(annotation_df.columns)}"
        )

    # Create mapping dictionary from annotation
    annotation_dict = annotation_df.set_index(annotation_id_col).to_dict(orient="index")
    annotation_keys = set(annotation_dict.keys())

    # Try exact matching first
    data_samples_set = set(data_samples)
    if data_samples_set.issubset(annotation_keys):
        # Perfect match
        matched_annotation = annotation_df.set_index(annotation_id_col).loc[data_samples]
        matched_annotation.index.name = "sample"
        return matched_annotation

    if not fuzzy:
        missing = data_samples_set - annotation_keys
        raise ValueError(
            f"Could not match samples. Missing samples in annotation: {missing}"
        )

    # Try fuzzy matching (clean filenames)
    cleaned_data = {clean_sample_name(s): s for s in data_samples}
    cleaned_annotation = {clean_sample_name(s): s for s in annotation_keys}

    # Build mapping: data_sample -> annotation_sample
    sample_mapping = {}
    unmatched = []

    for cleaned_data_name, original_data_name in cleaned_data.items():
        if cleaned_data_name in cleaned_annotation:
            annotation_name = cleaned_annotation[cleaned_data_name]
            sample_mapping[original_data_name] = annotation_name
        else:
            unmatched.append(original_data_name)

    if unmatched:
        raise ValueError(
            f"Could not match {len(unmatched)} samples to annotation:\n"
            f"Unmatched samples: {unmatched}\n"
            f"Available annotation samples: {list(annotation_keys)}\n"
            f"Hint: Ensure sample names match between data and annotation file."
        )

    # Create matched annotation dataframe
    matched_rows = []
    for data_sample in data_samples:
        annotation_sample = sample_mapping[data_sample]
        row = annotation_dict[annotation_sample].copy()
        matched_rows.append(row)

    matched_df = pd.DataFrame(matched_rows, index=data_samples)
    matched_df.index.name = "sample"

    return matched_df


def validate_annotation(
    annotation_df: pd.DataFrame,
    required_factor_cols: Optional[List[str]] = None,
    required_label_cols: Optional[List[str]] = None
) -> None:
    """
    Validate annotation dataframe structure.

    Parameters
    ----------
    annotation_df : pd.DataFrame
        Annotation dataframe to validate.
    required_factor_cols : list of str, optional
        Required factor columns (experimental conditions). If None, no validation.
    required_label_cols : list of str, optional
        Required label columns (sample identifiers). If None, no validation.

    Raises
    ------
    ValueError
        If validation fails.

    Examples
    --------
    >>> annotation = pd.DataFrame({
    ...     "sample": ["s1", "s2"],
    ...     "condition": ["A", "B"]
    ... })
    >>> validate_annotation(annotation, required_factor_cols=["condition"])
    """
    if required_factor_cols:
        missing_factors = set(required_factor_cols) - set(annotation_df.columns)
        if missing_factors:
            raise ValueError(
                f"Missing required factor columns: {missing_factors}. "
                f"Available columns: {list(annotation_df.columns)}"
            )

    if required_label_cols:
        missing_labels = set(required_label_cols) - set(annotation_df.columns)
        if missing_labels:
            raise ValueError(
                f"Missing required label columns: {missing_labels}. "
                f"Available columns: {list(annotation_df.columns)}"
            )


def load_annotation(
    annotation_path: str,
    id_col: str = "sample",
    **read_csv_kwargs
) -> pd.DataFrame:
    """
    Load annotation file from CSV/TSV.

    Parameters
    ----------
    annotation_path : str
        Path to annotation file (.csv or .tsv).
    id_col : str, optional
        Column name containing sample identifiers. Defaults to "sample".
    **read_csv_kwargs
        Additional arguments passed to pd.read_csv().

    Returns
    -------
    pd.DataFrame
        Loaded annotation dataframe.

    Raises
    ------
    FileNotFoundError
        If annotation file doesn't exist.
    ValueError
        If id_col not found in file.

    Examples
    --------
    >>> annotation = load_annotation("samples.csv")
    >>> annotation = load_annotation("samples.tsv", id_col="run_id")
    """
    # Auto-detect separator
    if annotation_path.endswith(".tsv"):
        sep = "\t"
    else:
        sep = ","

    kwargs = {"sep": sep, **read_csv_kwargs}
    annotation_df = pd.read_csv(annotation_path, **kwargs)

    if id_col not in annotation_df.columns:
        raise ValueError(
            f"ID column '{id_col}' not found in annotation file. "
            f"Available columns: {list(annotation_df.columns)}"
        )

    return annotation_df
