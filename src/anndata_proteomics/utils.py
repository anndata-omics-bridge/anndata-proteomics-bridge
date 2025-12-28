"""Utility functions for column name normalization."""

import re


def normalize_column_name(col: str) -> str:
    """
    Normalize column name to a valid Python identifier.

    Transforms column names to lowercase with underscores, preserving
    the original semantic meaning.

    Examples
    --------
    >>> normalize_column_name("Modified.Sequence")
    'modified_sequence'
    >>> normalize_column_name("Raw file")
    'raw_file'
    >>> normalize_column_name("Precursor.Quantity")
    'precursor_quantity'
    >>> normalize_column_name("Q.Value")
    'q_value'

    Parameters
    ----------
    col : str
        Original column name from software output.

    Returns
    -------
    str
        Normalized column name (lowercase, underscores, valid Python identifier).
    """
    # Convert to lowercase
    normalized = col.lower()
    # Replace whitespace, dots, hyphens with underscores
    normalized = re.sub(r"[\s\.\-]+", "_", normalized)
    # Remove other special characters
    normalized = re.sub(r"[^\w]", "", normalized)
    # Ensure doesn't start with a number
    if normalized and normalized[0].isdigit():
        normalized = "col_" + normalized
    # Handle empty result
    if not normalized:
        normalized = "column"
    return normalized


def normalize_all_columns(columns: list[str]) -> dict[str, str]:
    """
    Create a mapping from original to normalized column names.

    Parameters
    ----------
    columns : list of str
        List of original column names.

    Returns
    -------
    dict
        Mapping from original name to normalized name.
    """
    return {col: normalize_column_name(col) for col in columns}
