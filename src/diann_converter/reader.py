"""
File readers for different proteomics software outputs.

Adapted from ProteoBench (proteobench/io/parsing/parse_ion.py)
https://github.com/Proteobench/ProteoBench
"""

import pandas as pd


def load_diann(input_path: str) -> pd.DataFrame:
    """
    Load a DIA-NN output file.

    Automatically detects file format (TSV or Parquet) based on file extension.

    Parameters
    ----------
    input_path : str
        The path to the DIA-NN output file (.tsv or .parquet).

    Returns
    -------
    pd.DataFrame
        The loaded dataframe.

    Examples
    --------
    >>> df = load_diann("report.tsv")
    >>> df = load_diann("report.parquet")
    """
    if input_path.endswith(".parquet"):
        return pd.read_parquet(input_path)
    else:
        return pd.read_csv(input_path, low_memory=False, sep="\t")


def load_spectronaut(input_path: str) -> pd.DataFrame:
    """
    Load a Spectronaut output file.

    Handles both comma and period decimal separators.
    Strips underscores from labeled sequences.

    Parameters
    ----------
    input_path : str
        The path to the Spectronaut output file (.tsv).

    Returns
    -------
    pd.DataFrame
        The loaded dataframe.

    Examples
    --------
    >>> df = load_spectronaut("spectronaut_report.tsv")
    """
    # Try loading with standard decimal separator
    df = pd.read_csv(input_path, low_memory=False, sep="\t")

    # If FG.Quantity is object type, might be using comma as decimal
    if "FG.Quantity" in df.columns and df["FG.Quantity"].dtype == object:
        df = pd.read_csv(input_path, low_memory=False, sep="\t", decimal=",")

    # Strip underscores from labeled sequence (Spectronaut adds _ padding)
    if "FG.LabeledSequence" in df.columns:
        df["FG.LabeledSequence"] = df["FG.LabeledSequence"].str.strip("_")

    return df


def auto_detect_software(input_path: str) -> str:
    """
    Auto-detect software type from file structure.

    Checks for software-specific column names.

    Parameters
    ----------
    input_path : str
        The path to the input file.

    Returns
    -------
    str
        Software name: "diann", "spectronaut", or "unknown".

    Examples
    --------
    >>> software = auto_detect_software("report.tsv")
    >>> print(software)
    'diann'
    """
    # Read first few rows to check columns
    if input_path.endswith(".parquet"):
        df_sample = pd.read_parquet(input_path, nrows=10)
    else:
        df_sample = pd.read_csv(input_path, sep="\t", nrows=10)

    columns = set(df_sample.columns)

    # DIA-NN specific columns
    if "Precursor.Normalised" in columns or "Modified.Sequence" in columns:
        return "diann"

    # Spectronaut specific columns
    if "FG.LabeledSequence" in columns or "PG.ProteinGroups" in columns:
        return "spectronaut"

    return "unknown"


def load_file(input_path: str, software: str = None) -> pd.DataFrame:
    """
    Load proteomics file with optional auto-detection of software.

    Parameters
    ----------
    input_path : str
        The path to the input file.
    software : str, optional
        Software name ("diann" or "spectronaut"). If None, auto-detects.

    Returns
    -------
    pd.DataFrame
        The loaded dataframe.

    Raises
    ------
    ValueError
        If software type is unknown or not supported.

    Examples
    --------
    >>> df = load_file("report.tsv", software="diann")
    >>> df = load_file("report.tsv")  # Auto-detect
    """
    if software is None:
        software = auto_detect_software(input_path)

    loaders = {
        "diann": load_diann,
        "spectronaut": load_spectronaut,
    }

    if software not in loaders:
        raise ValueError(
            f"Unsupported or unknown software: {software}. "
            f"Supported software: {', '.join(loaders.keys())}"
        )

    return loaders[software](input_path)
