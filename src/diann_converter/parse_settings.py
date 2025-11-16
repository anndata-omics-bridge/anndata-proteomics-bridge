"""
Configuration loader for software-specific column mappings.

Simplified from ProteoBench's parse_settings.py to handle TOML configs.
"""

import os
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd
import toml

from .proforma import get_proforma_bracketed


class SoftwareConfig:
    """
    Configuration for a specific proteomics software tool.

    Loads column mappings and modification parsing settings from TOML files.

    Parameters
    ----------
    config_path : str
        Path to the TOML configuration file.

    Attributes
    ----------
    mapper : dict
        Maps software-specific column names to standardized names.
    modifications_parser : dict or None
        Configuration for converting modifications to ProForma.
    general : dict
        General processing flags (contaminant_flag, decoy_flag, data_format).
    """

    def __init__(self, config_path: str):
        """Initialize configuration from TOML file."""
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Configuration file not found: {config_path}")

        config = toml.load(config_path)

        self.mapper = config.get("mapper", {})
        self.modifications_parser = config.get("modifications_parser", None)
        self.general = config.get("general", {})

        # Validate required sections
        if not self.mapper:
            raise ValueError(f"Configuration file must contain [mapper] section: {config_path}")

    def rename_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Rename columns according to mapper configuration.

        Only renames columns that exist in the dataframe.

        Parameters
        ----------
        df : pd.DataFrame
            Input dataframe.

        Returns
        -------
        pd.DataFrame
            Dataframe with renamed columns.
        """
        # Only rename columns that exist
        rename_dict = {k: v for k, v in self.mapper.items() if k in df.columns}
        return df.rename(columns=rename_dict)

    def convert_modifications(self, df: pd.DataFrame, sequence_col: str = "sequence") -> pd.DataFrame:
        """
        Convert modifications to ProForma notation.

        Adds a 'proforma' column to the dataframe.

        Parameters
        ----------
        df : pd.DataFrame
            Input dataframe.
        sequence_col : str, optional
            Name of the column containing modified sequences. Defaults to "sequence".

        Returns
        -------
        pd.DataFrame
            Dataframe with added 'proforma' column.
        """
        if self.modifications_parser is None:
            # No modification conversion needed - use sequence as-is
            df["proforma"] = df[sequence_col]
            return df

        # Extract parameters
        pattern = self.modifications_parser["pattern"]
        modification_dict = self.modifications_parser["modification_dict"]
        before_aa = self.modifications_parser.get("before_aa", False)
        isalpha = self.modifications_parser.get("isalpha", True)
        isupper = self.modifications_parser.get("isupper", True)

        # Apply conversion
        df["proforma"] = df[sequence_col].apply(
            lambda seq: get_proforma_bracketed(
                seq,
                before_aa=before_aa,
                isalpha=isalpha,
                isupper=isupper,
                pattern=pattern,
                modification_dict=modification_dict,
            )
        )

        return df

    def filter_decoys(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Filter out decoy entries if applicable.

        Parameters
        ----------
        df : pd.DataFrame
            Input dataframe.

        Returns
        -------
        pd.DataFrame
            Dataframe with decoys removed.
        """
        decoy_flag = self.general.get("decoy_flag", False)
        if not decoy_flag:
            return df

        # Check for common decoy column names
        decoy_cols = ["decoy", "Decoy", "DECOY", "Reverse", "is_decoy"]
        for col in decoy_cols:
            if col in df.columns:
                return df[df[col] == False]  # noqa: E712

        return df

    def mark_contaminants(self, df: pd.DataFrame, protein_col: str = "proteins") -> pd.DataFrame:
        """
        Mark contaminant proteins.

        Adds a 'contaminant' boolean column.

        Parameters
        ----------
        df : pd.DataFrame
            Input dataframe.
        protein_col : str, optional
            Name of the column containing protein identifiers. Defaults to "proteins".

        Returns
        -------
        pd.DataFrame
            Dataframe with added 'contaminant' column.
        """
        contaminant_flag = self.general.get("contaminant_flag", "Cont_")

        if protein_col not in df.columns:
            df["contaminant"] = False
            return df

        df["contaminant"] = df[protein_col].str.contains(contaminant_flag, na=False)
        return df


def load_config(software: str, config_dir: Optional[str] = None) -> SoftwareConfig:
    """
    Load configuration for a software tool.

    Parameters
    ----------
    software : str
        Software name (e.g., "diann", "spectronaut").
    config_dir : str, optional
        Directory containing TOML config files. If None, uses default config directory.

    Returns
    -------
    SoftwareConfig
        Configuration object for the software.

    Raises
    ------
    FileNotFoundError
        If configuration file doesn't exist.
    ValueError
        If software is not supported.

    Examples
    --------
    >>> config = load_config("diann")
    >>> df_renamed = config.rename_columns(df)
    """
    if config_dir is None:
        # Use default config directory relative to this file
        package_dir = Path(__file__).parent.parent.parent
        config_dir = package_dir / "config"

    config_path = Path(config_dir) / f"{software}.toml"

    if not config_path.exists():
        available = [f.stem for f in Path(config_dir).glob("*.toml")]
        raise ValueError(
            f"Unsupported software: {software}. "
            f"Available configurations: {', '.join(available)}"
        )

    return SoftwareConfig(str(config_path))
