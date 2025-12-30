"""Spectronaut report strategy."""

from pathlib import Path

import pandas as pd

from ..proforma import ProFormaConverter
from ..utils import normalize_column_name


class SpectronautStrategy:
    """
    Strategy for Spectronaut export files.

    Handles:
    - TSV format with potential comma decimal separator
    - FG.LabeledSequence with underscore padding
    - Multiple quantity columns (FG.Quantity, EG.TotalQuantity)
    """

    name = "Spectronaut"
    level = "ion"

    # Columns that indicate this is a Spectronaut file
    DETECTION_COLUMNS = ["FG.LabeledSequence", "R.FileName", "FG.Quantity"]

    # IDs for obs (samples) and var (features)
    obs_id = "R.FileName"
    var_id = "precursor_id"  # Created from FG.LabeledSequence + FG.Charge

    # Columns for var metadata (precursor-level information)
    VAR_COLUMNS = [
        "FG.LabeledSequence",
        "EG.StrippedSequence",
        "FG.Charge",
        "PG.ProteinGroups",
        "PG.ProteinNames",
        "PG.Genes",
        "EG.IsDecoy",
        "EG.ModifiedSequence",
        "FG.PrecursorMz",
    ]

    # Columns for layers (first = default X)
    LAYER_COLUMNS = [
        "EG.TotalQuantity (Settings)",  # Default X
        "FG.Quantity",
        "EG.Qvalue",
        "FG.Qvalue",
        "PG.Qvalue",
        "EG.PEP",
    ]

    def __init__(self):
        # Load ProForma converter from TOML
        config_path = Path(__file__).parent.parent / "configs" / "spectronaut.toml"
        if config_path.exists():
            self.proforma = ProFormaConverter.from_toml(config_path)
        else:
            self.proforma = None

    def detect(self, path: Path) -> bool:
        """Check if this is a Spectronaut file."""
        path = Path(path)

        try:
            df = pd.read_csv(path, sep="\t", nrows=0)
            return all(col in df.columns for col in self.DETECTION_COLUMNS)
        except Exception:
            return False

    def load(self, path: Path) -> pd.DataFrame:
        """
        Load Spectronaut export file.

        Handles comma decimal separators and strips underscore padding.
        """
        path = Path(path)

        # Try loading with standard decimal separator
        df = pd.read_csv(path, sep="\t", low_memory=False)

        # If FG.Quantity is object type, might be using comma as decimal
        if "FG.Quantity" in df.columns and df["FG.Quantity"].dtype == object:
            df = pd.read_csv(path, sep="\t", low_memory=False, decimal=",")

        # Strip underscores from labeled sequence (Spectronaut adds _ padding)
        if "FG.LabeledSequence" in df.columns:
            df["FG.LabeledSequence"] = df["FG.LabeledSequence"].str.strip("_")

        # Create precursor_id from FG.LabeledSequence and FG.Charge
        if "FG.LabeledSequence" in df.columns and "FG.Charge" in df.columns:
            df["precursor_id"] = (
                df["FG.LabeledSequence"].astype(str) + "/" + df["FG.Charge"].astype(str)
            )

        return df

    def normalize_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Normalize all column names to valid Python identifiers.

        Preserves original column semantics, only changes formatting:
        - FG.LabeledSequence -> fg_labeledsequence
        - R.FileName -> r_filename
        """
        df.columns = [normalize_column_name(c) for c in df.columns]
        return df

    def get_obs(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Return obs metadata (unique per observation/sample).

        For Spectronaut, obs is primarily populated from the annotation file.
        This returns an empty DataFrame indexed by obs_id values.
        """
        unique_obs = df[self.obs_id].unique()
        return pd.DataFrame(index=unique_obs)

    def get_var(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Return var metadata (unique per variable/precursor).

        Extracts all VAR_COLUMNS that exist in the data, grouped by var_id.
        """
        # Find which VAR_COLUMNS exist in the dataframe
        available_cols = [c for c in self.VAR_COLUMNS if c in df.columns]

        # Group by var_id and take first value for each column
        var_df = df.groupby(self.var_id).first()[available_cols]

        # Add proforma column if converter is available
        if self.proforma and "FG.LabeledSequence" in df.columns:
            seq_values = df.groupby(self.var_id)["FG.LabeledSequence"].first()
            var_df["proforma"] = seq_values.apply(self.proforma.convert)
            var_df["stripped_sequence"] = seq_values.apply(self.proforma.get_stripped_sequence)

        return var_df

    def get_layers(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Return layers data with obs_id, var_id, and layer columns.

        The first layer column (EG.TotalQuantity) is the default X matrix.
        """
        # Start with obs_id and var_id
        cols = [self.obs_id, self.var_id]

        # Add available layer columns
        available_layers = [c for c in self.LAYER_COLUMNS if c in df.columns]
        cols.extend(available_layers)

        return df[cols].copy()

    def __repr__(self) -> str:
        return "SpectronautStrategy()"
