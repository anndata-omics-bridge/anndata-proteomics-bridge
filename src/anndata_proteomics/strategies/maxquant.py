"""MaxQuant evidence.txt strategy."""

from pathlib import Path

import pandas as pd

from ..proforma import ProFormaConverter
from ..utils import normalize_column_name


class MaxQuantStrategy:
    """
    Strategy for MaxQuant evidence.txt files.

    Handles MaxQuant-specific quirks:
    - Modified sequence format: _(Acetyl (Protein N-term))PEPTIDEK_
    - Fixed modifications from params
    - Multiple protein columns
    - Contaminant and reverse filtering
    """

    name = "MaxQuant"

    # Columns that indicate this is a MaxQuant file
    DETECTION_COLUMNS = ["Modified sequence", "Raw file", "Intensity", "Proteins"]

    # IDs for obs (samples) and var (features)
    obs_id = "Raw file"
    var_id = "precursor_id"  # Created from Modified sequence + Charge

    # Columns for var metadata (precursor-level information)
    VAR_COLUMNS = [
        "Modified sequence",
        "Sequence",
        "Charge",
        "Proteins",
        "Leading proteins",
        "Gene names",
        "Protein names",
        "Mass",
        "m/z",
        "Retention time",
        "Score",
    ]

    # Columns for layers (first = default X)
    LAYER_COLUMNS = [
        "Intensity",  # Default X
        "PEP",
        "Score",
        "Delta score",
        "MS/MS count",
    ]

    def __init__(self):
        # Load ProForma converter from TOML
        config_path = Path(__file__).parent.parent / "configs" / "maxquant.toml"
        if config_path.exists():
            self.proforma = ProFormaConverter.from_toml(config_path)
        else:
            self.proforma = None

    def detect(self, path: Path) -> bool:
        """Check if this is a MaxQuant evidence.txt file."""
        try:
            df = pd.read_csv(path, sep="\t", nrows=0)
            return all(col in df.columns for col in self.DETECTION_COLUMNS)
        except Exception:
            return False

    def load(self, path: Path) -> pd.DataFrame:
        """Load MaxQuant evidence.txt file."""
        df = pd.read_csv(path, sep="\t", low_memory=False)

        # Filter contaminants and reverse hits
        if "Reverse" in df.columns:
            df = df[df["Reverse"] != "+"]
        if "Potential contaminant" in df.columns:
            df = df[df["Potential contaminant"] != "+"]

        # Create precursor_id from Modified sequence and Charge
        df["precursor_id"] = (
            df["Modified sequence"].astype(str) + "/" + df["Charge"].astype(str)
        )

        return df

    def normalize_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Normalize all column names to valid Python identifiers.

        Preserves original column semantics, only changes formatting:
        - Modified sequence -> modified_sequence
        - Raw file -> raw_file
        """
        df.columns = [normalize_column_name(c) for c in df.columns]
        return df

    def get_obs(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Return obs metadata (unique per observation/sample).

        For MaxQuant, obs is primarily populated from the annotation file.
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
        if self.proforma and "Modified sequence" in df.columns:
            seq_values = df.groupby(self.var_id)["Modified sequence"].first()
            var_df["proforma"] = seq_values.apply(self.proforma.convert)
            var_df["stripped_sequence"] = seq_values.apply(self.proforma.get_stripped_sequence)

        return var_df

    def get_layers(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Return layers data with obs_id, var_id, and layer columns.

        The first layer column (Intensity) is the default X matrix.
        """
        # Start with obs_id and var_id
        cols = [self.obs_id, self.var_id]

        # Add available layer columns
        available_layers = [c for c in self.LAYER_COLUMNS if c in df.columns]
        cols.extend(available_layers)

        return df[cols].copy()

    def __repr__(self) -> str:
        return "MaxQuantStrategy()"
