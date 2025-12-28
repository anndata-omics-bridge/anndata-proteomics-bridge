"""DIA-NN report.tsv strategy."""

from pathlib import Path

import pandas as pd

from ..proforma import ProFormaConverter
from ..utils import normalize_column_name


class DIANNStrategy:
    """
    Strategy for DIA-NN report.tsv files.

    Handles:
    - TSV and Parquet formats
    - Modified.Sequence already in ProForma-like format
    - Multiple quantity columns (Precursor.Quantity, Precursor.Normalised, etc.)
    """

    name = "DIA-NN"

    # Columns that indicate this is a DIA-NN file
    DETECTION_COLUMNS = ["Modified.Sequence", "Run", "Precursor.Quantity", "Precursor.Id"]

    # IDs for obs (samples) and var (features)
    obs_id = "Run"
    var_id = "Precursor.Id"

    # Columns for var metadata (precursor-level information)
    VAR_COLUMNS = [
        "Modified.Sequence",
        "Stripped.Sequence",
        "Precursor.Charge",
        "Protein.Ids",
        "Protein.Names",
        "Genes",
        "Proteotypic",
        "Protein.Group",
        "PG.MaxLFQ",
    ]

    # Columns for layers (first = default X)
    LAYER_COLUMNS = [
        "Precursor.Quantity",  # Default X
        "Precursor.Normalised",
        "Ms1.Area",
        "Q.Value",
        "Global.Q.Value",
        "PG.Q.Value",
        "Lib.Q.Value",
        "PEP",
    ]

    def __init__(self):
        # Load ProForma converter from TOML
        config_path = Path(__file__).parent.parent / "configs" / "diann.toml"
        if config_path.exists():
            self.proforma = ProFormaConverter.from_toml(config_path)
        else:
            self.proforma = None

    def detect(self, path: Path) -> bool:
        """Check if this is a DIA-NN file."""
        path = Path(path)

        # Handle parquet files
        if path.suffix == ".parquet":
            try:
                df = pd.read_parquet(path, columns=self.DETECTION_COLUMNS[:1])
                return True
            except Exception:
                return False

        # Handle TSV files
        try:
            df = pd.read_csv(path, sep="\t", nrows=0)
            return all(col in df.columns for col in self.DETECTION_COLUMNS)
        except Exception:
            return False

    def load(self, path: Path) -> pd.DataFrame:
        """Load DIA-NN report file."""
        path = Path(path)

        if path.suffix == ".parquet":
            df = pd.read_parquet(path)
        else:
            df = pd.read_csv(path, sep="\t", low_memory=False)

        return df

    def normalize_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Normalize all column names to valid Python identifiers.

        Preserves original column semantics, only changes formatting:
        - Modified.Sequence -> modified_sequence
        - Precursor.Quantity -> precursor_quantity
        """
        df.columns = [normalize_column_name(c) for c in df.columns]
        return df

    def get_obs(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Return obs metadata (unique per observation/sample).

        For DIA-NN, obs is primarily populated from the annotation file.
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
        if self.proforma and "Modified.Sequence" in df.columns:
            seq_values = df.groupby(self.var_id)["Modified.Sequence"].first()
            var_df["proforma"] = seq_values.apply(self.proforma.convert)
            var_df["stripped_sequence"] = seq_values.apply(self.proforma.get_stripped_sequence)

        return var_df

    def get_layers(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Return layers data with obs_id, var_id, and layer columns.

        The first layer column (Precursor.Quantity) is the default X matrix.
        """
        # Start with obs_id and var_id
        cols = [self.obs_id, self.var_id]

        # Add available layer columns
        available_layers = [c for c in self.LAYER_COLUMNS if c in df.columns]
        cols.extend(available_layers)

        return df[cols].copy()

    def __repr__(self) -> str:
        return "DIANNStrategy()"
