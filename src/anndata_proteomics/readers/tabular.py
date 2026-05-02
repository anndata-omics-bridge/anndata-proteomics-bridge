"""Generic file → pandas.DataFrame readers (no vendor semantics)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def read_csv(path: Path | str) -> pd.DataFrame:
    """Read a comma-delimited file. UTF-8 with BOM tolerance."""
    return pd.read_csv(path, encoding="utf-8-sig")


def read_tsv(path: Path | str) -> pd.DataFrame:
    """Read a tab-delimited file. UTF-8 with BOM tolerance."""
    return pd.read_csv(path, sep="\t", encoding="utf-8-sig")


def read_parquet(path: Path | str) -> pd.DataFrame:
    """Read a parquet file (via pyarrow)."""
    return pd.read_parquet(path)
