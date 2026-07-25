"""Generic file → pandas.DataFrame readers (no vendor semantics)."""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

_TEXT_DELIMITERS = ",\t"


def detect_text_delimiter(path: Path | str) -> str:
    """Detect comma- versus tab-delimited text from a bounded content sample."""
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(65536)
    try:
        return csv.Sniffer().sniff(sample, delimiters=_TEXT_DELIMITERS).delimiter
    except csv.Error:
        # A one-column text file has no delimiter to detect. Preserve the historical
        # generic-.txt default so it still reads as one tabular column.
        return "\t"


def read_csv(path: Path | str) -> pd.DataFrame:
    """Read a comma-delimited file. UTF-8 with BOM tolerance."""
    return pd.read_csv(path, encoding="utf-8-sig")


def read_tsv(path: Path | str) -> pd.DataFrame:
    """Read a tab-delimited file. UTF-8 with BOM tolerance."""
    return pd.read_csv(path, sep="\t", encoding="utf-8-sig")


def read_detected_text(path: Path | str) -> pd.DataFrame:
    """Read comma- or tab-delimited text after content-based delimiter detection."""
    return pd.read_csv(
        path,
        sep=detect_text_delimiter(path),
        encoding="utf-8-sig",
    )


def read_delimited_columns(path: Path | str, *, delimiter: str) -> list[str]:
    """Read only a delimited text file's header."""
    return list(
        pd.read_csv(
            path,
            sep=delimiter,
            encoding="utf-8-sig",
            nrows=0,
        ).columns
    )


def read_parquet(path: Path | str) -> pd.DataFrame:
    """Read a parquet file (via pyarrow)."""
    return pd.read_parquet(path)
