"""Dispatch by file extension to the right tabular reader."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from anndata_proteomics.readers.tabular import (
    detect_text_delimiter,
    read_csv,
    read_delimited_columns,
    read_detected_text,
    read_parquet,
    read_tsv,
)


class UnknownFormat(ValueError):
    """Raised when a file extension has no registered reader."""


EXTENSION_TO_READER = {
    ".csv": read_csv,
    ".tsv": read_tsv,
    ".txt": read_detected_text,
    ".parquet": read_parquet,
}


def read_table(path: Path | str) -> pd.DataFrame:
    """Read a tabular file, dispatching by extension.

    Raises UnknownFormat if the extension is not registered.
    """
    p = Path(path)
    reader = EXTENSION_TO_READER.get(p.suffix.lower())
    if reader is None:
        raise UnknownFormat(
            f"unsupported extension {p.suffix!r} for {p}; known: {sorted(EXTENSION_TO_READER)}"
        )
    return reader(p)


def read_table_columns(path: Path | str) -> list[str]:
    """Read only column names using the same format rules as :func:`read_table`."""
    p = Path(path)
    suffix = p.suffix.lower()
    if suffix == ".parquet":
        return list(pq.read_schema(p).names)
    if suffix == ".csv":
        delimiter = ","
    elif suffix == ".tsv":
        delimiter = "\t"
    elif suffix == ".txt":
        delimiter = detect_text_delimiter(p)
    else:
        raise UnknownFormat(
            f"unsupported extension {p.suffix!r} for {p}; known: {sorted(EXTENSION_TO_READER)}"
        )
    return read_delimited_columns(p, delimiter=delimiter)
