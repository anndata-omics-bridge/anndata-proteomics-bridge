"""Dispatch by file extension to the right tabular reader."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from anndata_proteomics.readers.tabular import (
    detect_text_delimiter,
    read_csv,
    read_csv_preserving_strings,
    read_delimited_columns,
    read_detected_text,
    read_detected_text_preserving_strings,
    read_parquet,
    read_tsv,
    read_tsv_preserving_strings,
)


class UnknownFormat(ValueError):
    """Raised when a file extension has no registered reader."""


EXTENSION_TO_READER = {
    ".csv": read_csv,
    ".tsv": read_tsv,
    ".txt": read_detected_text,
    ".parquet": read_parquet,
}

EXTENSION_TO_STRING_PRESERVING_READER = {
    ".csv": read_csv_preserving_strings,
    ".tsv": read_tsv_preserving_strings,
    ".txt": read_detected_text_preserving_strings,
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


def read_table_preserving_strings(
    path: Path | str,
    string_columns: frozenset[str],
) -> pd.DataFrame:
    """Read a table while preserving exact text for selected vendor sources.

    Delimited text uses the rule-derived source set to prevent pandas from treating
    identifiers as numbers. Parquet already carries a physical schema and is read as-is;
    conversion applies the logical output contract afterward.
    """
    source = Path(path)
    suffix = source.suffix.lower()
    if suffix == ".parquet":
        return read_parquet(source)
    reader = EXTENSION_TO_STRING_PRESERVING_READER.get(suffix)
    if reader is None:
        raise UnknownFormat(
            f"unsupported extension {source.suffix!r} for {source}; "
            f"known: {sorted(EXTENSION_TO_READER)}"
        )
    return reader(source, string_columns)


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
