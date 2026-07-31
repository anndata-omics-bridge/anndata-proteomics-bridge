"""Generic file → pandas.DataFrame readers (no vendor semantics)."""

from __future__ import annotations

import csv
import re
from collections.abc import Hashable, Mapping
from itertools import islice
from pathlib import Path

import pandas as pd

_TEXT_DELIMITERS = ",\t"
_COMMA_DECIMAL_RE = re.compile(r"^-?\d+,(\d+)$")
_THOUSANDS_GROUP_WIDTH = 3
_DECIMAL_SAMPLE_LINES = 500


def detect_text_delimiter(path: Path | str) -> str:
    """Detect comma- versus tab-delimited text from a bounded content sample."""
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(65536)
    if not any(delimiter in sample for delimiter in _TEXT_DELIMITERS):
        # A one-column text file has no delimiter to detect. Preserve the historical
        # generic-.txt default so it still reads as one tabular column.
        return "\t"
    return csv.Sniffer().sniff(sample, delimiters=_TEXT_DELIMITERS).delimiter


def detect_decimal_separator(path: Path | str, *, delimiter: str) -> str:
    """Detect whether a delimited text file writes numbers with a comma decimal mark.

    Vendors export numbers in the regional format of the machine that produced the file,
    and nothing in the file declares which one was used, so this is inferred from content.
    Only the shape of the number distinguishes the two readings of ``1,234``: a thousands
    separator always groups exactly three digits, while a decimal comma is followed by a
    fraction of any other width. A field whose comma is followed by three digits is
    therefore left ambiguous and never counted as evidence.

    A comma-delimited file cannot carry bare comma decimals at all, so it is reported as
    dot-decimal without inspection.
    """
    if delimiter == ",":
        return "."
    decimal_like = 0
    # Tolerant decoding: this scan only looks for digits and commas, and must not turn a
    # file pandas can still read into a decode failure before pandas ever sees it.
    with Path(path).open(encoding="utf-8-sig", newline="", errors="replace") as handle:
        handle.readline()
        for line in islice(handle, _DECIMAL_SAMPLE_LINES):
            for field in line.rstrip("\n").split(delimiter):
                match = _COMMA_DECIMAL_RE.match(field)
                if match is not None and len(match.group(1)) != _THOUSANDS_GROUP_WIDTH:
                    decimal_like += 1
    return "," if decimal_like else "."


def _read_delimited(path: Path | str, delimiter: str) -> pd.DataFrame:
    """Read delimited text in the number format the file itself was written with."""
    return pd.read_csv(
        path,
        sep=delimiter,
        encoding="utf-8-sig",
        decimal=detect_decimal_separator(path, delimiter=delimiter),
    )


def _read_delimited_preserving_strings(
    path: Path | str,
    delimiter: str,
    string_columns: frozenset[str],
) -> pd.DataFrame:
    """Read delimited data while preserving rule-declared textual source tokens."""
    dtypes: Mapping[Hashable, str] = dict.fromkeys(string_columns, "string")
    return pd.read_csv(
        path,
        sep=delimiter,
        encoding="utf-8-sig",
        decimal=detect_decimal_separator(path, delimiter=delimiter),
        dtype=dtypes,
    )


def read_csv(path: Path | str) -> pd.DataFrame:
    """Read a comma-delimited file. UTF-8 with BOM tolerance."""
    return _read_delimited(path, ",")


def read_csv_preserving_strings(
    path: Path | str,
    string_columns: frozenset[str],
) -> pd.DataFrame:
    """Read comma-delimited data with declared textual sources preserved."""
    return _read_delimited_preserving_strings(path, ",", string_columns)


def read_tsv(path: Path | str) -> pd.DataFrame:
    """Read a tab-delimited file. UTF-8 with BOM tolerance."""
    return _read_delimited(path, "\t")


def read_tsv_preserving_strings(
    path: Path | str,
    string_columns: frozenset[str],
) -> pd.DataFrame:
    """Read tab-delimited data with declared textual sources preserved."""
    return _read_delimited_preserving_strings(path, "\t", string_columns)


def read_detected_text(path: Path | str) -> pd.DataFrame:
    """Read comma- or tab-delimited text after content-based delimiter detection."""
    return _read_delimited(path, detect_text_delimiter(path))


def read_detected_text_preserving_strings(
    path: Path | str,
    string_columns: frozenset[str],
) -> pd.DataFrame:
    """Read detected delimited text with declared textual sources preserved."""
    delimiter = detect_text_delimiter(path)
    return _read_delimited_preserving_strings(path, delimiter, string_columns)


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
