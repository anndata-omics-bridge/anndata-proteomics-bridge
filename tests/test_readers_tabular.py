"""Tests for readers/tabular.py."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from anndata_proteomics.readers.tabular import (
    detect_decimal_separator,
    read_csv,
    read_parquet,
    read_tsv,
)


def test_read_csv_basic(tmp_path: Path) -> None:
    p = tmp_path / "data.csv"
    p.write_text("a,b,c\n1,2,3\n4,5,6\n")
    df = read_csv(p)
    assert list(df.columns) == ["a", "b", "c"]
    assert df.shape == (2, 3)


def test_read_csv_quoted_embedded_commas(tmp_path: Path) -> None:
    p = tmp_path / "quoted.csv"
    p.write_text('name,description\n"Alpha","one,two,three"\n"Beta","x,y"\n')
    df = read_csv(p)
    assert df.iloc[0]["description"] == "one,two,three"
    assert df.shape == (2, 2)


def test_read_tsv_basic(tmp_path: Path) -> None:
    p = tmp_path / "data.tsv"
    p.write_text("a\tb\tc\n1\t2\t3\n4\t5\t6\n")
    df = read_tsv(p)
    assert list(df.columns) == ["a", "b", "c"]
    assert df.shape == (2, 3)


def test_read_tsv_strips_utf8_bom(tmp_path: Path) -> None:
    """Some MaxQuant-adjacent exports start with a UTF-8 BOM; utf-8-sig drops it."""
    p = tmp_path / "bom.tsv"
    p.write_bytes("﻿Sequence\tLength\nABC\t3\n".encode())
    df = read_tsv(p)
    assert list(df.columns) == ["Sequence", "Length"]
    assert df.iloc[0]["Sequence"] == "ABC"


def test_read_parquet_basic(tmp_path: Path) -> None:
    p = tmp_path / "data.parquet"
    pd.DataFrame({"a": [1, 2, 3], "b": ["x", "y", "z"]}).to_parquet(p)
    df = read_parquet(p)
    assert list(df.columns) == ["a", "b"]
    assert df.shape == (3, 2)


def test_detects_a_comma_decimal_export_from_fraction_width(tmp_path: Path) -> None:
    # Spectronaut writes numbers in the producing machine's locale; nothing in the file
    # declares it, so the fraction width is the only signal.
    path = tmp_path / "comma.tsv"
    path.write_text(
        "Run\tQuantity\nr1\t702,1904907226562\nr2\t168,5559844970703\n",
        encoding="utf-8",
    )

    assert detect_decimal_separator(path, delimiter="\t") == ","
    assert read_tsv(path)["Quantity"].dtype.kind == "f"


def test_dot_decimal_export_is_left_alone(tmp_path: Path) -> None:
    path = tmp_path / "dot.tsv"
    path.write_text("Run\tQuantity\nr1\t702.19\nr2\t168.55\n", encoding="utf-8")

    assert detect_decimal_separator(path, delimiter="\t") == "."
    assert read_tsv(path)["Quantity"].dtype.kind == "f"


def test_three_digit_groups_stay_ambiguous_so_thousands_do_not_trip_detection(
    tmp_path: Path,
) -> None:
    # 1,234 reads as either 1.234 or 1234; a thousands separator always groups exactly
    # three digits, so that width is never taken as evidence of a decimal comma.
    path = tmp_path / "thousands.tsv"
    path.write_text("Run\tCount\nr1\t1,234\nr2\t9,876\n", encoding="utf-8")

    assert detect_decimal_separator(path, delimiter="\t") == "."


def test_comma_delimited_files_are_reported_without_inspection(tmp_path: Path) -> None:
    path = tmp_path / "csv.csv"
    path.write_text("Run,Quantity\nr1,702.19\n", encoding="utf-8")

    assert detect_decimal_separator(path, delimiter=",") == "."


def test_detection_tolerates_bytes_that_are_not_utf8(tmp_path: Path) -> None:
    # One corpus input is not valid UTF-8; detection must not turn a file pandas can
    # still read into a decode failure before pandas sees it.
    path = tmp_path / "latin.tsv"
    path.write_bytes(b"Run\tName\tQuantity\nr1\tcaf\xe9\t702,19\n")

    assert detect_decimal_separator(path, delimiter="\t") == ","
