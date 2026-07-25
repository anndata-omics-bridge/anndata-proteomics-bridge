"""Tests for readers/dispatch.py."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from anndata_proteomics.readers.dispatch import (
    UnknownFormat,
    read_table,
    read_table_columns,
)


def test_dispatch_csv(tmp_path: Path) -> None:
    p = tmp_path / "data.csv"
    p.write_text("a,b\n1,2\n")
    df = read_table(p)
    assert list(df.columns) == ["a", "b"]


def test_dispatch_tsv(tmp_path: Path) -> None:
    p = tmp_path / "data.tsv"
    p.write_text("a\tb\n1\t2\n")
    df = read_table(p)
    assert list(df.columns) == ["a", "b"]


def test_dispatch_tab_delimited_txt(tmp_path: Path) -> None:
    """.txt preserves the common MaxQuant-style tab-delimited convention."""
    p = tmp_path / "data.txt"
    p.write_text("a\tb\n1\t2\n")
    df = read_table(p)
    assert list(df.columns) == ["a", "b"]


def test_dispatch_comma_delimited_txt(tmp_path: Path) -> None:
    """.txt also supports comma-delimited vendor exports such as PEAKS."""
    p = tmp_path / "data.txt"
    p.write_text('"a","b"\n1,2\n')
    df = read_table(p)
    assert list(df.columns) == ["a", "b"]


def test_dispatch_one_column_txt_uses_safe_fallback(tmp_path: Path) -> None:
    p = tmp_path / "data.txt"
    p.write_text("only\nvalue\n")

    assert read_table(p).columns.tolist() == ["only"]


def test_dispatch_parquet(tmp_path: Path) -> None:
    p = tmp_path / "data.parquet"
    pd.DataFrame({"a": [1], "b": [2]}).to_parquet(p)
    df = read_table(p)
    assert list(df.columns) == ["a", "b"]


def test_dispatch_unknown_extension_raises(tmp_path: Path) -> None:
    p = tmp_path / "data.xyz"
    p.write_text("anything")
    with pytest.raises(UnknownFormat, match="xyz"):
        read_table(p)


def test_dispatch_extension_case_insensitive(tmp_path: Path) -> None:
    p = tmp_path / "DATA.CSV"
    p.write_text("a,b\n1,2\n")
    df = read_table(p)
    assert list(df.columns) == ["a", "b"]


def test_read_table_columns_uses_dispatch_without_loading_rows(tmp_path: Path) -> None:
    csv_path = tmp_path / "data.csv"
    csv_path.write_text("a,b\n1,2\n")
    tsv_path = tmp_path / "data.tsv"
    tsv_path.write_text("c\td\n3\t4\n")

    assert read_table_columns(csv_path) == ["a", "b"]
    assert read_table_columns(tsv_path) == ["c", "d"]

    unknown = tmp_path / "data.xyz"
    unknown.write_text("a,b\n")
    with pytest.raises(UnknownFormat, match="xyz"):
        read_table_columns(unknown)
