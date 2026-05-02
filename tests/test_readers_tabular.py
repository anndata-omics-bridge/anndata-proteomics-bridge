"""Tests for readers/tabular.py."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from anndata_proteomics.readers.tabular import read_csv, read_parquet, read_tsv


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
    p.write_bytes("﻿Sequence\tLength\nABC\t3\n".encode("utf-8"))
    df = read_tsv(p)
    assert list(df.columns) == ["Sequence", "Length"]
    assert df.iloc[0]["Sequence"] == "ABC"


def test_read_parquet_basic(tmp_path: Path) -> None:
    p = tmp_path / "data.parquet"
    pd.DataFrame({"a": [1, 2, 3], "b": ["x", "y", "z"]}).to_parquet(p)
    df = read_parquet(p)
    assert list(df.columns) == ["a", "b"]
    assert df.shape == (3, 2)
