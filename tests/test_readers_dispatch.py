"""Tests for readers/dispatch.py."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from anndata_proteomics.readers.dispatch import UnknownFormat, read_table


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


def test_dispatch_txt_treats_as_tsv(tmp_path: Path) -> None:
    """.txt is read as tab-delimited — MaxQuant convention."""
    p = tmp_path / "data.txt"
    p.write_text("a\tb\n1\t2\n")
    df = read_table(p)
    assert list(df.columns) == ["a", "b"]


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
