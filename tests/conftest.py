"""Bridge loguru output into pytest's capsys/capfd.

Loguru's default sink captures `sys.stderr` at handler-registration time, which
runs once at module import. Pytest's `capsys` monkeypatches `sys.stderr` per
test, so without a bridge, loguru output bypasses capture and tests can't
assert on it. We replace the default sink with one whose writer callable
looks up `sys.stderr` at *write* time, picking up whatever pytest has patched
in for the current test.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
import pytest
from loguru import logger

from anndata_proteomics.converters.recognize import matches
from anndata_proteomics.rules.loader import load_rule
from anndata_proteomics.rules.registry import iter_packaged_rules


@pytest.fixture(autouse=True)
def _loguru_to_pytest_capsys():
    logger.remove()
    logger.add(
        lambda msg: sys.stderr.write(msg),
        format="{level: <7} | {message}",
        level="DEBUG",
    )
    yield
    # Drop every sink — including any added by code under test (e.g.
    # configure_default_sink() or per-run file sinks). Avoids leaking sinks
    # across tests and tolerates main()-style code that calls logger.remove().
    logger.remove()


# --- Shared DIA-NN multi-level test data -------------------------------------
# DIA-NN's report.tsv backs every quantification level, but its fragment/PG columns vary
# a lot across versions (some exports drop Fragment.Info or the PG.* quant columns). The
# helpers below glob the on-disk cache (gitignored; the curated index and the disk diverge)
# for the first file carrying every level's columns, so the multi-level/MuData tests can
# build all five levels from one file. They skip cleanly when no such file is cached.

_CACHE_DIR = Path(__file__).resolve().parent.parent / "test_data_download" / "json_dir"
_SUBSET_ROWS = 4000  # precursor rows; the fragment level explodes this ~12x


def _diann_rules() -> list:
    return [load_rule(p) for p in iter_packaged_rules() if p.parent.name == "diann"]


def _headers(path: Path) -> set[str]:
    if path.suffix == ".parquet":
        return set(pq.read_schema(path).names)
    return set(pd.read_csv(path, sep="\t", nrows=0).columns)


def _carries_all_levels(headers: set[str]) -> bool:
    rules = _diann_rules()
    if not all(matches(headers, r) for r in rules):
        return False
    frag = next(r for r in rules if r.fragments is not None)
    return frag.fragments.label_column in headers  # not covered by matches()


def find_full_diann_file() -> Path | None:
    """First cached DIA-NN file whose columns cover all five quantification levels."""
    if not _CACHE_DIR.exists():
        return None
    candidates = sorted(_CACHE_DIR.glob("**/input_file.tsv")) + sorted(
        _CACHE_DIR.glob("**/input_file.txt")
    )
    for path in candidates:
        try:
            if _carries_all_levels(_headers(path)):
                return path
        except (OSError, ValueError):
            continue
    return None


@pytest.fixture(scope="session")
def diann_full_subset() -> pd.DataFrame:
    """One run, capped at a few thousand precursor rows, from a full-column DIA-NN file."""
    path = find_full_diann_file()
    if path is None:
        pytest.skip("no cached DIA-NN file carrying all five level columns")
    from anndata_proteomics.readers.dispatch import read_table

    df = read_table(path)
    run0 = df["Run"].iloc[0]
    return df[df["Run"] == run0].head(_SUBSET_ROWS).copy()
