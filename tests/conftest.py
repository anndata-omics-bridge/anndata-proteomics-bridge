"""Bridge loguru output into pytest's capsys/capfd.

Loguru's default sink captures `sys.stderr` at handler-registration time, which
runs once at module import. Pytest's `capsys` monkeypatches `sys.stderr` per
test, so without a bridge, loguru output bypasses capture and tests can't
assert on it. We replace the default sink with one whose writer callable
looks up `sys.stderr` at *write* time, picking up whatever pytest has patched
in for the current test.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest
from loguru import logger


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


# --- Shared cached ProteoBench test data (catalog-free) ----------------------
# The ProteoBench browser/catalog now lives in apb_studio; apb tests read the cached index
# (test_data.DOWNLOADED_DB) directly so they stay self-contained. Conversion is param-driven: the
# co-located param file gives the software version, which selects the rule variants.

_SUBSET_ROWS = 4000  # precursor rows; the fragment level explodes this ~12x


def _read_headers(path: Path) -> set[str]:
    """Column names of a cached input (cheap; tsv via pandas, parquet via the arrow schema)."""
    import pandas as pd
    import pyarrow.parquet as pq

    if path.suffix == ".parquet":
        return set(pq.read_schema(path).names)
    return set(pd.read_csv(path, sep="\t", nrows=0).columns)


@pytest.fixture(scope="session")
def cached_datasets():
    """Callable ``slug -> list[dict]`` over cached ProteoBench inputs for a vendor.

    Each dict has ``input_path``, ``param_path`` (or None), ``version``, ``headers``, ``targets``.
    Reads ``test_data.DOWNLOADED_DB`` directly (no apb_studio catalog) and derives targets from the
    packaged rules via ``converters.pipeline``. Empty list when the gitignored cache is absent.
    """
    from anndata_proteomics.converters import pipeline
    from anndata_proteomics.test_data import DOWNLOADED_DB, TEST_DATA_DIR

    def lookup(slug: str) -> list[dict]:
        if not DOWNLOADED_DB.exists():
            return []
        rows: list[dict] = []
        seen: set[str] = set()
        with open(DOWNLOADED_DB) as f:
            for row in csv.DictReader(f):
                if row.get("status") != "ok":
                    continue
                if pipeline.software_slug(row["software_name"]) != slug:
                    continue
                rel = row["input_file_path"]
                if rel in seen:
                    continue
                seen.add(rel)
                input_path = TEST_DATA_DIR / "json_dir" / rel
                params = sorted(input_path.parent.glob("param_0.*"))
                param_path = params[0] if params else None
                try:
                    headers = _read_headers(input_path)
                except OSError:
                    continue
                version = pipeline._param_version(param_path, slug) if param_path else None
                rows.append(
                    {
                        "input_file_path": rel,
                        "input_path": input_path,
                        "param_path": param_path,
                        "version": version,
                        "headers": headers,
                        "targets": pipeline.available_targets(slug, version, headers),
                    }
                )
        return rows

    return lookup


@pytest.fixture(scope="session")
def spectronaut_datasets(cached_datasets) -> list[dict]:
    """Cached Spectronaut datasets (see ``cached_datasets``)."""
    return cached_datasets("spectronaut")


@pytest.fixture(scope="session")
def diann_full_subset(cached_datasets) -> dict:
    """`{df, version, slug}` for a DIA-NN dataset with report-backed ion/protein/fragment levels."""
    from anndata_proteomics.readers.dispatch import read_table

    required_levels = {"ion", "protein", "fragment"}
    matches = [d for d in cached_datasets("diann") if required_levels <= set(d["targets"])]
    if not matches:
        pytest.skip("no cached DIA-NN dataset with ion, protein, and fragment levels")
    dataset = matches[0]
    df = read_table(dataset["input_path"])
    run0 = df["Run"].iloc[0]
    return {
        "df": df[df["Run"] == run0].head(_SUBSET_ROWS).copy(),
        "version": dataset["version"],
        "slug": "diann",
    }
