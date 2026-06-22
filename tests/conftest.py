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


# --- Shared DIA-NN report-backed multi-level test data -----------------------
# Conversion is param-driven: the param file gives the DIA-NN version, which selects the rule
# variants. A DIA-NN 1.9.x export supports ion, protein, and positional fragment levels, so the
# fixture finds a cached DIA-NN dataset whose version resolves those levels and returns a small row
# subset plus that version. Skips cleanly when no such dataset is cached.

_SUBSET_ROWS = 4000  # precursor rows; the fragment level explodes this ~12x


@pytest.fixture(scope="session")
def diann_full_subset() -> dict:
    """`{df, version, slug}` for a DIA-NN dataset with report-backed levels."""
    from anndata_proteomics.readers.dispatch import read_table
    from anndata_proteomics.scripts import _ui_support as ui

    catalog = ui.load_catalog()
    if catalog.empty:
        pytest.skip("no cached test-data catalog")
    required_levels = {"ion", "protein", "fragment"}
    diann = catalog[
        (catalog["slug"] == "diann")
        & catalog["targets"].apply(lambda targets: required_levels <= set(targets))
    ]
    if diann.empty:
        pytest.skip("no cached DIA-NN dataset with ion, protein, and fragment levels")
    row = diann.iloc[0]
    version = ui._param_version(Path(row["param_path"]), "diann")
    df = read_table(ui._dataset_path(row["input_file_path"]))
    run0 = df["Run"].iloc[0]
    return {
        "df": df[df["Run"] == run0].head(_SUBSET_ROWS).copy(),
        "version": version,
        "slug": "diann",
    }
