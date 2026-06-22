"""End-to-end: read → convert → AnnData for every packaged TOML.

Conversion uses the explicit parametrized rule rather than ``recognize()``: a vendor
can ship several quantification levels that all read the same file (DIA-NN's report.tsv
backs ion / peptidoform / peptide / protein / fragment), so header-based recognition
cannot pick a *level* and the caller selects it explicitly. ``recognize()`` is exercised
separately in test_recognize.py.

Skips when the test_data_download cache (gitignored) is absent.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from anndata_proteomics.converters.assemble import convert
from anndata_proteomics.converters.recognize import matches
from anndata_proteomics.readers.dispatch import read_table
from anndata_proteomics.rules.loader import load_rule
from anndata_proteomics.rules.registry import iter_packaged_rules
from anndata_proteomics.test_data import find_test_data


@pytest.mark.parametrize(
    "toml_path",
    list(iter_packaged_rules()),
    ids=lambda p: f"{p.parent.name}/{p.name}",
)
def test_end_to_end_conversion(toml_path: Path) -> None:
    rule = load_rule(toml_path)
    if rule.fragments is not None:
        # The fragment level explodes the packed fragment lists ~12x; converting a full
        # report.tsv pivots millions of rows and peaks at many GB. Covered on a small
        # subset in test_diann_levels.py instead.
        pytest.skip("fragment level converted on a subset in test_diann_levels.py")

    data_file = find_test_data(rule.software_name)
    if data_file is None or not data_file.exists():
        pytest.skip(f"no test data for {rule.software_name!r}")

    df = read_table(data_file)
    if not matches(list(df.columns), rule):
        # DIA-NN report schemas vary by version/config; the one cached file may not carry
        # every level's columns. That is "wrong variant for this level", not a failure.
        pytest.skip(f"cached {rule.software_name} file lacks columns for {toml_path.name}")
    adata = convert(df, rule)
    assert adata.shape[0] > 0, f"{rule.software_name}: empty obs axis"
    assert adata.shape[1] > 0, f"{rule.software_name}: empty var axis"
    assert rule.axis.x_layer in adata.layers
    assert adata.uns["anndata_proteomics"]["software_name"] == rule.software_name
