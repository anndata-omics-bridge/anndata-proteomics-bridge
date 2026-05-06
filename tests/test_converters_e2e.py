"""End-to-end: read → recognize → convert → AnnData for every packaged TOML.

Skips when the test_data_download cache (gitignored) is absent.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from anndata_proteomics.converters.assemble import convert
from anndata_proteomics.converters.recognize import recognize
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
    data_file = find_test_data(rule.software_name)
    if data_file is None or not data_file.exists():
        pytest.skip(f"no test data for {rule.software_name!r}")

    df = read_table(data_file)
    recognised = recognize(list(df.columns))
    assert recognised is not None
    assert recognised.software_name == rule.software_name

    adata = convert(df, recognised)
    assert adata.shape[0] > 0, f"{rule.software_name}: empty obs axis"
    assert adata.shape[1] > 0, f"{rule.software_name}: empty var axis"
    assert recognised.axis.x_layer in adata.layers
    assert adata.uns["anndata_proteomics"]["software_name"] == rule.software_name
