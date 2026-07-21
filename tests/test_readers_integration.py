"""End-to-end: every packaged JSON rule reads its corresponding test-data file.

Skipped when the test_data_download cache (gitignored, regenerable) is absent.
"""

from __future__ import annotations

import pytest

from anndata_proteomics.readers.dispatch import read_table
from anndata_proteomics.rules.loader import load_rule
from anndata_proteomics.rules.registry import iter_packaged_rules
from anndata_proteomics.test_data import find_test_data


@pytest.mark.parametrize(
    "locator",
    list(iter_packaged_rules()),
    ids=lambda item: f"{item.path.parent.name}/{item.level}",
)
def test_reader_loads_test_data_for_packaged_rule(locator) -> None:
    rule = load_rule(locator)
    data_file = find_test_data(rule.software_name)
    if data_file is None or not data_file.exists():
        pytest.skip(
            f"no downloaded test data for {rule.software_name!r}; "
            f"regenerate via apb-testdata catalog/select/download"
        )
    df = read_table(data_file)
    expected_min_cols = (
        len(rule.columns.var.select) + len(rule.columns.obs.select) + len(rule.layers)
    )
    assert not df.empty, f"{data_file} produced an empty DataFrame"
    assert len(df.columns) >= expected_min_cols, (
        f"{data_file}: got {len(df.columns)} columns, rule expects at least {expected_min_cols}"
    )
