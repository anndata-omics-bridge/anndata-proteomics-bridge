"""End-to-end: read → convert → AnnData for every packaged JSON rule.

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
from anndata_proteomics.converters.pipeline import software_slug
from anndata_proteomics.converters.recognize import matches
from anndata_proteomics.params.registry import parse_params
from anndata_proteomics.readers.dispatch import read_table
from anndata_proteomics.rules.loader import load_rule, load_rule_document
from anndata_proteomics.rules.registry import RuleLocator, iter_packaged_rules
from anndata_proteomics.test_data import find_test_data


def _level_available_for(locator: RuleLocator, data_file: Path) -> bool:
    """Whether ``locator.level`` is offered for the parameters beside ``data_file``.

    Ungated levels are always available. A gated level needs the cached parameter file that
    sits next to the vendor input; if there is none, the level cannot be confirmed.
    """
    document = load_rule_document(locator.path)
    if not document.levels[locator.level].requires_search_parameters:
        return True
    param_paths = sorted(data_file.parent.glob("param_0.*"))
    if not param_paths:
        return False
    slug = software_slug(document.software_name)
    return document.level_is_available(locator.level, parse_params(param_paths[0], slug))


@pytest.mark.parametrize(
    "locator",
    list(iter_packaged_rules()),
    ids=lambda item: f"{item.path.parent.name}/{item.level}",
)
def test_end_to_end_conversion(locator: RuleLocator) -> None:
    rule = load_rule(locator)
    if rule.fragments is not None:
        # The fragment level explodes the packed fragment lists ~12x; converting a full
        # report.tsv pivots millions of rows and peaks at many GB. Covered on a small
        # subset in test_diann_levels.py instead.
        pytest.skip("fragment level converted on a subset in test_diann_levels.py")

    data_file = find_test_data(rule.software_name, rule.software_version)
    if data_file is None or not data_file.exists():
        pytest.skip(f"no test data for {rule.software_name!r} {rule.software_version!r}")

    if not _level_available_for(locator, data_file):
        # A parameter-gated level (Sage's charge-collapsed vs charge-resolved lfq.tsv) is not
        # selectable from version or headers, so the first cached file for this vendor may be
        # the other level's. Production resolves this through the parsed parameters; so do we.
        pytest.skip(f"cached {rule.software_name} file is not {locator.level}-level")

    df = read_table(data_file)
    if not matches(list(df.columns), rule):
        # DIA-NN report schemas vary by version/config; the one cached file may not carry
        # every level's columns. That is "wrong variant for this level", not a failure.
        pytest.skip(f"cached {rule.software_name} file lacks columns for {locator.level}")
    adata = convert(df, rule)
    assert adata.shape[0] > 0, f"{rule.software_name}: empty obs axis"
    assert adata.shape[1] > 0, f"{rule.software_name}: empty var axis"
    assert rule.axis.x_layer in adata.layers
    assert adata.uns["anndata_proteomics"]["software_name"] == rule.software_name
