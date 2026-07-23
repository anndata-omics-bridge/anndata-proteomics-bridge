"""Verify packaged software-version documents and every effective level."""

from __future__ import annotations

import pytest

from anndata_proteomics.rules.loader import load_packaged_rule, load_rule, load_rule_document
from anndata_proteomics.rules.registry import (
    RuleLocator,
    iter_packaged_documents,
    iter_packaged_rules,
)
from anndata_proteomics.rules.validate import validate_all_packaged


def _locator_id(locator: RuleLocator) -> str:
    return f"{locator.path.parent.name}/{locator.level}"


def test_all_packaged_documents_validate() -> None:
    results = validate_all_packaged()
    failed = [result for result in results if not result.ok]
    assert not failed, "\n".join(f"{result.path}: {result.error}" for result in failed)


def test_at_least_one_long_and_one_wide_rule() -> None:
    shapes = {load_rule(locator).input_shape for locator in iter_packaged_rules()}
    assert shapes == {"long", "wide"}


def test_all_documents_declare_software_version() -> None:
    for path in iter_packaged_documents():
        assert load_rule_document(path).software_version


@pytest.mark.parametrize(
    "locator",
    list(iter_packaged_rules()),
    ids=_locator_id,
)
def test_locator_level_matches_effective_rule(locator: RuleLocator) -> None:
    assert load_rule(locator).quantification_level == locator.level


def test_documents_use_uniform_base_and_levels_shape() -> None:
    for path in iter_packaged_documents():
        source = path.read_text(encoding="utf-8")
        assert '"base"' in source
        assert '"levels"' in source
        assert '"$extends"' not in source


def test_fragpipe_intensity_declares_zero_as_missing() -> None:
    rule = load_packaged_rule("fragpipe", "ion", "22.1-build02")
    by_name = {layer.name: layer for layer in rule.layers}
    assert by_name["Intensity"].missing_values == [0.0]
    assert by_name["Spectral_Count"].missing_values == []
    assert by_name["Match_Type"].missing_values == []
