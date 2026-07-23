"""Tests for software-version rule-document discovery."""

from __future__ import annotations

import pytest

from anndata_proteomics.rules.registry import (
    RuleNotFound,
    document_paths_for_software,
    find_rule,
    iter_packaged_documents,
    iter_packaged_rules,
    packaged_rules_root,
)


def test_packaged_rules_root_exists() -> None:
    assert packaged_rules_root().exists()


def test_iter_packaged_documents_returns_seven_sorted() -> None:
    documents = list(iter_packaged_documents())
    assert len(documents) == 7
    assert documents == sorted(documents)
    assert all(path.name == "rules.json" for path in documents)


def test_iter_packaged_rules_returns_twelve_document_levels() -> None:
    locators = list(iter_packaged_rules())
    assert len(locators) == 12
    assert sum(item.level == "ion" for item in locators) == 6


def test_diann_has_two_version_documents() -> None:
    paths = document_paths_for_software("diann")
    assert [path.parent.name for path in paths] == ["v1", "v2"]


def test_find_rule_resolves_existing_version_group() -> None:
    locator = find_rule("diann", "protein", "1.9.2")
    assert locator.path.parent.name == "v1"
    assert locator.level == "protein"
    assert find_rule("diann", "protein", "2.3.0").path.parent.name == "v2"


def test_find_rule_without_version_resolves_identical_diann_ion() -> None:
    locator = find_rule("diann", "ion")
    assert locator.path.parent.name == "v1"


def test_find_rule_unknown_software() -> None:
    with pytest.raises(RuleNotFound, match="nope"):
        find_rule("nope", "ion")
