"""Tests for rules/registry.py."""

from __future__ import annotations

import pytest

from anndata_proteomics.rules.registry import (
    RuleNotFound,
    find_rule,
    iter_packaged_rules,
    packaged_rules_root,
)


def test_packaged_rules_root_exists() -> None:
    assert packaged_rules_root().exists()


def test_iter_packaged_rules_returns_ten_sorted() -> None:
    rules = list(iter_packaged_rules())
    assert len(rules) == 10
    vendors = [p.parent.name for p in rules]
    assert vendors == sorted(vendors)


def test_find_rule_happy() -> None:
    p = find_rule("diann", "ion")
    assert p.name == "parse_diann_ion_1.toml"
    assert p.parent.name == "diann"
    assert p.exists()


def test_find_rule_unknown_software() -> None:
    with pytest.raises(RuleNotFound, match="nope"):
        find_rule("nope", "ion")


def test_find_rule_unknown_level_lists_available() -> None:
    with pytest.raises(RuleNotFound) as exc_info:
        find_rule("diann", "psm")  # not a shipped level
    assert "parse_diann_ion_1.toml" in str(exc_info.value)
