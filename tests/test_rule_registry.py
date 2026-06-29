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


def test_iter_packaged_rules_returns_eleven_sorted() -> None:
    rules = list(iter_packaged_rules())
    # 4 DIA-NN (ion + v1 fragment/protein + v2 protein) + 3 Spectronaut + maxquant/fragpipe/peaks/wombat
    assert len(rules) == 11
    assert rules == sorted(rules)  # path-sorted


def test_iter_packaged_rules_excludes_vendor_base_files() -> None:
    # diann/diann.toml and spectronaut/spectronaut.toml exist as inheritance bases but are not
    # rules — the parse_*.toml glob must skip them.
    root = packaged_rules_root()
    assert (root / "diann" / "diann.toml").exists()
    assert (root / "spectronaut" / "spectronaut.toml").exists()
    names = {p.name for p in iter_packaged_rules()}
    assert "diann.toml" not in names
    assert "spectronaut.toml" not in names
    assert all(p.name.startswith("parse_") for p in iter_packaged_rules())


def test_find_rule_happy() -> None:
    p = find_rule("diann", "ion")  # version-agnostic → vendor root
    assert p.name == "parse_diann_ion.toml"
    assert p.parent.name == "diann"
    assert p.exists()


def test_find_rule_version_folder() -> None:
    assert find_rule("diann", "protein", "1.9.2").parent.name == "v1"
    assert find_rule("diann", "protein", "2.3.0").parent.name == "v2"


def test_find_rule_unknown_software() -> None:
    with pytest.raises(RuleNotFound, match="nope"):
        find_rule("nope", "ion")


def test_find_rule_unknown_level() -> None:
    with pytest.raises(RuleNotFound, match="psm"):
        find_rule("diann", "psm")  # not a shipped level
