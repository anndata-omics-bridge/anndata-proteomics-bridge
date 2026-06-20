"""Consistency between the param-parser registry, the packaged conversion
rules, and the ProteoBench param fixtures.

These three surfaces have **different coverage by design**: the param-parser
registry (``params.registry._REGISTRY``) is the superset — it includes vendors
ported from ProteoBench that are not yet packaged with conversion rules (e.g.
alphapept / metamorpheus / msaid / sage). The ``parsing_rules/`` directory and
``test_data._PROTEOBENCH_PARAM_FIXTURES`` cover the same smaller set of
*packaged* tools.

So we do **not** assert all three are equal. Instead we assert the invariants
that must hold, resolving each surface through ``get_parser`` so that naming /
casing differences (``"DIA-NN"`` vs the ``diann`` rule dir vs the ``dia-nn``
registry alias) don't cause false failures:

* every param fixture maps to a registered parser (else a runtime ``KeyError``),
* every packaged rule directory maps to a registered parser,
* the packaged rules and the param fixtures describe the *same* set of tools.

This catches the drift the June code review flagged (a rule or fixture for a
vendor with no parser; a packaged rule with no fixture or vice versa).
"""

from __future__ import annotations

from pathlib import Path

from anndata_proteomics.params import registry
from anndata_proteomics.test_data import _PROTEOBENCH_PARAM_FIXTURES

_PARSING_RULES_DIR = Path(registry.__file__).resolve().parent.parent / "parsing_rules"


def _rule_vendor_dirs() -> set[str]:
    """Vendor slugs that ship a packaged conversion-rule directory."""
    return {
        p.name
        for p in _PARSING_RULES_DIR.iterdir()
        if p.is_dir() and not p.name.startswith(("_", "."))
    }


def test_every_param_fixture_has_a_registered_parser() -> None:
    for software_name in _PROTEOBENCH_PARAM_FIXTURES:
        registry.get_parser(software_name)  # raises KeyError if unregistered


def test_every_packaged_rule_has_a_registered_parser() -> None:
    for vendor in _rule_vendor_dirs():
        registry.get_parser(vendor)  # raises KeyError if unregistered


def test_packaged_rules_and_param_fixtures_cover_the_same_tools() -> None:
    """Packaged conversion rules and param fixtures must agree on the tool set.

    Resolved through the registry so that ``"DIA-NN"`` (fixture key) and
    ``"diann"`` (rule dir) collapse to the same parser.
    """
    fixture_parsers = {
        registry.get_parser(name) for name in _PROTEOBENCH_PARAM_FIXTURES
    }
    rule_parsers = {registry.get_parser(vendor) for vendor in _rule_vendor_dirs()}
    assert fixture_parsers == rule_parsers
