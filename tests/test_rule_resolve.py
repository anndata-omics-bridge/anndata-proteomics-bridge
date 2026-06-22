"""Version-folder rule resolution + column validation.

DIA-NN report columns vary by version, so version-dependent levels live in version subfolders
(``diann/v1/``, ``diann/v2/``) selected by the software version parsed from the param file;
version-agnostic levels stay at the vendor root. Synthetic header sets (derived from the shipped
rules) keep these tests data-free.
"""

from __future__ import annotations

from anndata_proteomics.converters.recognize import _expected_long_columns
from anndata_proteomics.rules.loader import resolve_rule_for_version
from anndata_proteomics.rules.registry import resolve_rule_path
from anndata_proteomics.scripts import _ui_support as ui

_V19 = "1.9.2"
_V23 = "2.3.0 Academia "  # messy real catalog string


def _headers_for(rule) -> set[str]:
    cols = set(_expected_long_columns(rule))
    if rule.fragments is not None and rule.fragments.label_column:
        cols.add(rule.fragments.label_column)
    return cols


def _diann_headers(version: str) -> set[str]:
    cols: set[str] = set()
    for level in ui.LEVELS:
        rule = resolve_rule_for_version("diann", level, version)
        if rule is not None:
            cols |= _headers_for(rule)
    return cols


def test_resolve_path_picks_version_folder() -> None:
    assert resolve_rule_path("diann", "protein", _V19).parent.name == "v1"
    assert resolve_rule_path("diann", "protein", _V23).parent.name == "v2"
    assert resolve_rule_path("diann", "fragment", _V19).parent.name == "v1"
    assert resolve_rule_path("diann", "fragment", _V23) is None  # no fragment for 2.x
    # version-agnostic levels live at the vendor root
    assert resolve_rule_path("diann", "ion", _V19).parent.name == "diann"
    assert resolve_rule_path("diann", "ion", _V23).parent.name == "diann"


def test_resolve_flat_vendor_and_unknown() -> None:
    assert resolve_rule_path("maxquant", "ion", None).name == "parse_maxquant_ion_1.toml"
    assert resolve_rule_path("nope", "ion", "1.0") is None


def test_protein_variants_differ_by_version() -> None:
    v1_layers = {layer.name for layer in resolve_rule_for_version("diann", "protein", _V19).layers}
    v2_layers = {layer.name for layer in resolve_rule_for_version("diann", "protein", _V23).layers}
    assert "PG_Normalised" in v1_layers
    assert "PG_Normalised" not in v2_layers  # dropped in DIA-NN 2.x


def test_fragment_v1_is_positional() -> None:
    frag = resolve_rule_for_version("diann", "fragment", _V19)
    assert frag.fragments is not None
    assert frag.fragments.label_column is None  # positional labels (no Fragment.Info)


def test_convertible_levels_by_version() -> None:
    assert ui.convertible_levels("diann", _V19, _diann_headers(_V19)) == ui.LEVELS  # all 5
    assert ui.convertible_levels("diann", _V23, _diann_headers(_V23)) == [
        "ion", "peptidoform", "peptide", "protein",
    ]
    assert "mudata" in ui.available_targets("diann", _V23, _diann_headers(_V23))


def test_select_rule_errors() -> None:
    headers = _diann_headers(_V23)
    # fragment has no rule covering 2.x
    try:
        ui.select_rule("diann", "fragment", _V23, headers)
        raise AssertionError("expected ValueError (no rule covers version)")
    except ValueError as exc:
        assert "no rule covers" in str(exc)
    # columns missing for the version-selected rule → mismatch error
    try:
        ui.select_rule("diann", "protein", _V23, headers - {"PG.MaxLFQ"})
        raise AssertionError("expected ValueError (columns don't match)")
    except ValueError as exc:
        assert "don't match" in str(exc)
