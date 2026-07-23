"""Version-folder rule resolution + column validation.

DIA-NN report columns vary by version, so version-dependent levels live in version subfolders
(``diann/v1/``, ``diann/v2/``) selected by the software version parsed from the param file;
version-agnostic levels stay at the vendor root. Synthetic header sets (derived from the shipped
rules) keep these tests data-free.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import pandas as pd
import pytest
from anndata import AnnData

from anndata_proteomics.converters import pipeline as ui
from anndata_proteomics.rules.loader import (
    load_packaged_rule,
    resolve_rule_for_version,
    resolve_rule_locator,
)
from anndata_proteomics.rules.schema import ParseRule, QuantificationLevel

_V19 = "1.9.2"
_V23 = "2.3.0 Academia "  # messy real catalog string


def _headers_for(rule: ParseRule) -> set[str]:
    """Return a complete synthetic vendor header set for a long-form rule."""
    assert rule.input_shape == "long"
    cols = set(rule.columns.obs.select.values())
    cols.update(rule.columns.var.select.values())
    cols.update(layer.source for layer in rule.layers)
    if rule.fragments is not None and rule.fragments.label_strategy == "column":
        cols.add(rule.fragments.label_column)
    return cols


def _diann_headers(version: str) -> set[str]:
    cols: set[str] = set()
    for level in ui.LEVELS:
        rule = resolve_rule_for_version("diann", level, version)
        if rule is not None:
            cols |= _headers_for(rule)
    return cols


def test_resolve_locator_picks_version_document() -> None:
    protein_v1 = resolve_rule_locator("diann", "protein", _V19)
    protein_v2 = resolve_rule_locator("diann", "protein", _V23)
    fragment_v1 = resolve_rule_locator("diann", "fragment", _V19)
    ion_v1 = resolve_rule_locator("diann", "ion", _V19)
    ion_v2 = resolve_rule_locator("diann", "ion", _V23)
    assert protein_v1 is not None
    assert protein_v2 is not None
    assert fragment_v1 is not None
    assert ion_v1 is not None
    assert ion_v2 is not None
    assert protein_v1.path.parent.name == "v1"
    assert protein_v2.path.parent.name == "v2"
    assert fragment_v1.path.parent.name == "v1"
    assert resolve_rule_locator("diann", "fragment", _V23) is None
    assert ion_v1.path.parent.name == "v1"
    assert ion_v2.path.parent.name == "v2"


def test_rule_software_version_regex_must_match_params_version() -> None:
    ion_v1 = resolve_rule_for_version("diann", "ion", _V19)
    protein_v1 = resolve_rule_for_version("diann", "protein", "1.9.2")
    protein_v2 = resolve_rule_for_version("diann", "protein", "2.3.0 Academia ")
    assert ion_v1 is not None
    assert protein_v1 is not None
    assert protein_v2 is not None
    assert ion_v1.software_version == "^1\\..*"
    assert resolve_rule_for_version("diann", "ion", "3.0.0") is None
    assert protein_v1.software_version == "^1\\..*"
    assert protein_v2.software_version == "^2\\..*"
    with pytest.raises(ValueError, match="no packaged rule"):
        load_packaged_rule("diann", "ion", "3.0.0")


def test_resolve_flat_vendor_and_unknown() -> None:
    maxquant = resolve_rule_locator("maxquant", "ion", None)
    assert maxquant is not None
    assert maxquant.path.name == "rules.json"
    assert resolve_rule_locator("nope", "ion", "1.0") is None


def test_protein_variants_differ_by_version() -> None:
    v1 = resolve_rule_for_version("diann", "protein", _V19)
    v2 = resolve_rule_for_version("diann", "protein", _V23)
    assert v1 is not None
    assert v2 is not None
    v1_layers = {layer.name for layer in v1.layers}
    v2_layers = {layer.name for layer in v2.layers}
    assert "PG_Normalised" in v1_layers
    assert "PG_Normalised" not in v2_layers  # dropped in DIA-NN 2.x


def test_fragment_v1_is_positional() -> None:
    frag = resolve_rule_for_version("diann", "fragment", _V19)
    assert frag is not None
    assert frag.fragments is not None
    assert frag.fragments.label_strategy == "positional"  # positional labels (no Fragment.Info)


def test_convertible_levels_by_version() -> None:
    assert ui.convertible_levels("diann", _V19, _diann_headers(_V19)) == [
        "ion",
        "protein",
        "fragment",
    ]
    assert ui.convertible_levels("diann", _V23, _diann_headers(_V23)) == [
        "ion",
        "protein",
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


def test_convert_level_passes_params_path(monkeypatch: pytest.MonkeyPatch) -> None:
    # convert_level must thread params_path through to converters.assemble.convert.
    import numpy as np

    from anndata_proteomics.converters import assemble

    captured: dict[str, str | Path | None] = {}

    def fake_convert(
        df: pd.DataFrame,
        rule: ParseRule,
        *,
        params_path: str | Path | None = None,
    ) -> AnnData:
        captured["params_path"] = params_path
        return AnnData(
            X=np.array([[1.0]]),
            obs=pd.DataFrame(index=["run1"]),
            var=pd.DataFrame(index=["feature1"]),
        )

    def fake_select_rule(
        slug: str,
        level: QuantificationLevel,
        version: str | None,
        headers: Iterable[str],
    ) -> ParseRule:
        return load_packaged_rule("diann", "ion", "1.9.2")

    monkeypatch.setattr(ui, "select_rule", fake_select_rule)
    monkeypatch.setattr(assemble, "convert", fake_convert)

    adata = ui.convert_level(
        pd.DataFrame({"x": [1]}),
        "diann",
        "ion",
        "1.9.2",
        params_path="/tmp/param_0..txt",
    )

    assert adata.shape == (1, 1)
    assert captured["params_path"] == "/tmp/param_0..txt"
