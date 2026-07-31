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
from anndata_proteomics.params.model import Parameters
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


def test_compound_parameters_resolve_rule_version_from_quantification_software() -> None:
    resolution = ui.ParameterResolution(
        source_path=Path("fragpipe.workflow"),
        parameters=Parameters(
            software_name="FragPipe",
            software_version="24.0",
            quantification_software="DIA-NN",
            quantification_software_version="1.8.2 beta 8",
        ),
        version="24.0",
        version_status="present",
    )

    assert ui.resolve_rule_version(resolution, "fragpipe") == ("24.0", "present")
    assert ui.resolve_rule_version(resolution, "diann") == ("1.8.2 beta 8", "present")
    assert ui.resolve_rule_version(resolution, "spectronaut") == (None, "missing")


def test_effective_rule_version_decision_table() -> None:
    mismatch = ui.ParameterResolution(
        source_path=Path("fragpipe.workflow"),
        parameters=Parameters(
            software_name="FragPipe",
            software_version="24.0",
        ),
        version="24.0",
        version_status="present",
    )

    assert ui._effective_rule_version("diann", None, None) == (None, None)
    assert ui._effective_rule_version("diann", "1.9.2", mismatch) == (
        "1.9.2",
        "present",
    )
    assert ui._effective_rule_version("diann", None, mismatch) == (None, "missing")


@pytest.mark.parametrize("version", ["22.0", "22.1-build02", "23.0"])
def test_fragpipe_known_major_versions_resolve(version: str) -> None:
    assert resolve_rule_for_version("fragpipe", "ion", version) is not None


def test_fragpipe_unknown_major_stays_uncovered() -> None:
    assert resolve_rule_for_version("fragpipe", "ion", "24.0") is None


@pytest.mark.parametrize("version", ["13", "13 20250515", "13 20250520", "13.1"])
def test_peaks_known_major_versions_resolve(version: str) -> None:
    assert resolve_rule_for_version("peaks", "ion", version) is not None


def test_peaks_unknown_major_stays_uncovered() -> None:
    assert resolve_rule_for_version("peaks", "ion", "14.0") is None


def test_genuinely_missing_version_selects_unique_rule_by_columns() -> None:
    rule = load_packaged_rule("peaks", "ion")
    headers = set(rule.columns.var.select.values())
    headers.add("LFQ_Run_1 Normalized Area")

    selected = ui.select_rule(
        "peaks",
        "ion",
        None,
        headers,
        version_status="missing",
    )

    assert selected.software_name == "PEAKS"


def test_missing_version_rejects_zero_or_multiple_column_matches(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(ui, "_column_matching_rule_variants", lambda *_args: [])
    with pytest.raises(ValueError, match="no rule matches"):
        ui.select_rule("peaks", "ion", None, [], version_status="missing")

    rule = load_packaged_rule("peaks", "ion")
    monkeypatch.setattr(
        ui,
        "_column_matching_rule_variants",
        lambda *_args: [rule, rule],
    )
    with pytest.raises(ValueError, match="2 rules match"):
        ui.select_rule("peaks", "ion", None, [], version_status="missing")


def test_missing_version_enumerates_multiple_matching_packaged_locators(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rule = load_packaged_rule("peaks", "ion")
    locator = resolve_rule_locator("peaks", "ion", "13")
    assert locator is not None
    monkeypatch.setattr(ui, "iter_packaged_rules", lambda: iter([locator, locator]))
    monkeypatch.setattr(
        ui,
        "load_rule",
        lambda _locator, search_parameters=None: rule,
    )
    headers = set(rule.columns.var.select.values()) | {"LFQ_Run_1 Normalized Area"}

    with pytest.raises(ValueError, match="2 rules match"):
        ui.select_rule("peaks", "ion", None, headers, version_status="missing")

    monkeypatch.setattr(ui, "iter_packaged_rules", lambda: iter([locator]))
    with pytest.raises(ValueError, match="no rule matches"):
        ui.select_rule("peaks", "ion", None, [], version_status="missing")


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


def test_convertible_levels_filters_only_expected_rule_unavailability(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    ion = load_packaged_rule("diann", "ion", _V19)

    def select(
        _slug: str,
        level: QuantificationLevel,
        _version: str | None,
        _headers: Iterable[str],
        **_kwargs: object,
    ) -> tuple[ParseRule, ui.RuleSelectionMethod]:
        if level == "ion":
            return ion, "software_version"
        raise ui.RuleUnavailableError(f"{level} is unavailable")

    monkeypatch.setattr(ui, "_select_rule", select)

    assert ui.convertible_levels("diann", _V19, _diann_headers(_V19)) == ["ion"]


def test_convertible_levels_propagates_unexpected_rule_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail(*_args: object, **_kwargs: object) -> None:
        raise ValueError("malformed packaged rule")

    monkeypatch.setattr(ui, "resolve_rule_for_version", fail)

    with pytest.raises(ValueError, match="malformed packaged rule"):
        ui.convertible_levels("diann", _V19, _diann_headers(_V19))


def test_pipeline_rule_selection_applies_search_parameter_override() -> None:
    headers = _diann_headers(_V23)

    dda = ui.select_rule(
        "diann",
        "ion",
        _V23,
        headers,
        search_parameters=Parameters(acquisition_method="DDA"),
    )
    dia = ui.select_rule(
        "diann",
        "ion",
        _V23,
        headers,
        search_parameters=Parameters(acquisition_method="DIA"),
    )

    assert dda.axis.x_layer == "Ms1_Normalised"
    assert dia.axis.x_layer == "Precursor_Normalised"


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


def test_convert_level_materializes_rule_from_params_once(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import numpy as np

    from anndata_proteomics.converters import assemble
    from anndata_proteomics.params.model import Parameters

    captured: dict[str, object] = {}

    def fake_convert(
        df: pd.DataFrame,
        rule: ParseRule,
        *,
        params_path: str | Path | None = None,
        strict: bool = False,
    ) -> AnnData:
        captured["params_path"] = params_path
        captured["strict"] = strict
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
        *,
        version_status: str | None = None,
        search_parameters: Parameters | None = None,
    ) -> tuple[ParseRule, str]:
        captured["search_parameters"] = search_parameters
        return load_packaged_rule("diann", "ion", "1.9.2"), "software_version"

    resolution = ui.ParameterResolution(
        source_path=Path("/tmp/param_0..txt"),
        parameters=Parameters(
            software_version="1.9.2",
            acquisition_method="DDA",
        ),
        version="1.9.2",
        version_status="present",
    )
    monkeypatch.setattr(ui, "resolve_parameters", lambda *_args: resolution)
    monkeypatch.setattr(ui, "_select_rule", fake_select_rule)
    monkeypatch.setattr(assemble, "convert", fake_convert)

    adata = ui.convert_level(
        pd.DataFrame({"x": [1]}),
        "diann",
        "ion",
        "1.9.2",
        params_path="/tmp/param_0..txt",
    )

    assert adata.shape == (1, 1)
    assert captured["params_path"] is None
    assert captured["search_parameters"] == resolution.parameters
