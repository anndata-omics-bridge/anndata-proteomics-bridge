"""Version-folder rule resolution + column validation.

DIA-NN report columns vary by version, so version-dependent levels live in version subfolders
(``diann/v1/``, ``diann/v2/``) selected by the software version parsed from the param file;
version-agnostic levels stay at the vendor root. Synthetic header sets (derived from the shipped
rules) keep these tests data-free.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from anndata_proteomics.converters import pipeline as ui
from anndata_proteomics.converters._pieces import ConversionPieces
from anndata_proteomics.params.model import Parameters
from anndata_proteomics.rules.loader import (
    PackagedRuleUnavailable,
    RuleLocatorUnavailable,
    load_packaged_rule,
    load_packaged_rule_for_version,
    resolve_rule_for_version,
    resolve_rule_locator_for_version,
    resolve_rule_locator_without_version,
)
from anndata_proteomics.rules.registry import RuleLocator
from anndata_proteomics.rules.schema import ParseRule, QuantificationLevel
from anndata_proteomics.workflows import conversion as conversion_workflow

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
        if isinstance(rule, ParseRule):
            cols |= _headers_for(rule)
    return cols


def test_resolve_locator_picks_version_document() -> None:
    protein_v1 = resolve_rule_locator_for_version("diann", "protein", _V19)
    protein_v2 = resolve_rule_locator_for_version("diann", "protein", _V23)
    fragment_v1 = resolve_rule_locator_for_version("diann", "fragment", _V19)
    ion_v1 = resolve_rule_locator_for_version("diann", "ion", _V19)
    ion_v2 = resolve_rule_locator_for_version("diann", "ion", _V23)
    assert isinstance(protein_v1, RuleLocator)
    assert isinstance(protein_v2, RuleLocator)
    assert isinstance(fragment_v1, RuleLocator)
    assert isinstance(ion_v1, RuleLocator)
    assert isinstance(ion_v2, RuleLocator)
    assert protein_v1.path.parent.name == "v1"
    assert protein_v2.path.parent.name == "v2"
    assert fragment_v1.path.parent.name == "v1"
    unavailable = resolve_rule_locator_for_version("diann", "fragment", _V23)
    assert isinstance(unavailable, RuleLocatorUnavailable)
    assert ion_v1.path.parent.name == "v1"
    assert ion_v2.path.parent.name == "v2"


def test_rule_software_version_regex_must_match_params_version() -> None:
    ion_v1 = resolve_rule_for_version("diann", "ion", _V19)
    protein_v1 = resolve_rule_for_version("diann", "protein", "1.9.2")
    protein_v2 = resolve_rule_for_version("diann", "protein", "2.3.0 Academia ")
    assert isinstance(ion_v1, ParseRule)
    assert isinstance(protein_v1, ParseRule)
    assert isinstance(protein_v2, ParseRule)
    assert ion_v1.software_version == "^1\\..*"
    assert isinstance(
        resolve_rule_for_version("diann", "ion", "3.0.0"),
        PackagedRuleUnavailable,
    )
    assert protein_v1.software_version == "^1\\..*"
    assert protein_v2.software_version == "^2\\..*"
    with pytest.raises(ValueError, match="no packaged rule"):
        load_packaged_rule_for_version("diann", "ion", "3.0.0")


def test_compound_parameters_resolve_rule_version_from_quantification_software() -> None:
    resolution = ui.ParameterResolution(
        source_path=Path("fragpipe.workflow"),
        parameters=Parameters(
            software_name="FragPipe",
            software_version="24.0",
            quantification_software="DIA-NN",
            quantification_software_version="1.8.2 beta 8",
        ),
        version=ui.PresentRuleVersion("24.0"),
    )

    assert ui.resolve_rule_version(resolution, "fragpipe") == ui.PresentRuleVersion("24.0")
    assert ui.resolve_rule_version(resolution, "diann") == ui.PresentRuleVersion("1.8.2 beta 8")
    assert ui.resolve_rule_version(resolution, "spectronaut") == ui.MissingRuleVersion()


def test_mismatched_parameter_software_has_no_applicable_rule_version() -> None:
    mismatch = ui.ParameterResolution(
        source_path=Path("fragpipe.workflow"),
        parameters=Parameters(
            software_name="FragPipe",
            software_version="24.0",
        ),
        version=ui.PresentRuleVersion("24.0"),
    )

    assert ui.resolve_rule_version(mismatch, "diann") == ui.MissingRuleVersion()


@pytest.mark.parametrize("version", ["22.0", "22.1-build02", "23.0"])
def test_fragpipe_known_major_versions_resolve(version: str) -> None:
    assert isinstance(resolve_rule_for_version("fragpipe", "ion", version), ParseRule)


def test_fragpipe_unknown_major_stays_uncovered() -> None:
    assert isinstance(
        resolve_rule_for_version("fragpipe", "ion", "24.0"),
        PackagedRuleUnavailable,
    )


@pytest.mark.parametrize("version", ["13", "13 20250515", "13 20250520", "13.1"])
def test_peaks_known_major_versions_resolve(version: str) -> None:
    assert isinstance(resolve_rule_for_version("peaks", "ion", version), ParseRule)


def test_peaks_unknown_major_stays_uncovered() -> None:
    assert isinstance(
        resolve_rule_for_version("peaks", "ion", "14.0"),
        PackagedRuleUnavailable,
    )


def test_genuinely_missing_version_selects_unique_rule_by_columns() -> None:
    rule = load_packaged_rule("peaks", "ion")
    headers = set(rule.columns.var.select.values())
    headers.add("LFQ_Run_1 Normalized Area")

    selected = ui.select_rule_by_columns(
        "peaks",
        "ion",
        headers,
    )

    assert selected.rule.software_name == "PEAKS"


def test_missing_version_rejects_zero_or_multiple_column_matches(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(ui, "_column_matching_rules", lambda *_args: [])
    with pytest.raises(ValueError, match="no rule matches"):
        ui.select_rule_by_columns("peaks", "ion", [])

    rule = load_packaged_rule("peaks", "ion")
    monkeypatch.setattr(
        ui,
        "_column_matching_rules",
        lambda *_args: [rule, rule],
    )
    with pytest.raises(ValueError, match="2 rules match"):
        ui.select_rule_by_columns("peaks", "ion", [])


def test_missing_version_enumerates_multiple_matching_packaged_locators(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rule = load_packaged_rule("peaks", "ion")
    locator = resolve_rule_locator_for_version("peaks", "ion", "13")
    assert isinstance(locator, RuleLocator)
    monkeypatch.setattr(ui, "iter_packaged_rules", lambda: iter([locator, locator]))
    monkeypatch.setattr(
        ui,
        "load_rule",
        lambda _locator: rule,
    )
    headers = set(rule.columns.var.select.values()) | {"LFQ_Run_1 Normalized Area"}

    with pytest.raises(ValueError, match="2 rules match"):
        ui.select_rule_by_columns("peaks", "ion", headers)

    monkeypatch.setattr(ui, "iter_packaged_rules", lambda: iter([locator]))
    with pytest.raises(ValueError, match="no rule matches"):
        ui.select_rule_by_columns("peaks", "ion", [])


def test_resolve_flat_vendor_and_unknown() -> None:
    maxquant = resolve_rule_locator_without_version("maxquant", "ion")
    assert isinstance(maxquant, RuleLocator)
    assert maxquant.path.name == "rules.json"
    unknown = resolve_rule_locator_for_version("nope", "ion", "1.0")
    assert isinstance(unknown, RuleLocatorUnavailable)


def test_protein_variants_differ_by_version() -> None:
    v1 = resolve_rule_for_version("diann", "protein", _V19)
    v2 = resolve_rule_for_version("diann", "protein", _V23)
    assert isinstance(v1, ParseRule)
    assert isinstance(v2, ParseRule)
    v1_layers = {layer.name for layer in v1.layers}
    v2_layers = {layer.name for layer in v2.layers}
    assert "PG_Normalised" in v1_layers
    assert "PG_Normalised" not in v2_layers  # dropped in DIA-NN 2.x


def test_fragment_v1_is_positional() -> None:
    frag = resolve_rule_for_version("diann", "fragment", _V19)
    assert isinstance(frag, ParseRule)
    assert frag.fragments is not None
    assert frag.fragments.label_strategy == "positional"  # positional labels (no Fragment.Info)


def test_convertible_levels_by_version() -> None:
    assert list(ui.available_rules_for_version("diann", _V19, _diann_headers(_V19))) == [
        "ion",
        "protein",
        "fragment",
    ]
    levels = list(ui.available_rules_for_version("diann", _V23, _diann_headers(_V23)))
    assert levels == [
        "ion",
        "protein",
    ]
    assert levels


def test_convertible_levels_filters_only_expected_rule_unavailability(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    ion = load_packaged_rule_for_version("diann", "ion", _V19)

    def find(
        _slug: str,
        level: QuantificationLevel,
        _version: str,
        _headers: set[str],
    ) -> ui.RuleLookup:
        if level == "ion":
            return ui.RuleSelection(ion, "software_version")
        return ui.RuleUnavailable(f"{level} is unavailable")

    monkeypatch.setattr(ui, "find_rule_for_version", find)

    assert list(ui.available_rules_for_version("diann", _V19, _diann_headers(_V19))) == ["ion"]


def test_convertible_levels_propagates_unexpected_rule_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail(*_args: object, **_kwargs: object) -> None:
        raise ValueError("malformed packaged rule")

    monkeypatch.setattr(ui, "resolve_rule_for_version", fail)

    with pytest.raises(ValueError, match="malformed packaged rule"):
        ui.available_rules_for_version("diann", _V19, _diann_headers(_V19))


def test_pipeline_rule_selection_applies_search_parameter_override() -> None:
    headers = _diann_headers(_V23)

    dda = ui.select_parameterized_rule_for_version(
        "diann",
        "ion",
        _V23,
        headers,
        Parameters(acquisition_method="DDA"),
    )
    dia = ui.select_parameterized_rule_for_version(
        "diann",
        "ion",
        _V23,
        headers,
        Parameters(acquisition_method="DIA"),
    )

    assert dda.rule.axis.x_layer == "Ms1_Normalised"
    assert dia.rule.axis.x_layer == "Precursor_Normalised"


def test_select_rule_errors() -> None:
    headers = _diann_headers(_V23)
    # fragment has no rule covering 2.x
    with pytest.raises(ValueError, match="no rule covers"):
        ui.select_rule_for_version("diann", "fragment", _V23, headers)
    # columns missing for the version-selected rule → mismatch error
    with pytest.raises(ValueError, match="don't match"):
        ui.select_rule_for_version("diann", "protein", _V23, headers - {"PG.MaxLFQ"})


def test_convert_level_materializes_rule_from_params_once(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}
    rule = load_packaged_rule_for_version("diann", "ion", "1.9.2")
    pieces = ConversionPieces(
        X=np.array([[1.0]]),
        obs=pd.DataFrame(index=["run1"]),
        var=pd.DataFrame(index=["feature1"]),
        layers={rule.axis.x_layer: np.array([[1.0]])},
    )

    def fake_select_rule_from_parameters(
        _headers: pd.Index[str],
        _slug: str,
        level: QuantificationLevel,
        resolution: ui.ParameterResolution,
    ) -> ui.RuleSelection:
        captured["level"] = level
        captured["parameters"] = resolution.parameters
        return ui.RuleSelection(rule, "software_version")

    resolution = ui.ParameterResolution(
        source_path=Path("/tmp/param_0..txt"),
        parameters=Parameters(
            software_version="1.9.2",
            acquisition_method="DDA",
        ),
        version=ui.PresentRuleVersion("1.9.2"),
    )
    monkeypatch.setattr(
        conversion_workflow,
        "select_rule_from_parameters",
        fake_select_rule_from_parameters,
    )
    monkeypatch.setattr(conversion_workflow, "convert_table", lambda *_args, **_kwargs: pieces)

    conversion = conversion_workflow.convert_level_from_parameters(
        pd.DataFrame({"x": [1]}),
        "diann",
        "ion",
        resolution,
    )

    assert conversion.pieces is pieces
    assert captured == {"level": "ion", "parameters": resolution.parameters}
