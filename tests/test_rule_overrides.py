"""Tests for search-parameter-conditional rule-axis overrides."""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest
from pydantic import ValidationError

from anndata_proteomics.params.model import Parameters
from anndata_proteomics.rules.loader import (
    load_parameterized_packaged_rule_for_version,
    load_parameterized_rules,
    load_single_parameterized_rule,
    parse_rule_document,
    resolve_parameterized_rule_for_version,
)
from anndata_proteomics.rules.schema import (
    ParseRule,
    ParseRuleDocument,
    SearchParameterCondition,
    _search_parameter_conditions_are_compatible,
)


def _document(
    overrides: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    return {
        "schema_version": "0.1",
        "file_version": "1",
        "software_name": "Synthetic",
        "software_version": "^1$",
        "base": {
            "input_shape": "long",
            "axis": {"obs_keys": ["Run"], "duplicates": {"mode": "error"}},
            "columns": {"obs": {"select": {"Run": "Run"}}},
        },
        "levels": {
            "ion": {
                "axis": {"var_keys": ["Ion"], "x_layer": "Intensity"},
                "columns": {
                    "var": {
                        "select": {
                            "Ion": "Ion",
                            "OtherIon": "OtherIon",
                        }
                    }
                },
                "layers": [
                    {"name": "Intensity", "source": "Intensity"},
                    {"name": "Alternate", "source": "Alternate"},
                ],
                "search_parameter_overrides": overrides or [],
            }
        },
    }


def _parse(data: dict[str, Any]) -> ParseRuleDocument:
    return parse_rule_document(json.dumps(data))


def test_override_selects_axis_from_normalized_search_parameters() -> None:
    document = _parse(
        _document(
            [
                {
                    "when_search_parameters": {"enable_match_between_runs": "yes"},
                    "axis": {"x_layer": "Alternate"},
                }
            ]
        )
    )

    override = document.levels["ion"].search_parameter_overrides[0]
    assert isinstance(override.when_search_parameters, SearchParameterCondition)
    assert override.when_search_parameters.model_dump(exclude_unset=True) == {
        "enable_match_between_runs": True
    }
    assert document.effective_rule("ion").axis.x_layer == "Intensity"
    assert document.parameterized_effective_rule("ion", Parameters()).axis.x_layer == "Intensity"
    assert (
        document.parameterized_effective_rule(
            "ion",
            Parameters(enable_match_between_runs=True),
        ).axis.x_layer
        == "Alternate"
    )


def test_all_matching_overrides_apply_in_source_order() -> None:
    document = _parse(
        _document(
            [
                {
                    "when_search_parameters": {"acquisition_method": "DDA"},
                    "axis": {"x_layer": "Alternate"},
                },
                {
                    "when_search_parameters": {
                        "acquisition_method": "DDA",
                        "enable_match_between_runs": True,
                    },
                    "axis": {"x_layer": "Intensity"},
                },
            ]
        )
    )

    assert (
        document.parameterized_effective_rule(
            "ion",
            Parameters(acquisition_method="DDA"),
        ).axis.x_layer
        == "Alternate"
    )
    assert (
        document.parameterized_effective_rule(
            "ion",
            Parameters(
                acquisition_method="DDA",
                enable_match_between_runs=True,
            ),
        ).axis.x_layer
        == "Intensity"
    )
    assert (
        document.parameterized_effective_rules(Parameters(acquisition_method="DIA"))[
            "ion"
        ].axis.x_layer
        == "Intensity"
    )


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (
            lambda data: data["base"].update({"search_parameter_overrides": []}),
            "extra_forbidden",
        ),
        (
            lambda data: data["levels"]["ion"]["search_parameter_overrides"].append(
                {
                    "when_search_parameters": {},
                    "axis": {"x_layer": "Alternate"},
                }
            ),
            "must declare at least one field",
        ),
        (
            lambda data: data["levels"]["ion"]["search_parameter_overrides"].append(
                {
                    "when_search_parameters": {"not_a_parameter": "value"},
                    "axis": {"x_layer": "Alternate"},
                }
            ),
            "extra_forbidden",
        ),
        (
            lambda data: data["levels"]["ion"]["search_parameter_overrides"].append(
                {
                    "when_search_parameters": {"acquisition_method": "SWATH"},
                    "axis": {"x_layer": "Alternate"},
                }
            ),
            "acquisition_method",
        ),
        (
            lambda data: data["levels"]["ion"]["search_parameter_overrides"].append(
                {
                    "when_search_parameters": {"acquisition_method": "DDA"},
                    "axis": {},
                }
            ),
            "must declare at least one field",
        ),
        (
            lambda data: data["levels"]["ion"]["search_parameter_overrides"].append(
                {
                    "when_search_parameters": {"acquisition_method": "DDA"},
                    "axis": {"x_layer": "Alternate"},
                    "layers": [],
                }
            ),
            "extra_forbidden",
        ),
    ],
)
def test_override_schema_rejects_invalid_or_recursive_shapes(
    mutation: Callable[[dict[str, Any]], None],
    message: str,
) -> None:
    data = _document()
    mutation(data)

    with pytest.raises(ValidationError, match=message):
        _parse(data)


def test_document_validation_checks_conditional_axis() -> None:
    data = _document(
        [
            {
                "when_search_parameters": {"acquisition_method": "DDA"},
                "axis": {"x_layer": "Missing"},
            }
        ]
    )

    with pytest.raises(ValidationError, match="x_layer"):
        _parse(data)


def test_cross_field_invalid_condition_combination_is_unreachable() -> None:
    data = _document(
        [
            {
                "when_search_parameters": {"min_precursor_mz": 900},
                "axis": {"x_layer": "Alternate"},
            },
            {
                "when_search_parameters": {"max_precursor_mz": 300},
                "axis": {"var_keys": ["OtherIon"]},
            },
        ]
    )

    document = _parse(data)
    overrides = tuple(document.levels["ion"].search_parameter_overrides)
    assert not _search_parameter_conditions_are_compatible(overrides)
    assert document.effective_rule("ion").axis.x_layer == "Intensity"


def test_conflicting_equalities_are_incompatible() -> None:
    document = _parse(
        _document(
            [
                {
                    "when_search_parameters": {"acquisition_method": "DDA"},
                    "axis": {"x_layer": "Alternate"},
                },
                {
                    "when_search_parameters": {"acquisition_method": "DIA"},
                    "axis": {"x_layer": "Intensity"},
                },
            ]
        )
    )
    overrides = tuple(document.levels["ion"].search_parameter_overrides)

    assert not _search_parameter_conditions_are_compatible(overrides)
    with pytest.raises(KeyError, match="protein"):
        document.validate_effective_rule_variants("protein")


def test_loader_functions_propagate_search_parameters(tmp_path: Path) -> None:
    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            _document(
                [
                    {
                        "when_search_parameters": {"acquisition_method": "DDA"},
                        "axis": {"x_layer": "Alternate"},
                    }
                ]
            )
        ),
        encoding="utf-8",
    )
    parameters = Parameters(acquisition_method="DDA")

    assert load_single_parameterized_rule(path, parameters).axis.x_layer == "Alternate"
    assert load_parameterized_rules(path, parameters)["ion"].axis.x_layer == "Alternate"


def test_packaged_loader_functions_propagate_search_parameters() -> None:
    parameters = Parameters(acquisition_method="DDA")

    packaged = load_parameterized_packaged_rule_for_version(
        "diann",
        "ion",
        "2.6.0",
        parameters,
    )
    resolved = resolve_parameterized_rule_for_version(
        "diann",
        "ion",
        "2.6.0",
        parameters,
    )

    assert packaged.axis.x_layer == "Ms1_Normalised"
    assert isinstance(resolved, ParseRule)
    assert resolved.axis.x_layer == "Ms1_Normalised"


def _gated_document(
    ion_requires: dict[str, Any] | None = None,
    extra_levels: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """The synthetic document with an availability gate on its ion level."""
    data = _document()
    if ion_requires is not None:
        data["levels"]["ion"]["requires_search_parameters"] = ion_requires
    if extra_levels:
        data["levels"].update(extra_levels)
    return data


def test_requires_search_parameters_normalizes_and_gates_availability() -> None:
    document = _parse(_gated_document({"combine_charge_states": False}))

    requirements = document.levels["ion"].requires_search_parameters
    assert isinstance(requirements, SearchParameterCondition)
    assert requirements.model_dump(exclude_unset=True) == {"combine_charge_states": False}
    assert document.level_is_available_for("ion", Parameters(combine_charge_states=False))
    assert not document.level_is_available_for("ion", Parameters(combine_charge_states=True))
    # A gate needs parameters to compare against; without them the level is not offered.
    assert not document.level_is_unconditionally_available("ion")


def test_ungated_level_is_always_available() -> None:
    document = _parse(_document())

    assert document.level_is_unconditionally_available("ion")
    assert document.level_is_available_for("ion", Parameters(combine_charge_states=True))


def test_level_is_available_is_false_for_an_undeclared_level() -> None:
    document = _parse(_document())

    assert not document.level_is_available_for("protein", Parameters())


def test_available_levels_filters_by_gate_in_source_order() -> None:
    document = _parse(_gated_document({"combine_charge_states": False}))

    assert document.available_levels_for(Parameters(combine_charge_states=False)) == ["ion"]
    assert document.available_levels_for(Parameters(combine_charge_states=True)) == []
    assert document.unconditionally_available_levels() == []
    assert _parse(_document()).unconditionally_available_levels() == ["ion"]


def test_requires_search_parameters_rejects_unknown_fields() -> None:
    with pytest.raises(ValidationError, match="extra_forbidden"):
        _parse(_gated_document({"not_a_parameter_field": True}))


def test_requires_search_parameters_does_not_reach_the_effective_rule() -> None:
    """The gate decides availability only; it must not leak into the merged rule body."""
    document = _parse(_gated_document({"combine_charge_states": False}))

    rule = document.parameterized_effective_rule("ion", Parameters(combine_charge_states=False))

    assert not hasattr(rule, "requires_search_parameters")
    assert document.parameterized_effective_rules(
        Parameters(combine_charge_states=True)
    ).keys() == {"ion"}
