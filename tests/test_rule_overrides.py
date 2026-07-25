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
    load_packaged_rule,
    load_rule,
    load_rules,
    parse_rule_document,
    resolve_rule_for_version,
)
from anndata_proteomics.rules.schema import (
    ParseRuleDocument,
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
    assert override.when_search_parameters == {"enable_match_between_runs": True}
    assert document.effective_rule("ion").axis.x_layer == "Intensity"
    assert document.effective_rule("ion", Parameters()).axis.x_layer == "Intensity"
    assert (
        document.effective_rule(
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
        document.effective_rule(
            "ion",
            Parameters(acquisition_method="DDA"),
        ).axis.x_layer
        == "Alternate"
    )
    assert (
        document.effective_rule(
            "ion",
            Parameters(
                acquisition_method="DDA",
                enable_match_between_runs=True,
            ),
        ).axis.x_layer
        == "Intensity"
    )
    assert (
        document.effective_rules(Parameters(acquisition_method="DIA"))["ion"].axis.x_layer
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
            "too_short",
        ),
        (
            lambda data: data["levels"]["ion"]["search_parameter_overrides"].append(
                {
                    "when_search_parameters": {"not_a_parameter": "value"},
                    "axis": {"x_layer": "Alternate"},
                }
            ),
            "unknown search-parameter",
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

    assert load_rule(path, search_parameters=parameters).axis.x_layer == "Alternate"
    assert load_rules(path, search_parameters=parameters)["ion"].axis.x_layer == "Alternate"


def test_packaged_loader_functions_propagate_search_parameters() -> None:
    parameters = Parameters(acquisition_method="DDA")

    packaged = load_packaged_rule(
        "diann",
        "ion",
        "2.6.0",
        search_parameters=parameters,
    )
    resolved = resolve_rule_for_version(
        "diann",
        "ion",
        "2.6.0",
        search_parameters=parameters,
    )

    assert packaged.axis.x_layer == "Ms1_Normalised"
    assert resolved is not None
    assert resolved.axis.x_layer == "Ms1_Normalised"
