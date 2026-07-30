"""Tests for ``columns.*.optional_select``: configuration-dependent vendor columns.

``select`` is required and gates recognition; ``optional_select`` is captured when the
vendor export carries the column and skipped when it does not, mirroring
``Layer.required = false``.
"""

from __future__ import annotations

import copy
from typing import Any

import pandas as pd
import pytest
from pydantic import ValidationError
from test_rule_models import LONG_EXAMPLE, WIDE_EXAMPLE

from anndata_proteomics.converters.assemble import convert
from anndata_proteomics.converters.recognize import matches
from anndata_proteomics.rules.schema import ParseRule

_INPUT = pd.DataFrame(
    {
        "Run": ["r1", "r1", "r2", "r2"],
        "File.Name": ["/data/r1.raw", "/data/r1.raw", "/data/r2.raw", "/data/r2.raw"],
        "Modified.Sequence": ["PEPTIDE", "PEPTIK", "PEPTIDE", "PEPTIK"],
        "Protein.Ids": ["P1", "P2", "P1", "P2"],
        "Precursor.Charge": [2, 2, 2, 2],
        "Genes": ["G1", "G2", "G1", "G2"],
        "Precursor.Normalised": [10.0, 20.0, 30.0, 40.0],
    }
)


def _optional_rule(**optional_var: str) -> dict[str, Any]:
    """LONG_EXAMPLE with an optional var group, and no packed-list computes to worry about."""
    document = copy.deepcopy(LONG_EXAMPLE)
    document["columns"]["var"]["optional_select"] = dict(optional_var)
    return document


def test_optional_select_column_is_captured_when_present() -> None:
    rule = ParseRule.model_validate(_optional_rule(Lib_Q_Value="Lib.Q.Value"))
    frame = _INPUT.assign(**{"Lib.Q.Value": [0.001, 0.002, 0.001, 0.002]})

    adata = convert(frame, rule)

    assert "Lib_Q_Value" in adata.var.columns


def test_optional_select_column_is_skipped_when_absent() -> None:
    """The same rule converts an export that omits the optional column."""
    rule = ParseRule.model_validate(_optional_rule(Lib_Q_Value="Lib.Q.Value"))

    adata = convert(_INPUT, rule)

    assert "Lib_Q_Value" not in adata.var.columns
    assert adata.shape == (2, 2)


def test_optional_select_does_not_gate_recognition() -> None:
    rule = ParseRule.model_validate(_optional_rule(Lib_Q_Value="Lib.Q.Value"))

    assert matches(_INPUT.columns, rule)


def test_required_select_source_still_raises_when_absent() -> None:
    rule = ParseRule.model_validate(_optional_rule(Lib_Q_Value="Lib.Q.Value"))

    with pytest.raises(ValueError, match="cannot select column 'Genes'"):
        convert(_INPUT.drop(columns=["Genes"]), rule)


def test_coalesce_falls_through_a_skipped_optional_source() -> None:
    document = _optional_rule(Leading_Proteins="Leading.Proteins")
    document["columns"]["var"]["compute"].insert(
        0,
        {
            "name": "Proteins",
            "from": ["Leading_Proteins", "Protein_Ids"],
            "how": "coalesce",
        },
    )
    rule = ParseRule.model_validate(document)

    adata = convert(_INPUT, rule)

    assert list(adata.var["Proteins"]) == ["P1", "P2"]


def test_coalesce_with_every_source_skipped_raises() -> None:
    document = _optional_rule(
        Leading_Proteins="Leading.Proteins",
        Leading_Proteins_Legacy="Leading Proteins",
    )
    document["columns"]["var"]["compute"].insert(
        0,
        {
            "name": "Proteins",
            "from": ["Leading_Proteins", "Leading_Proteins_Legacy"],
            "how": "coalesce",
        },
    )
    rule = ParseRule.model_validate(document)

    with pytest.raises(ValueError, match="every source column is an optional_select"):
        convert(_INPUT, rule)


def test_name_in_both_select_and_optional_select_is_rejected() -> None:
    document = _optional_rule(Genes="Genes.Alternate")

    with pytest.raises(ValidationError, match="declared in both select and optional_select"):
        ParseRule.model_validate(document)


def test_axis_var_key_may_not_be_optional() -> None:
    document = copy.deepcopy(LONG_EXAMPLE)
    document["columns"]["var"]["optional_select"] = {"Precursor_Id": "Precursor.Id"}
    document["axis"]["var_keys"] = ["Precursor_Id"]

    with pytest.raises(ValidationError, match="var_keys must not name optional_select"):
        ParseRule.model_validate(document)


def test_axis_obs_key_may_not_be_optional() -> None:
    document = copy.deepcopy(LONG_EXAMPLE)
    document["columns"]["obs"]["optional_select"] = {"Experiment": "Experiment"}
    document["axis"]["obs_keys"] = ["Experiment"]

    with pytest.raises(ValidationError, match="obs_keys must not name optional_select"):
        ParseRule.model_validate(document)


def test_wide_obs_optional_select_is_rejected() -> None:
    document = copy.deepcopy(WIDE_EXAMPLE)
    document["columns"]["obs"]["optional_select"] = {"Fraction": "Fraction"}

    with pytest.raises(ValidationError, match="not valid for wide rules"):
        ParseRule.model_validate(document)


def test_derived_modification_column_may_not_be_optionally_selected() -> None:
    document = copy.deepcopy(LONG_EXAMPLE)
    document["columns"]["var"]["optional_select"] = {"Sneaky": "proforma_sequence"}

    with pytest.raises(ValidationError, match="must be declared in columns.var.compute"):
        ParseRule.model_validate(document)
