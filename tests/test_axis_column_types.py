"""Logical typing contracts for selected observation and feature columns."""

from __future__ import annotations

import copy
import json
from pathlib import Path

import anndata as ad
import mudata
import pandas as pd
import pytest
from pydantic import ValidationError
from test_rule_models import LONG_EXAMPLE

from anndata_proteomics.adapters.anndata import conversion as conversion_adapter
from anndata_proteomics.converters import pipeline as conversion_pipeline
from anndata_proteomics.converters.assemble import convert_table
from anndata_proteomics.rules.loader import parse_rule_document
from anndata_proteomics.rules.schema import ParseRule
from anndata_proteomics.workflows.conversion import LevelConversion


def _typed_long_rule() -> ParseRule:
    document = copy.deepcopy(LONG_EXAMPLE)
    document["columns"]["obs"]["types"] = {"Run": "string"}
    document["columns"]["var"]["types"] = {
        "Precursor_Charge": "integer",
        "Genes": "boolean",
    }
    return ParseRule.model_validate(document)


def _typed_level_conversion() -> LevelConversion:
    rule = _typed_long_rule()
    frame = pd.DataFrame(
        {
            "File.Name": ["run.raw", "run.raw"],
            "Run": ["001", "001"],
            "Modified.Sequence": ["PEPTIDE", "PEPTIDER"],
            "Protein.Ids": ["0007", "9007199254740993"],
            "Precursor.Charge": [2, 3],
            "Genes": [1, None],
            "Precursor.Normalised": [10.0, 20.0],
            "Q.Value": [0.01, 0.02],
            "RT": [1.0, 2.0],
            "Ms1.Area": [3.0, 4.0],
        }
    )
    selection = conversion_pipeline.RuleSelection(rule, "rule_config")
    return LevelConversion("ion", selection, convert_table(frame, rule))


def _assert_logical_values_survive_round_trip(restored: ad.AnnData) -> None:
    assert restored.obs["Run"].astype("string").tolist() == ["001"]
    assert restored.var["Protein_Ids"].astype("string").tolist() == [
        "0007",
        "9007199254740993",
    ]
    assert restored.var["Precursor_Charge"].tolist() == [2, 3]
    genes = restored.var["Genes"].astype("boolean")
    assert bool(genes.iloc[0]) is True
    assert pd.isna(genes.iloc[1])


def test_selected_column_types_default_to_string() -> None:
    rule = _typed_long_rule()

    assert rule.columns.obs.type_for("Run") == "string"
    assert rule.columns.var.type_for("Modified_Sequence") == "string"
    assert rule.columns.var.type_for("Precursor_Charge") == "integer"
    assert rule.columns.var.type_for("Genes") == "boolean"


def test_type_declaration_must_name_selected_output() -> None:
    document = copy.deepcopy(LONG_EXAMPLE)
    document["columns"]["var"]["types"] = {"Missing": "number"}

    with pytest.raises(ValidationError, match="types must name selected columns"):
        ParseRule.model_validate(document)


def test_type_declaration_rejects_unknown_logical_type() -> None:
    document = copy.deepcopy(LONG_EXAMPLE)
    document["columns"]["var"]["types"] = {"Precursor_Charge": "uint8"}

    with pytest.raises(ValidationError, match="string.*integer.*number.*boolean"):
        ParseRule.model_validate(document)


def test_document_composes_base_and_level_column_types() -> None:
    document = {
        "schema_version": "0.2",
        "file_version": "1",
        "software_name": "Synthetic",
        "software_version": "1",
        "base": {
            "input_shape": "long",
            "axis": {"obs_keys": ["Run"]},
            "columns": {
                "obs": {
                    "select": {"Run": "Run"},
                    "types": {"Run": "string"},
                }
            },
        },
        "levels": {
            "ion": {
                "axis": {
                    "var_keys": ["Feature"],
                    "x_layer": "Intensity",
                },
                "columns": {
                    "var": {
                        "select": {
                            "Feature": "Feature",
                            "Score": "Score",
                        },
                        "types": {
                            "Feature": "integer",
                            "Score": "number",
                        },
                    }
                },
                "layers": [{"name": "Intensity", "source": "Intensity"}],
            }
        },
    }

    rule = parse_rule_document(json.dumps(document)).effective_rule("ion")

    assert rule.columns.obs.types == {"Run": "string"}
    assert rule.columns.var.types == {
        "Feature": "integer",
        "Score": "number",
    }


def test_proforma_ion_requires_integer_declared_charge() -> None:
    document = copy.deepcopy(LONG_EXAMPLE)
    document["columns"]["var"]["types"].pop("Precursor_Charge")

    with pytest.raises(ValidationError, match="proforma_ion.*integer"):
        ParseRule.model_validate(document)


def test_selected_columns_are_coerced_before_computed_columns() -> None:
    rule = _typed_long_rule()
    frame = pd.DataFrame(
        {
            "File.Name": ["run.raw", "run.raw"],
            "Run": [1, 1],
            "Modified.Sequence": ["PEPTIDE", "PEPTIDER"],
            "Protein.Ids": ["P1", "P2"],
            "Precursor.Charge": [2.0, 3.0],
            "Genes": [1, 0],
            "Precursor.Normalised": [10.0, 20.0],
            "Q.Value": [0.01, 0.02],
            "RT": [1.0, 2.0],
            "Ms1.Area": [3.0, 4.0],
        }
    )

    pieces = convert_table(frame, rule)

    assert str(pieces.obs["Run"].dtype) == "string"
    assert str(pieces.var["Precursor_Charge"].dtype) == "Int64"
    assert str(pieces.var["Genes"].dtype) == "boolean"
    assert pieces.var["ProForma_ion"].tolist() == ["PEPTIDE/2", "PEPTIDER/3"]


def test_invalid_declared_number_raises_with_column_context() -> None:
    document = copy.deepcopy(LONG_EXAMPLE)
    document["columns"]["var"]["types"] = {
        "Precursor_Charge": "integer",
        "Genes": "number",
    }
    rule = ParseRule.model_validate(document)
    frame = pd.DataFrame(
        {
            "File.Name": ["run.raw"],
            "Run": ["run"],
            "Modified.Sequence": ["PEPTIDE"],
            "Protein.Ids": ["P1"],
            "Precursor.Charge": [2],
            "Genes": ["not-a-number"],
            "Precursor.Normalised": [10.0],
            "Q.Value": [0.01],
            "RT": [1.0],
            "Ms1.Area": [3.0],
        }
    )

    with pytest.raises(
        ValueError,
        match="var column 'Genes'.*vendor source 'Genes'.*logical type 'number'",
    ):
        convert_table(frame, rule)


@pytest.mark.parametrize(
    ("logical_type", "value"),
    [
        ("integer", 2.5),
        ("boolean", "yes"),
    ],
)
def test_invalid_declared_values_never_become_missing(
    logical_type: str,
    value: object,
) -> None:
    document = copy.deepcopy(LONG_EXAMPLE)
    document["columns"]["var"]["types"]["Genes"] = logical_type
    rule = ParseRule.model_validate(document)
    frame = pd.DataFrame(
        {
            "File.Name": ["run.raw"],
            "Run": ["run"],
            "Modified.Sequence": ["PEPTIDE"],
            "Protein.Ids": ["P1"],
            "Precursor.Charge": [2],
            "Genes": [value],
            "Precursor.Normalised": [10.0],
            "Q.Value": [0.01],
            "RT": [1.0],
            "Ms1.Area": [3.0],
        }
    )

    with pytest.raises(ValueError, match=f"logical type '{logical_type}'"):
        convert_table(frame, rule)


def test_integer_direct_key_has_stable_identifier_spelling() -> None:
    rule = ParseRule.model_validate(
        {
            "schema_version": "0.2",
            "file_version": "1",
            "software_name": "Synthetic",
            "software_version": "1",
            "input_shape": "long",
            "quantification_level": "fragment",
            "axis": {
                "obs_keys": ["Run"],
                "var_keys": ["Fragment_Number"],
                "x_layer": "Intensity",
            },
            "columns": {
                "obs": {"select": {"Run": "Run"}},
                "var": {
                    "select": {"Fragment_Number": "Fragment.Number"},
                    "types": {"Fragment_Number": "integer"},
                },
            },
            "layers": [{"name": "Intensity", "source": "Intensity"}],
        }
    )
    frame = pd.DataFrame(
        {
            "Run": ["r1", "r1"],
            "Fragment.Number": [1.0, 2.0],
            "Intensity": [10.0, 20.0],
        }
    )

    pieces = convert_table(frame, rule)

    assert pieces.var.index.tolist() == ["1", "2"]
    assert str(pieces.var["Fragment_Number"].dtype) == "Int64"


def test_rule_set_derives_exact_string_sources_for_text_reading() -> None:
    rule = _typed_long_rule()

    sources = conversion_pipeline.string_sources_for_rules((rule,))

    assert sources == frozenset(
        {
            "File.Name",
            "Run",
            "Modified.Sequence",
            "Protein.Ids",
        }
    )


def test_rule_set_rejects_conflicting_types_for_one_vendor_source() -> None:
    integer_rule = _typed_long_rule()
    document = copy.deepcopy(LONG_EXAMPLE)
    document["columns"]["var"]["select"] = {"Alias": "Precursor.Charge"}
    document["columns"]["var"]["types"] = {"Alias": "string"}
    document["columns"]["var"]["compute"] = []
    document["axis"]["var_keys"] = ["Alias"]
    string_rule = ParseRule.model_validate(document)

    with pytest.raises(ValueError, match="conflicting logical types.*Precursor.Charge"):
        conversion_pipeline.string_sources_for_rules((integer_rule, string_rule))


def test_logical_axis_values_survive_h5ad_round_trip(tmp_path: Path) -> None:
    adata = conversion_adapter.to_anndata(_typed_level_conversion())
    path = tmp_path / "typed.h5ad"

    adata.write_h5ad(path)

    _assert_logical_values_survive_round_trip(ad.read_h5ad(path))


def test_logical_axis_values_survive_h5mu_round_trip(tmp_path: Path) -> None:
    container = conversion_adapter.to_mudata({"ion": _typed_level_conversion()})
    path = tmp_path / "typed.h5mu"

    container.write_h5mu(path)
    with mudata.set_options(pull_on_update=False):
        restored = mudata.read_h5mu(path)

    ion = restored.mod["ion"]
    if not isinstance(ion, ad.AnnData):
        raise TypeError(f"ion modality must be AnnData, got {type(ion).__name__}")
    _assert_logical_values_survive_round_trip(ion)
