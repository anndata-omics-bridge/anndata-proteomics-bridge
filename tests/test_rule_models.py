"""Tests for the pydantic ParseRule schema."""

from __future__ import annotations

import tomllib

import pytest
from pydantic import ValidationError

from anndata_proteomics.rules.schema import ParseRule


LONG_EXAMPLE = """
schema_version = "0.1"
file_version = "1"
software_name = "DIA-NN"
software_version = "1.9.1"
input_shape = "long"
quantification_level = "ion"

[axis]
obs_keys = ["Run"]
var_keys = ["ProForma_ion"]
x_layer = "Precursor_Normalised"

[axis.duplicates]
mode = "error"

[columns.obs.select]
File_Name = "File.Name"
Run = "Run"

[columns.var.select]
Modified_Sequence = "Modified.Sequence"
Protein_Ids = "Protein.Ids"
Precursor_Charge = "Precursor.Charge"
Genes = "Genes"

[[columns.var.compute]]
name = "ProForma_peptidoform"
from = ["Modified_Sequence"]
how = "proforma_sequence"

[[columns.var.compute]]
name = "ProForma_ion"
from = ["ProForma_peptidoform", "Precursor_Charge"]
how = "proforma_ion"

[[layers]]
name = "Precursor_Normalised"
source = "Precursor.Normalised"

[[layers]]
name = "Q_Value"
source = "Q.Value"

[[layers]]
name = "RT"
source = "RT"

[[layers]]
name = "Ms1_Area"
source = "Ms1.Area"

[modifications]
source_column = "Modified.Sequence"
parser = "already_proforma"
output_column = "proforma_sequence"
"""


WIDE_EXAMPLE = """
schema_version = "0.1"
file_version = "1"
software_name = "FragPipe"
software_version = "23.0"
input_shape = "wide"
quantification_level = "ion"

[axis]
obs_keys = ["sample"]
var_keys = ["Modified_Sequence", "Charge"]
x_layer = "Intensity"

[axis.duplicates]
mode = "error"

[columns.obs.select]
sample = "<sample>"

[columns.var.select]
Peptide_Sequence = "Peptide Sequence"
Modified_Sequence = "Modified Sequence"
Charge = "Charge"
Protein_ID = "Protein ID"
Gene = "Gene"

[[layers]]
name = "Intensity"
source = "^(?P<sample>.+) Intensity$"

[[layers]]
name = "Spectral_Count"
source = "^(?P<sample>.+) Spectral Count$"

[[layers]]
name = "Match_Type"
source = "^(?P<sample>.+) Match Type$"
encoding_mode = "factor"
categories = { "unmatched" = 0, "MS/MS" = 1, "MBR" = 2 }

[[layers]]
name = "Localization"
source = "^(?P<sample>.+) Localization$"
encoding_mode = "factor"
categories = { "Localized" = 1, "Ambiguous" = 0 }

[sample_name_cleanup]
pattern = ""
"""


def _parse(toml_str: str) -> ParseRule:
    return ParseRule.model_validate(tomllib.loads(toml_str))


def test_long_example_validates():
    rule = _parse(LONG_EXAMPLE)
    assert rule.input_shape == "long"
    assert rule.software_name == "DIA-NN"
    assert len(rule.layers) == 4
    assert rule.axis.x_layer == "Precursor_Normalised"
    assert rule.axis.duplicates.mode == "error"


def test_wide_example_validates():
    rule = _parse(WIDE_EXAMPLE)
    assert rule.input_shape == "wide"
    match_type = next(layer for layer in rule.layers if layer.name == "Match_Type")
    assert match_type.encoding_mode == "factor"
    assert match_type.categories == {"unmatched": 0, "MS/MS": 1, "MBR": 2}


def test_long_layer_missing_source():
    bad = LONG_EXAMPLE.replace('source = "RT"', "")
    with pytest.raises(ValidationError, match="source"):
        _parse(bad)


def test_layer_column_pattern_is_unknown_field():
    # column_pattern was removed from the model; it is now just an extra (forbidden) key.
    bad = LONG_EXAMPLE.replace(
        'name = "RT"\nsource = "RT"',
        'name = "RT"\nsource = "RT"\ncolumn_pattern = "^.+ RT$"',
    )
    with pytest.raises(ValidationError, match="column_pattern"):
        _parse(bad)


def test_wide_layer_missing_source():
    bad = WIDE_EXAMPLE.replace('source = "^(?P<sample>.+) Intensity$"', "")
    with pytest.raises(ValidationError, match="source"):
        _parse(bad)


def test_layer_source_column_is_unknown_field():
    # source_column was removed from the layer model; it is now just an extra (forbidden) key.
    bad = WIDE_EXAMPLE.replace(
        'name = "Intensity"\nsource = "^(?P<sample>.+) Intensity$"',
        'name = "Intensity"\nsource = "^(?P<sample>.+) Intensity$"\nsource_column = "Intensity"',
    )
    with pytest.raises(ValidationError, match="source_column"):
        _parse(bad)


def test_wide_source_requires_sample_group():
    bad = WIDE_EXAMPLE.replace('source = "^(?P<sample>.+) Intensity$"', 'source = "^.+ Intensity$"')
    with pytest.raises(ValidationError, match="sample"):
        _parse(bad)


def test_factor_requires_categories():
    bad = WIDE_EXAMPLE.replace('categories = { "unmatched" = 0, "MS/MS" = 1, "MBR" = 2 }', "")
    with pytest.raises(ValidationError, match="categories"):
        _parse(bad)


def test_x_layer_must_exist():
    bad = LONG_EXAMPLE.replace('x_layer = "Precursor_Normalised"', 'x_layer = "DoesNotExist"')
    with pytest.raises(ValidationError, match="x_layer"):
        _parse(bad)


def test_invalid_duplicates_mode():
    bad = LONG_EXAMPLE.replace('mode = "error"', 'mode = "wrong"')
    with pytest.raises(ValidationError):
        _parse(bad)


def test_top_level_duplicates_rejected():
    bad = LONG_EXAMPLE + '\n[duplicates]\nmode = "error"\n'
    with pytest.raises(ValidationError, match="duplicates"):
        _parse(bad)


def test_unknown_top_level_key_rejected():
    bad = LONG_EXAMPLE + '\nfoo = "bar"\n'
    with pytest.raises(ValidationError, match="foo"):
        _parse(bad)


def test_sample_name_cleanup_rejected_for_long():
    bad = LONG_EXAMPLE + '\n[sample_name_cleanup]\npattern = "(.+)"\n'
    with pytest.raises(ValidationError, match="sample_name_cleanup"):
        _parse(bad)


def test_json_schema_export_has_expected_top_level_properties():
    schema = ParseRule.model_json_schema()
    expected = {
        "schema_version",
        "file_version",
        "software_name",
        "software_version",
        "input_shape",
        "quantification_level",
        "axis",
        "columns",
        "layers",
        "sample_name_cleanup",
        "modifications",
        "fragments",
    }
    assert set(schema["properties"]) == expected


def test_invalid_quantification_level():
    bad = LONG_EXAMPLE.replace('quantification_level = "ion"', 'quantification_level = "wrong"')
    with pytest.raises(ValidationError):
        _parse(bad)


def test_fragment_level_can_be_native_row_level_without_fragments_block():
    good = LONG_EXAMPLE.replace('quantification_level = "ion"', 'quantification_level = "fragment"')
    rule = _parse(good)
    assert rule.quantification_level == "fragment"
    assert rule.fragments is None


def test_missing_quantification_level():
    bad = LONG_EXAMPLE.replace('quantification_level = "ion"\n', "")
    with pytest.raises(ValidationError, match="quantification_level"):
        _parse(bad)


def test_missing_software_version():
    bad = LONG_EXAMPLE.replace('software_version = "1.9.1"\n', "")
    with pytest.raises(ValidationError, match="software_version"):
        _parse(bad)


def test_proforma_ion_requires_two_sources():
    bad = LONG_EXAMPLE.replace(
        'from = ["ProForma_peptidoform", "Precursor_Charge"]',
        'from = ["ProForma_peptidoform"]',
    )
    with pytest.raises(ValidationError, match="exactly two"):
        _parse(bad)


def test_proforma_ion_must_be_var_axis_key():
    bad = LONG_EXAMPLE.replace('var_keys = ["ProForma_ion"]', 'var_keys = ["ProForma_peptidoform"]')
    with pytest.raises(ValidationError, match="axis.var_keys"):
        _parse(bad)


def test_proforma_compute_names_are_pinned():
    bad = LONG_EXAMPLE.replace('name = "ProForma_peptidoform"', 'name = "MyPeptidoform"')
    with pytest.raises(ValidationError, match="ProForma_peptidoform"):
        _parse(bad)


def test_apb_derived_columns_cannot_be_selected():
    bad = LONG_EXAMPLE.replace(
        'Modified_Sequence = "Modified.Sequence"',
        'Modified_Sequence = "Modified.Sequence"\nBad = "proforma_sequence"',
    )
    with pytest.raises(ValidationError, match="derived"):
        _parse(bad)


def test_proforma_sequence_compute_requires_modifications():
    bad = LONG_EXAMPLE.replace(
        '\n[modifications]\nsource_column = "Modified.Sequence"\n'
        'parser = "already_proforma"\noutput_column = "proforma_sequence"\n',
        "\n",
    )
    with pytest.raises(ValidationError, match="modifications"):
        _parse(bad)
