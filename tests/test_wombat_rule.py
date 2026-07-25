"""Regression tests for the packaged WOMBAT ion rule."""

from __future__ import annotations

import pandas as pd
from numpy.testing import assert_array_equal

from anndata_proteomics.converters.assemble import convert
from anndata_proteomics.converters.recognize import matches
from anndata_proteomics.rules.loader import load_packaged_rule


def test_wombat_preserves_charge_states_as_distinct_ions() -> None:
    rule = load_packaged_rule("wombat", "ion", "0.9.11")
    frame = pd.DataFrame(
        {
            "modified_peptide": [
                "[Acetyl]-AAAAAAAGDSDSWDADAFSVEDPVR",
                "[Acetyl]-AAAAAAAGDSDSWDADAFSVEDPVR",
            ],
            "protein_group": ["O75822", "O75822"],
            "charge": [2, 3],
            "abundance_A_1": [10.0, 20.0],
            "abundance_A_2": [11.0, 21.0],
        }
    )

    result = convert(frame, rule)

    assert result.shape == (2, 2)
    assert result.var_names.tolist() == [
        "[UNIMOD:1]-AAAAAAAGDSDSWDADAFSVEDPVR/2",
        "[UNIMOD:1]-AAAAAAAGDSDSWDADAFSVEDPVR/3",
    ]
    assert result.var["charge"].tolist() == [2, 3]


def test_wombat_levels_match_their_native_output_schemas() -> None:
    ion = load_packaged_rule("wombat", "ion", "0.9.11")
    peptidoform = load_packaged_rule("wombat", "peptidoform", "0.9.11")
    ion_headers = {
        "modified_peptide",
        "protein_group",
        "charge",
        "abundance_A_1",
    }
    peptidoform_headers = {
        "modified_peptide",
        "protein_group",
        "number_of_psms_A_1",
        "abundance_A_1",
    }

    assert matches(ion_headers, ion)
    assert not matches(ion_headers, peptidoform)
    assert matches(peptidoform_headers, peptidoform)
    assert not matches(peptidoform_headers, ion)


def test_wombat_peptidoform_uses_reported_abundance_and_psm_layers() -> None:
    rule = load_packaged_rule("wombat", "peptidoform", "0.9.11")
    frame = pd.DataFrame(
        {
            "modified_peptide": ["PEPTIDE"],
            "protein_group": ["P12345"],
            "abundance_A_1": [10.0],
            "number_of_psms_A_1": [3],
        }
    )

    result = convert(frame, rule)

    assert result.var_names.tolist() == ["PEPTIDE"]
    assert_array_equal(result.X, [[10.0]])
    assert_array_equal(result.layers["Number_Of_Psms"], [[3.0]])
