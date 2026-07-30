"""Focused failure and alternate-format paths for vendor parameter parsers."""

from __future__ import annotations

import json
from io import StringIO
from pathlib import Path
from typing import Any, cast

import pandas as pd
import pytest
import yaml

from anndata_proteomics.params.parsers import (
    alphapept,
    diann,
    fragpipe,
    maxquant,
    metamorpheus,
    peaks,
    sage,
    spectronaut,
    wombat,
)


def test_alphapept_and_wombat_non_trypsin_inputs() -> None:
    alpha_document = {
        "summary": {"version": "1"},
        "fasta": {
            "protease": "lys-c",
            "mods_fixed": [],
            "mods_fixed_terminal": [],
            "mods_fixed_terminal_prot": [],
            "mods_variable": [],
            "mods_variable_terminal": [],
            "mods_variable_terminal_prot": [],
            "n_missed_cleavages": 1,
            "n_modifications_max": 2,
            "pep_length_min": 6,
            "pep_length_max": 30,
        },
        "search": {
            "ppm": False,
            "prec_tol": 0.5,
            "frag_tol": 0.2,
            "protein_fdr": 0.01,
            "peptide_fdr": 0.01,
        },
        "features": {"iso_charge_min": 1, "iso_charge_max": 4},
        "workflow": {"match": False},
    }
    alpha = alphapept.extract_params(StringIO(yaml.safe_dump(alpha_document)))
    assert alpha.enzyme == "Lys-C"
    assert alpha.precursor_mass_tolerance is not None
    assert alpha.precursor_mass_tolerance.unit == "Da"

    wombat_document = {
        "version": "1",
        "params": {
            "enzyme": "lys-c",
            "miscleavages": 1,
            "fixed_mods": "plain",
            "variable_mods": "plain",
            "max_mods": 2,
            "min_peptide_length": 6,
            "max_peptide_length": 30,
            "precursor_mass_tolerance": "10 ppm",
            "fragment_mass_tolerance": "20 ppm",
            "ident_fdr_protein": 0.01,
            "ident_fdr_peptide": 0.01,
            "ident_fdr_psm": 0.01,
            "min_precursor_charge": 1,
            "max_precursor_charge": 4,
            "enable_match_between_runs": False,
            "normalization_method": "none",
        },
    }
    assert wombat.extract_params(StringIO(yaml.safe_dump(wombat_document))).enzyme == "Lys-C"
    wombat_params = cast(dict[str, Any], wombat_document["params"])
    wombat_params["enzyme"] = "trypsin"
    assert wombat.extract_params(StringIO(yaml.safe_dump(wombat_document))).enzyme == "Trypsin"


def test_diann_command_line_and_helper_guards() -> None:
    assert diann._parse_cmdline("", "1.7")["mod"] == []
    with pytest.raises(ValueError, match="invalid `unimod`"):
        diann._parse_cmdline("unimod4 value", "1.7")
    assert diann._parse_cmdline("unimod999", "1.7")["mod"] == []
    below = diann._parse_cmdline("unimod4 --unimod35", "1.7")
    assert below["mod"] == ["Carbamidomethyl (C)"]
    assert below["var-mod"] == ["Oxidation (M)"]
    assert diann._parse_cmdline("mod value", "2.0")["mod"] == ["value"]

    assert diann._extract_cfg(["not-an-int"], r"(.*)", int, default=7) == 7
    assert diann._extract_modifications(["x", "y"], [r"(x\n)"]) == "x"
    assert diann._protein_inference({"no-prot-inf": True}) == "Disabled"
    assert diann._quantification_strategy({"direct-quant": True}) == "Legacy"
    with pytest.raises(ValueError, match="requires a value"):
        diann._protein_inference({"pg-level": []})
    with pytest.raises(TypeError, match="must contain arguments"):
        diann._from_cmdline({"qvalue": "wrong-type"})
    with pytest.raises(TypeError, match="enzyme setting"):
        diann._from_cmdline({"cut": True})
    assert diann._from_cfg(["Normalisation disabled"])["abundance_normalization_ions"] == "None"


def test_diann_rejects_non_numeric_tolerance(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        diann,
        "_defaults",
        lambda: {"fragment_mass_tolerance": object()},
    )
    monkeypatch.setattr(diann, "_from_cmdline", lambda _settings: {})
    monkeypatch.setattr(diann, "_from_log_regex", lambda _lines, *, have: {})
    with pytest.raises(TypeError, match="must be numeric"):
        diann.extract_params(StringIO("diann --flag\n"))


def test_fragpipe_modification_and_workflow_variants() -> None:
    assert fragpipe._parse_fixed_mods("") == ""
    fixed = fragpipe._parse_fixed_mods(
        "bad; "
        "57.02146,C (cysteine),false,-1; "
        "57.02146,N-Term,true,-1; "
        "57.02146,C-Term,true,-1; "
        "12.3,custom,true,-1"
    )
    assert "N-term[Carbamidomethyl]" in fixed
    assert "C-term[Carbamidomethyl]" in fixed
    assert "custom[12.3]" in fixed

    assert fragpipe._parse_variable_mods("") == ""
    variable = fragpipe._parse_variable_mods(
        "bad; 0.0,M,true,1; 42.0106,nM,true,1; 42.0106,n?,true,1"
    )
    assert "N-term M[Acetyl]" in variable
    assert "N-term[Acetyl]" in variable

    parsed = fragpipe._parse_lines(["flag", "key=value # comment"])
    assert parsed[0].value is None
    _header, version, _fragpipe_version, _diann_version, _records = fragpipe._read_workflow(
        "# Header\nfragpipe-config.bin-msfragger=/path/no-version.jar"
    )
    assert version is None

    assert (
        fragpipe._resolve_enzyme(
            pd.Series(
                {
                    "msfragger.search_enzyme_name_1": "enzyme-a",
                    "msfragger.search_enzyme_name_2": "enzyme-b",
                }
            )
        )
        == "enzyme-a|enzyme-b"
    )
    assert (
        fragpipe._resolve_enzyme(
            pd.Series(
                {
                    "msfragger.search_enzyme_name_1": "other",
                    "msfragger.search_enzyme_name_2": "null",
                }
            )
        )
        == "other"
    )
    assert fragpipe._charge_range(
        pd.Series(
            {
                "msfragger.override_charge": "true",
                "msfragger.misc.fragger.precursor-charge-lo": "2",
                "msfragger.misc.fragger.precursor-charge-hi": "5",
            }
        )
    ) == (2, 5)
    assert (
        fragpipe._protein_inference(
            pd.Series(
                {
                    "protein-prophet.run-protein-prophet": "true",
                    "protein-prophet.cmd-opts": "--best",
                }
            )
        )
        == "ProteinProphet: --best"
    )
    assert (
        fragpipe._protein_inference(pd.Series({"protein-prophet.run-protein-prophet": "false"}))
        is None
    )
    neither = pd.Series(
        {
            "diann.run-dia-nn": "false",
            "phi-report.filter": "",
            "quantitation.run-label-free-quant": "false",
        }
    )
    assert fragpipe._fdr_and_mbr(neither)[4:] == (None, None)

    workflow = (Path(__file__).parent / "params" / "fragpipe.workflow").read_text(encoding="utf-8")
    lines = [
        "# Unversioned workflow" if index == 0 else line
        for index, line in enumerate(workflow.splitlines())
        if not line.startswith("# FragPipe version")
    ]
    assert fragpipe.extract_params(StringIO("\n".join(lines))).software_version is None


def test_maxquant_xml_structure_and_scalar_guards(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    assert maxquant._homogenize_mods("") == ""
    data: dict[str, maxquant.XmlValue] = {}
    maxquant._add_record(data, "tag", "one")
    maxquant._add_record(data, "tag", "two")
    maxquant._add_record(data, "tag", "three")
    assert data["tag"] == ["one", "two", "three"]

    with pytest.raises(ValueError, match="root did not parse"):
        maxquant._read_xml(StringIO("<root/>"))
    with pytest.raises(ValueError, match="tuple too long"):
        maxquant._extend(("a", "b"), 1)
    flattened = maxquant._flatten(
        cast(dict[str, maxquant.XmlValue], {"items": [cast(Any, 1), "kept", None]})
    )
    assert flattened == [(("items",), "kept"), (("items",), None)]
    with pytest.raises(TypeError, match="one text"):
        maxquant._text(1, "field")
    with pytest.raises(TypeError, match="entries must be text"):
        maxquant._joined_text(pd.Series(["one", 2]), "field")
    with pytest.raises(TypeError, match="must contain text"):
        maxquant._joined_text(1, "field")

    for record, message in [
        ({}, "must be a list"),
        ({"msmsParamsArray": [1]}, "entries must be mappings"),
        ({"msmsParamsArray": [{}]}, "entry must be a mapping"),
    ]:
        monkeypatch.setattr(maxquant, "_read_xml", lambda _source, value=record: value)
        with pytest.raises(ValueError, match=message):
            maxquant.extract_params(StringIO("<root/>"))


def test_metamorpheus_helpers_reject_invalid_shapes_and_render_termini() -> None:
    with pytest.raises(TypeError, match="must be a table"):
        metamorpheus._mapping([], "field")
    with pytest.raises(TypeError, match="must be text"):
        metamorpheus._text(1, "field")
    assert metamorpheus._homogenize_mod("plain") == "plain"
    assert metamorpheus._homogenize_mod("Acetyl on X (Prot N-Term)") == ("Protein N-term[Acetyl]")
    assert metamorpheus._homogenize_mod("Acetyl on X (Pep N-Term)") == ("N-term[Acetyl]")
    assert metamorpheus._homogenize_mod("Acetyl on X (Prot C-Term)") == ("Protein C-term[Acetyl]")
    assert metamorpheus._homogenize_mod("Acetyl on X (Pep C-Term)") == ("C-term[Acetyl]")
    assert metamorpheus._parse_modifications("ignored") == ""
    with pytest.raises(ValueError, match="expected one TOML"):
        metamorpheus._load_pair(StringIO("first"), StringIO("second"))


def test_peaks_sage_and_spectronaut_error_paths() -> None:
    with pytest.raises(ValueError, match="is missing"):
        peaks._required_value([], "Required:")
    assert peaks._between(
        ["Start", "- first", "End", "Start", "- second", "End"],
        "Start",
        "End",
    ) == ["first", "second"]
    assert peaks._between(
        ["Start", "- unfinished"],
        "Start",
        "End",
        only_last=True,
    ) == ["unfinished"]
    peaks_minimal = "\n".join(
        [
            "PEAKS Version: 1",
            "Peptide Length between: 7,30",
            "Precursor Charge between: 2,4",
            "Precursor M/Z between: 300,1500",
            "Max Missed Cleavage: 2",
            "Max Variable PTM per Peptide: 3",
        ]
    )
    parsed_peaks = peaks.extract_params(StringIO(peaks_minimal))
    assert parsed_peaks.min_precursor_mz == 300
    assert parsed_peaks.min_fragment_mz is None

    sage_document = {
        "database": {
            "enzyme": {
                "cleave_at": "X",
                "semi_enzymatic": "invalid",
            }
        }
    }
    with pytest.raises(ValueError, match="unknown semi_enzymatic"):
        sage.extract_params(StringIO(json.dumps(sage_document)))

    assert spectronaut._homogenize_mods("") == ""
    assert spectronaut._value(["unrelated"], "Missing:") is None
    with pytest.raises(ValueError, match="is missing"):
        spectronaut._required_value([], "Required:")
    assert spectronaut._value_regex(["unrelated"], r"^Missing:") is None
    static_incomplete = [
        "Pulsar Search\\Tolerances",
        "Thermo Orbitrap",
        "Main Search: Static",
        "MS1 Tolerance (Th): 5",
    ]
    assert spectronaut._extract_tolerances(
        static_incomplete,
        "Thermo Orbitrap",
    ) == (None, None)
    reverse_tolerance_order = [
        "Pulsar Search\\Tolerances",
        "Thermo Orbitrap",
        "Main Search: Static",
        "MS2 Tolerance (Th): 6",
        "MS1 Tolerance (Th): 5",
    ]
    assert spectronaut._extract_tolerances(
        reverse_tolerance_order,
        "Thermo Orbitrap",
    ) == ("[-5 Th, 5 Th]", "[-6 Th, 6 Th]")

    with pytest.raises(ValueError, match="unknown Spectronaut vendor"):
        spectronaut.extract_params(StringIO("Version 1\nVendor: Unknown\n"))
    spectronaut_minimal = "\n".join(
        [
            "Spectronaut 1",
            "Vendor: Thermo",
            "Peptide Charge: 3",
            "Missed Cleavages: 2",
            "Max Peptide Length: 30",
            "Min Peptide Length: 7",
            "Max Variable Modifications: 2",
        ]
    )
    parsed_spectronaut = spectronaut.extract_params(StringIO(spectronaut_minimal))
    assert parsed_spectronaut.min_precursor_charge == 3
