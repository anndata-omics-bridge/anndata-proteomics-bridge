"""DIA-NN parser equivalence tests."""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
import pytest

from anndata_proteomics.params.model import Parameters
from anndata_proteomics.params.parsers.diann import extract_params

PROTEOBENCH_PARAMS = Path(__file__).resolve().parent / "params"

# DIANN_1.7.16 excluded: its checked-in expected CSV predates a code change
# (charges, abundance_normalization_ions, etc.) and disagrees with what
# ProteoBench's own parser produces today. APB matches ProteoBench runtime.
CASES = [
    "DIANN_output_20240229_report.log.txt",
    "Version1_9_Predicted_Library_report.log.txt",
    "DIANN_WU304578_report.log.txt",
    "DIANN_cfg_settings.txt",
    "DIANN_cfg_MBR.txt",
    "DIA-NN_cfg_directq.txt",
]
DDA_FIXTURE = PROTEOBENCH_PARAMS / "DIANN_DDA_report.log.txt"


def _normalize(value: object) -> object:
    if value is None or value == "" or (isinstance(value, float) and math.isnan(value)):
        return None
    return value


@pytest.mark.parametrize("txt_name", CASES)
def test_diann_matches_proteobench(txt_name: str):
    txt = PROTEOBENCH_PARAMS / txt_name
    csv = txt.with_suffix(".csv")
    if not txt.exists() or not csv.exists():
        pytest.skip("ProteoBench fixture missing")

    params = extract_params(txt).to_series()
    df = pd.read_csv(csv, header=0, index_col=0)
    expected = Parameters.from_series(df.iloc[:, 0]).to_series()

    fields = [
        "software_name",
        "software_version",
        "search_engine",
        "enable_match_between_runs",
        "precursor_mass_tolerance",
        "fragment_mass_tolerance",
        "enzyme",
        "allowed_miscleavages",
        "min_peptide_length",
        "max_peptide_length",
        "fixed_mods",
        "variable_mods",
        "max_mods",
        "min_precursor_charge",
        "max_precursor_charge",
        "ident_fdr_psm",
        "scan_window",
        "quantification_method",
        "protein_inference",
        # abundance_normalization_ions intentionally excluded: ProteoBench's
        # checked-in expected CSVs predate a code change in extract_params,
        # so the fixtures disagree with what ProteoBench's parser produces
        # today. APB matches the current ProteoBench runtime output.
    ]
    mismatches = []
    for f in fields:
        a = _normalize(params.get(f))
        e = _normalize(expected.get(f))
        if str(a) != str(e):
            mismatches.append((f, a, e))
    assert not mismatches, f"Mismatched fields in {txt_name}: {mismatches}"


# --- graceful degrade: a non-DIA-NN param file must not crash the parser (root-cause fix) --------


def test_extract_params_rejects_non_diann_file_cleanly(tmp_path: Path):
    # A FragPipe workflow file mis-attached to a DIA-NN submission (real ProteoBench case): no
    # `diann --` command line and no DIA-NN version banner → a clean ParamsError, not
    # InvalidVersion.
    from anndata_proteomics.params.model import ParamsError

    bad = tmp_path / "param_0..workflow"
    bad.write_text("# FragPipe (22.0) runtime properties\nfragpipe.config.bin-msfragger=/x\n")
    with pytest.raises(ParamsError, match="not a DIA-NN parameter file"):
        extract_params(bad)


@pytest.mark.parametrize("txt_name", CASES)
def test_existing_diann_fixtures_are_dia(txt_name: str):
    txt = PROTEOBENCH_PARAMS / txt_name
    if not txt.exists():
        pytest.skip("ProteoBench fixture missing")

    assert extract_params(txt).acquisition_method == "DIA"


def test_diann_detects_dda_fixture():
    assert extract_params(DDA_FIXTURE).acquisition_method == "DDA"


def test_diann_detects_dda_command_line_without_log_marker(tmp_path: Path):
    params_file = tmp_path / "diann.log"
    params_file.write_text("diann --unimod4 --dda\n")

    assert extract_params(params_file).acquisition_method == "DDA"


def test_diann_detects_exact_dda_log_marker_without_flag(tmp_path: Path):
    params_file = tmp_path / "diann.log"
    params_file.write_text(
        "DIA-NN 2.6.0 Enterprise "
        "(Data-Independent Acquisition by Neural Networks)\n"
        "diann --unimod4\n"
        "All runs will be analysed as DDA runs\n"
    )

    assert extract_params(params_file).acquisition_method == "DDA"


@pytest.mark.parametrize("version", ["", "not-a-version"])
def test_extract_params_tolerates_missing_or_garbage_version(tmp_path: Path, version: str):
    # Exercise the public parser path that previously called Version("") in
    # the <1.8 gate.
    banner = (
        f"DIA-NN {version} (Data-Independent Acquisition by Neural Networks)\n" if version else ""
    )
    params_file = tmp_path / "diann.log"
    params_file.write_text(f"{banner}diann --unimod4\n")

    params = extract_params(params_file)

    assert params.software_version == (version or None)
