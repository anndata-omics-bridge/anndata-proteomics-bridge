"""Tests for ``test_data.find_fasta`` and FASTA-vs-module pairing.

Exercises the single-cell module (HY FASTA) alongside an HYE module to
prove the lookup distinguishes the two and that the
``fasta_to_dataframe`` pipeline produces species-distinct frames.
"""

from __future__ import annotations

import pytest

from anndata_proteomics.fasta.annotation import fasta_to_dataframe
from anndata_proteomics.test_data import (
    DOWNLOADED_DB,
    FASTA_DIR,
    TEST_DATA_DIR,
    find_fasta,
    find_test_data,
)

_HYE_NAME = "ProteoBenchFASTA_MixedSpecies_HYE.fasta"
_HY_NAME = "ProteoBenchFASTA_DDAQuantification_noecoli.fasta"


def _require_fasta_cache():
    if not FASTA_DIR.exists() or not any(FASTA_DIR.glob("*.fasta")):
        pytest.skip("FASTA cache not downloaded (run `make -C test_data_download fasta`)")


def test_find_fasta_returns_hye_for_dda_modules():
    _require_fasta_cache()
    fasta = find_fasta(module="dda_qexactive")
    assert fasta is not None
    assert fasta.name == _HYE_NAME


def test_find_fasta_returns_hy_for_singlecell():
    _require_fasta_cache()
    fasta = find_fasta(module="dia_singlecell")
    assert fasta is not None
    assert fasta.name == _HY_NAME


def test_find_fasta_returns_none_for_unknown_module():
    assert find_fasta(module="nonexistent_module") is None


def test_find_fasta_returns_none_without_arguments():
    assert find_fasta() is None


def test_find_fasta_resolves_module_from_dataset_dir():
    _require_fasta_cache()
    if not DOWNLOADED_DB.exists():
        pytest.skip("test_data cache index missing")
    # Use the canonical DIA-NN dataset path (AIF) which we already use in
    # generate_report.py to make sure the dataset-dir branch resolves
    # correctly through the index lookup.
    dataset = find_test_data("DIA-NN")
    if dataset is None:
        pytest.skip("DIA-NN test data not downloaded")
    fasta = find_fasta(dataset_dir=dataset)
    assert fasta is not None
    assert fasta.name == _HYE_NAME


def test_hye_fasta_contains_all_three_species():
    _require_fasta_cache()
    fasta = find_fasta(module="dia_aif")
    df = fasta_to_dataframe(fasta, include_sequence=False)
    suffixes = df["fasta.id"].str.extract(r"_(HUMAN|YEAST|ECOLI)$")[0].dropna()
    counts = suffixes.value_counts().to_dict()
    assert counts.get("HUMAN", 0) > 1000
    assert counts.get("YEAST", 0) > 1000
    assert counts.get("ECOLI", 0) > 1000


def test_hy_fasta_omits_ecoli_proteome():
    _require_fasta_cache()
    fasta = find_fasta(module="dia_singlecell")
    df = fasta_to_dataframe(fasta, include_sequence=False)
    suffixes = df["fasta.id"].str.extract(r"_(HUMAN|YEAST|ECOLI)$")[0].dropna()
    counts = suffixes.value_counts().to_dict()
    assert counts.get("HUMAN", 0) > 1000
    assert counts.get("YEAST", 0) > 1000
    # HY: the only ECOLI entry is the curated BGAL contaminant
    # (`Cont_P00722|BGAL_ECOLI`). The single-cell FASTA must not carry
    # the full E. coli proteome.
    assert counts.get("ECOLI", 0) <= 5


def test_singlecell_diann_test_data_is_available():
    """Companion to the HY-fasta tests: the input file we'd join against."""
    if not DOWNLOADED_DB.exists():
        pytest.skip("test_data cache index missing")
    import csv

    found = False
    with open(DOWNLOADED_DB) as f:
        for row in csv.DictReader(f):
            if (
                row.get("status") == "ok"
                and row.get("module") == "dia_singlecell"
                and row.get("software_name") == "DIA-NN"
            ):
                found = True
                dataset = TEST_DATA_DIR / "json_dir" / row["input_file_path"]
                assert dataset.exists()
                break
    if not found:
        pytest.skip("DIA-NN single-cell input not downloaded")
