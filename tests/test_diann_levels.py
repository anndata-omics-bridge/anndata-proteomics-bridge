"""Convert the DIA-NN quantification levels backed by real report layers.

The ``diann_full_subset`` fixture (conftest.py) supplies a small slice of a cached DIA-NN dataset
whose param-parsed version supports ion, protein, and fragment layers (a 1.9.x export:
positional fragment + PG.*), plus that version. Conversion goes through the param-driven
``converters.pipeline`` core. Skips when no such dataset is cached.
"""

from __future__ import annotations

import pytest

from anndata_proteomics.converters import pipeline as ui

_LEVELS = ["ion", "protein", "fragment"]


def _convert(fixture, level):
    return ui._convert_level(fixture["df"].copy(), fixture["slug"], level, fixture["version"])


def test_each_level_converts_to_nonempty_anndata(diann_full_subset) -> None:
    for level in _LEVELS:
        adata = _convert(diann_full_subset, level)
        assert adata.shape[0] == 1, f"{level}: expected one run"
        assert adata.shape[1] > 0, f"{level}: empty var axis"


def test_var_counts_follow_report_backed_levels(diann_full_subset) -> None:
    n = {level: _convert(diann_full_subset, level).shape[1] for level in _LEVELS}
    # many fragments -> ions; protein groups are an independent report-backed level.
    assert n["fragment"] > n["ion"]
    assert n["protein"] > 0


def test_ion_var_preserves_peptide_peptidoform_and_protein_metadata(diann_full_subset) -> None:
    adata = _convert(diann_full_subset, "ion")
    assert {
        "ProForma_peptidoform",
        "ProForma_peptide",
        "Protein_Group",
        "Protein_Ids",
        "Protein_Names",
        "Genes",
    } <= set(adata.var.columns)


def test_fragment_keys_are_positional(diann_full_subset) -> None:
    adata = _convert(diann_full_subset, "fragment")
    sample = list(adata.var_names[:20])
    # ProForma_fragment = "{peptidoform}/{charge}/frag_{i}" -> two slashes + a positional label
    assert all(name.count("/") >= 2 and "/frag_" in name for name in sample)


def test_protein_quant_is_kept_not_summed(diann_full_subset) -> None:
    """PG.MaxLFQ is pre-aggregated and repeated per (Run, Protein.Group); keep_first must
    take that value, not the sum over the precursor rows."""
    sub = diann_full_subset["df"]
    counts = sub.groupby("Protein.Group").size()
    multi_pg = counts[counts > 1].index[0]  # a protein with several precursor rows
    rows = sub[sub["Protein.Group"] == multi_pg]
    raw = float(rows["PG.MaxLFQ"].iloc[0])
    assert len(rows) > 1 and raw > 0

    adata = _convert(diann_full_subset, "protein")
    value = adata[:, multi_pg].layers["PG_MaxLFQ"][0, 0]
    assert value == pytest.approx(raw)  # keep_first, not len(rows) * raw
