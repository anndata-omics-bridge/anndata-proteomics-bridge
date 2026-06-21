"""Convert the five DIA-NN quantification levels from one report.tsv.

DIA-NN's report.tsv backs every level (ion / peptidoform / peptide / protein / fragment).
The ``diann_full_subset`` fixture (conftest.py) supplies a small slice of a cached full-column
DIA-NN file; the fragment level explodes the packed fragment lists, so the slice keeps memory
bounded. Skips when no such file is cached.
"""

from __future__ import annotations

import pytest

from anndata_proteomics.converters.assemble import convert
from anndata_proteomics.rules.loader import load_packaged_rule

_LEVELS = ["ion", "peptidoform", "peptide", "protein", "fragment"]


def test_each_level_converts_to_nonempty_anndata(diann_full_subset) -> None:
    for level in _LEVELS:
        rule = load_packaged_rule("diann", level)
        adata = convert(diann_full_subset.copy(), rule)
        assert adata.shape[0] == 1, f"{level}: expected one run"
        assert adata.shape[1] > 0, f"{level}: empty var axis"
        assert rule.axis.x_layer in adata.layers


def test_var_counts_follow_the_quantification_hierarchy(diann_full_subset) -> None:
    n = {
        level: convert(diann_full_subset.copy(), load_packaged_rule("diann", level)).shape[1]
        for level in _LEVELS
    }
    # many fragments -> ions -> peptidoforms -> peptides -> protein groups
    assert n["fragment"] > n["ion"] > n["peptidoform"] >= n["peptide"] > n["protein"]


def test_fragment_keys_use_the_ion_slash_label_grammar(diann_full_subset) -> None:
    adata = convert(diann_full_subset.copy(), load_packaged_rule("diann", "fragment"))
    sample = list(adata.var_names[:20])
    # ProForma_fragment = "{peptidoform}/{charge}/{fragment_label}" -> at least two slashes
    assert all(name.count("/") >= 2 for name in sample)


def test_protein_quant_is_kept_not_summed(diann_full_subset) -> None:
    """PG.MaxLFQ is pre-aggregated and repeated per (Run, Protein.Group); keep_first must
    take that value, not the sum over the precursor rows."""
    sub = diann_full_subset
    counts = sub.groupby("Protein.Group").size()
    multi_pg = counts[counts > 1].index[0]  # a protein with several precursor rows
    rows = sub[sub["Protein.Group"] == multi_pg]
    raw = float(rows["PG.MaxLFQ"].iloc[0])
    assert len(rows) > 1 and raw > 0

    adata = convert(sub.copy(), load_packaged_rule("diann", "protein"))
    value = adata[:, multi_pg].layers["PG_MaxLFQ"][0, 0]
    assert value == pytest.approx(raw)  # keep_first, not len(rows) * raw
