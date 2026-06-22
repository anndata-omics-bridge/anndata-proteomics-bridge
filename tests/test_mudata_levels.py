"""MuData proof: wrap DIA-NN report-backed levels on a shared sample axis.

This is the concrete answer to TODO_to_mu_data.md: each quantification level is a normal
AnnData, and MuData(axis=0) is a thin container over them sharing the run (obs) axis. It is
test-only — no public API ships (per the TODO's step 7).

Two design points it validates:
- var_names are prefixed per level so modalities stay globally unique.
- vendor identifiers remain in .var metadata even when they do not create standalone modalities.
"""

from __future__ import annotations

import mudata
from mudata import MuData

from anndata_proteomics.scripts import _ui_support as ui

# Per-level var_names prefix; the bare id stays available in the .var columns.
_PREFIX = {
    "fragment": "frg:",
    "ion": "ion:",
    "protein": "prt:",
}


def _convert(fixture, level):
    return ui._convert_level(fixture["df"].copy(), fixture["slug"], level, fixture["version"])


def _build_levels(fixture) -> dict:
    levels = {}
    for level, prefix in _PREFIX.items():
        adata = _convert(fixture, level)
        adata.var_names = [prefix + str(v) for v in adata.var_names]
        levels[level] = adata
    return levels


def _wire_foreign_keys(levels: dict) -> None:
    """Add prefixed FK columns pointing from each child level into its parent's var_names."""
    frg = levels["fragment"]
    frg.var["ion_fk"] = "ion:" + frg.var["ProForma_ion"].astype(str)


def test_ion_carries_peptide_peptidoform_and_protein_metadata(diann_full_subset) -> None:
    ion = _convert(diann_full_subset, "ion")
    assert {
        "ProForma_peptidoform",
        "ProForma_peptide",
        "Protein_Group",
        "Protein_Ids",
        "Protein_Names",
        "Genes",
    } <= set(ion.var.columns)


def test_foreign_keys_resolve_into_parent_var_names(diann_full_subset) -> None:
    levels = _build_levels(diann_full_subset)
    _wire_foreign_keys(levels)
    ion, prt, frg = levels["ion"], levels["protein"], levels["fragment"]
    # the core chain uses computed keys (never null): fragment -> ion.
    assert set(frg.var["ion_fk"]) <= set(ion.var_names)

    # Ion -> protein-group is carried as vendor metadata; non-blank groups resolve.
    pg = ion.var["Protein_Group"].astype("string")
    real = pg[pg.str.len().fillna(0) > 0].dropna()
    assert set("prt:" + real) <= set(prt.var_names)


def test_mudata_wraps_all_levels_and_round_trips(diann_full_subset, tmp_path) -> None:
    levels = _build_levels(diann_full_subset)
    _wire_foreign_keys(levels)

    mdata = MuData(levels, axis=0)  # axis=0: shared observations (runs)
    assert mdata.n_obs == 1  # the subset is one run
    assert mdata.n_vars == sum(a.n_vars for a in levels.values())

    path = tmp_path / "levels.h5mu"
    mdata.write(path)
    back = mudata.read_h5mu(path)

    assert set(back.mod) == set(_PREFIX)
    # per-level .var survives the round-trip
    for level, adata in levels.items():
        assert back[level].n_vars == adata.n_vars
    # prefixing kept the global axis unique, so the merged .var is not silently emptied
    assert back.var_names.is_unique
    assert back.n_vars == mdata.n_vars
