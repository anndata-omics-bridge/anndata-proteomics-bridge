"""Spectronaut TOML semantics: report-backed levels and correct layer placement."""

from __future__ import annotations

import pandas as pd

from anndata_proteomics.converters.assemble import convert
from anndata_proteomics.converters.recognize import matches
from anndata_proteomics.readers.dispatch import read_table
from anndata_proteomics.rules.loader import load_packaged_rule
from anndata_proteomics.rules.registry import resolve_rule_path
from anndata_proteomics.scripts import _ui_support as ui


def _spectronaut_catalog() -> pd.DataFrame:
    catalog = ui.load_catalog()
    return catalog[catalog["slug"] == "spectronaut"].reset_index(drop=True)


def test_spectronaut_has_report_backed_ion_protein_and_fragment_rules() -> None:
    assert resolve_rule_path("spectronaut", "ion") is not None
    assert resolve_rule_path("spectronaut", "protein") is not None
    assert resolve_rule_path("spectronaut", "fragment") is not None
    assert resolve_rule_path("spectronaut", "peptidoform") is None
    assert resolve_rule_path("spectronaut", "peptide") is None


def test_spectronaut_rule_matches_cached_common_headers() -> None:
    rules = [
        load_packaged_rule("spectronaut", "ion"),
        load_packaged_rule("spectronaut", "protein"),
    ]
    for _, row in _spectronaut_catalog().iterrows():
        headers = pd.read_csv(ui._dataset_path(row["input_file_path"]), sep="\t", nrows=0).columns
        for rule in rules:
            assert matches(headers, rule), f"{rule.quantification_level}: {row['input_file_path']}"
        assert not matches(headers, load_packaged_rule("spectronaut", "fragment"))


def test_spectronaut_catalog_offers_mudata() -> None:
    catalog = _spectronaut_catalog()
    if catalog.empty:
        return
    assert catalog["targets"].apply(lambda targets: {"ion", "protein", "mudata"} <= set(targets)).all()
    assert not catalog["targets"].apply(lambda targets: "fragment" in set(targets)).any()


def test_spectronaut_ion_conversion_preserves_metadata_and_layers() -> None:
    catalog = _spectronaut_catalog()
    if catalog.empty:
        return
    row = catalog.iloc[0]
    df = read_table(ui._dataset_path(row["input_file_path"]))
    run = df["R.FileName"].iloc[0]
    subset = df[df["R.FileName"] == run].head(2000).copy()

    adata = convert(subset, load_packaged_rule("spectronaut", "ion"))

    assert {
        "PG_ProteinGroups",
        "PG_ProteinAccessions",
        "ProForma_peptidoform",
        "ProForma_peptide",
        "ProForma_ion",
        "EG_ModifiedSequence",
    } <= set(adata.var.columns)
    assert {
        "FG_Quantity",
        "EG_TargetQuantity_Settings",
        "FG_PrecMz",
    } <= set(adata.layers.keys())
    assert "PG_Quantity" not in adata.layers


def test_spectronaut_protein_conversion_uses_pg_layers() -> None:
    catalog = _spectronaut_catalog()
    if catalog.empty:
        return
    row = catalog.iloc[0]
    df = read_table(ui._dataset_path(row["input_file_path"]))
    run = df["R.FileName"].iloc[0]
    subset = df[df["R.FileName"] == run].head(2000).copy()

    adata = convert(subset, load_packaged_rule("spectronaut", "protein"))

    assert {"PG_ProteinGroups", "PG_ProteinAccessions"} <= set(adata.var.columns)
    assert {
        "PG_Quantity",
        "PG_Cscore",
        "PG_PEP",
        "PG_Qvalue",
        "PG_RunEvidenceCount",
    } <= set(adata.layers.keys())
    assert adata.n_vars == subset["PG.ProteinGroups"].nunique()
