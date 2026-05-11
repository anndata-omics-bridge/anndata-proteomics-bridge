"""WOMBAT-P parameter-file parser (YAML)."""

from __future__ import annotations

from pathlib import Path
from typing import IO, Union

import yaml

from anndata_proteomics.params.model import Parameters


def extract_params(source: Union[str, Path, IO[bytes], IO[str]]) -> Parameters:
    """Parse a WOMBAT-P YAML configuration into :class:`Parameters`.

    Mirrors ``proteobench.io.params.wombat.extract_params``.
    """
    record = _load_yaml(source)
    p = record["params"]

    enzyme = p["enzyme"]
    if enzyme == "trypsin":
        enzyme = "Trypsin"

    return Parameters(
        software_name="Wombat",
        software_version=record["version"],
        search_engine="various",
        enzyme=enzyme,
        allowed_miscleavages=p["miscleavages"],
        fixed_mods=p["fixed_mods"],
        variable_mods=p["variable_mods"],
        max_mods=p["max_mods"],
        min_peptide_length=p["min_peptide_length"],
        max_peptide_length=p["max_peptide_length"],
        precursor_mass_tolerance=p["precursor_mass_tolerance"],
        fragment_mass_tolerance=p["fragment_mass_tolerance"],
        ident_fdr_protein=p["ident_fdr_protein"],
        ident_fdr_peptide=p["ident_fdr_peptide"],
        ident_fdr_psm=p["ident_fdr_psm"],
        min_precursor_charge=p["min_precursor_charge"],
        max_precursor_charge=p["max_precursor_charge"],
        enable_match_between_runs=p["enable_match_between_runs"],
        abundance_normalization_ions=p["normalization_method"],
    )


def _load_yaml(source: Union[str, Path, IO[bytes], IO[str]]) -> dict:
    if hasattr(source, "read"):
        return yaml.safe_load(source)
    with open(source, "rb") as handle:
        return yaml.safe_load(handle)
