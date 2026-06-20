"""WOMBAT-P parameter-file parser (YAML)."""

from __future__ import annotations

from pathlib import Path
from typing import IO, Union

import yaml

from anndata_proteomics.params._common import read_text
from anndata_proteomics.params.model import Parameters


def _homogenize_mod_xtandem(mod_str: str) -> str:
    """Convert a WOMBAT-P X!Tandem modification spec to ProForma-like notation.

    Format: ``{modname} of {residue}``, e.g. ``Oxidation of M``,
    ``Acetyl of Protein N-term``.
    """
    mod_str = mod_str.strip()
    if " of " not in mod_str:
        return mod_str
    name, residue_part = mod_str.split(" of ", 1)
    residue_part = residue_part.strip()
    lower = residue_part.lower()
    if "protein n-term" in lower:
        return f"Protein N-term[{name}]"
    if "n-term" in lower:
        return f"N-term[{name}]"
    if "protein c-term" in lower:
        return f"Protein C-term[{name}]"
    if "c-term" in lower:
        return f"C-term[{name}]"
    return f"{residue_part.upper()}[{name}]"


def extract_params(source: Union[str, Path, IO[bytes], IO[str]]) -> Parameters:
    """Parse a WOMBAT-P YAML configuration into :class:`Parameters`.

    Mirrors ``proteobench.io.params.wombat.extract_params``.
    """
    record = yaml.safe_load(read_text(source))
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
        fixed_mods=", ".join(_homogenize_mod_xtandem(m) for m in p["fixed_mods"].split(",")),
        variable_mods=", ".join(_homogenize_mod_xtandem(m) for m in p["variable_mods"].split(",")),
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
