"""AlphaPept parameter-file parser (YAML)."""

from __future__ import annotations

from pathlib import Path
from typing import IO

import yaml

from anndata_proteomics.params.model import MassTolerance, Parameters
from anndata_proteomics.params.parsers._common import read_text

MODIFICATION_MAPPING = {
    "cC": "C[Carbamidomethyl]",
    "oxM": "M[Oxidation]",
    "a<^": "N-term[Acetyl]",
}


def _map_modifications(modifications: list[object]) -> str:
    """Map a validated list of AlphaPept modification names."""
    mapped: list[str] = []
    for modification in modifications:
        if not isinstance(modification, str):
            raise TypeError("AlphaPept modification names must be strings")
        name = modification.strip()
        mapped.append(MODIFICATION_MAPPING.get(name, name))
    return ", ".join(mapped)


def extract_params(source: str | Path | IO[bytes] | IO[str]) -> Parameters:
    """Parse an AlphaPept YAML configuration file into :class:`Parameters`.

    Mirrors ``proteobench.io.params.alphapept.extract_params``.
    """
    record = yaml.safe_load(read_text(source))
    summary = record["summary"]
    fasta = record["fasta"]
    search = record["search"]
    features = record["features"]
    workflow = record["workflow"]

    enzyme = fasta["protease"]
    if enzyme == "trypsin":
        enzyme = "Trypsin"

    fixed = (
        list(fasta["mods_fixed"])
        + list(fasta["mods_fixed_terminal"])
        + list(fasta["mods_fixed_terminal_prot"])
    )
    variable = (
        list(fasta["mods_variable"])
        + list(fasta["mods_variable_terminal"])
        + list(fasta["mods_variable_terminal_prot"])
    )

    unit = "ppm" if search["ppm"] else "Da"
    prec_tol = MassTolerance(mode="absolute", value=float(search["prec_tol"]), unit=unit)
    frag_tol = MassTolerance(mode="absolute", value=float(search["frag_tol"]), unit=unit)

    return Parameters.model_validate(
        {
            "software_name": "AlphaPept",
            "software_version": summary["version"],
            "search_engine": "AlphaPept",
            "search_engine_version": summary["version"],
            "enzyme": enzyme,
            "allowed_miscleavages": fasta["n_missed_cleavages"],
            "fixed_mods": _map_modifications(fixed),
            "variable_mods": _map_modifications(variable),
            "max_mods": fasta["n_modifications_max"],
            "min_peptide_length": fasta["pep_length_min"],
            "max_peptide_length": fasta["pep_length_max"],
            "precursor_mass_tolerance": prec_tol,
            "fragment_mass_tolerance": frag_tol,
            "ident_fdr_protein": search["protein_fdr"],
            "ident_fdr_psm": search["peptide_fdr"],
            "min_precursor_charge": features["iso_charge_min"],
            "max_precursor_charge": features["iso_charge_max"],
            "enable_match_between_runs": workflow["match"],
        }
    )
