"""Sage parameter-file parser."""

from __future__ import annotations

import json
from pathlib import Path
from typing import IO, Union

from anndata_proteomics.params._common import format_tolerance_range, lookup_mass_mod, read_text
from anndata_proteomics.params.model import Parameters

# Mass shift (Da) -> human-readable modification name, matched within MASS_TOLERANCE.
MASS_TO_MOD_MAPPING = {
    57.021464: "Carbamidomethyl",
    15.9949: "Oxidation",
    42.0106: "Acetyl",
}
MASS_TOLERANCE = 0.001

# Sage uses "[" for N-terminal and "]" for C-terminal modifications.
RESIDUE_MAP = {"[": "Protein N-term", "]": "Protein C-term", "^": "N-term", "$": "C-term"}


def _lookup_mod_name(mass: float) -> str:
    """Return a modification name for a mass shift within tolerance, else the raw mass."""
    return lookup_mass_mod(mass, MASS_TO_MOD_MAPPING, tol=MASS_TOLERANCE) or str(mass)


def _parse_static_mods(mods: dict) -> str:
    """Render Sage ``static_mods`` ({residue: mass}) as a ProForma-like string."""
    results = []
    for residue, mass in mods.items():
        res = RESIDUE_MAP.get(residue, residue)
        results.append(f"{res}[{_lookup_mod_name(mass)}]")
    return ", ".join(results)


def _parse_variable_mods(mods: dict) -> str:
    """Render Sage ``variable_mods`` ({residue: [masses]}) as a ProForma-like string."""
    results = []
    for residue, masses in mods.items():
        res = RESIDUE_MAP.get(residue, residue)
        for mass in masses:
            results.append(f"{res}[{_lookup_mod_name(mass)}]")
    return ", ".join(results)


def extract_params(source: Union[str, Path, IO[bytes], IO[str]]) -> Parameters:
    """Parse a Sage JSON parameter file into a :class:`Parameters` record.

    Accepts a filesystem path or an open file-like object (bytes or text).
    Field mapping mirrors ProteoBench's ``proteobench.io.params.sage.extract_params``
    so existing expected-output CSVs are reproduced unchanged.
    """
    data = json.loads(read_text(source))

    enzyme = data["database"]["enzyme"]["cleave_at"]
    if enzyme in ("KR", "RK"):
        if "restrict" not in data["database"]["enzyme"]:
            enzyme = "Trypsin/P"
        elif data["database"]["enzyme"]["restrict"] == "P":
            enzyme = "Trypsin"
        # restrict present but not "P" (e.g. null) → keep raw KR/RK

    semi = data["database"]["enzyme"].get("semi_enzymatic")
    if semi is None or semi is False:
        semi_enzymatic = False
    elif semi is True:
        semi_enzymatic = True
    else:
        raise ValueError(f"unknown semi_enzymatic value: {semi!r}")

    max_len = data["database"]["enzyme"]["max_len"]

    return Parameters(
        software_name="Sage",
        software_version=data["version"],
        search_engine="Sage",
        search_engine_version=data["version"],
        enzyme=enzyme,
        semi_enzymatic=semi_enzymatic,
        allowed_miscleavages=data["database"]["enzyme"]["missed_cleavages"],
        fixed_mods=_parse_static_mods(data["database"]["static_mods"]),
        variable_mods=_parse_variable_mods(data["database"]["variable_mods"]),
        precursor_mass_tolerance=format_tolerance_range(data["precursor_tol"]),
        fragment_mass_tolerance=format_tolerance_range(data["fragment_tol"]),
        min_peptide_length=int(data["database"]["enzyme"]["min_len"]),
        max_peptide_length=int(max_len) if max_len is not None else None,
        max_mods=int(data["database"]["max_variable_mods"]),
        min_precursor_charge=int(data["precursor_charge"][0]),
        max_precursor_charge=int(data["precursor_charge"][1]),
        enable_match_between_runs=True,
    )
