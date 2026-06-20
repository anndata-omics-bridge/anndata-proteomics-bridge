"""MetaMorpheus parameter-file parser (TOML + version text)."""

from __future__ import annotations

import tomllib
from pathlib import Path
from typing import IO, Union

from anndata_proteomics.params._common import read_text
from anndata_proteomics.params.model import Parameters

_Source = Union[str, Path, IO]


def _format_tolerance(tolerance: str) -> str:
    """Format ``"±20.0000 PPM"`` → ``"[-20.00 PPM, 20.00 PPM]"``."""
    value, unit = tolerance.split()
    value = float(value.strip("±"))
    return f"[-{value:.2f} {unit}, {value:.2f} {unit}]"


def _homogenize_mod(mod_str: str) -> str:
    """Convert a MetaMorpheus modification spec to ProForma-like notation.

    MetaMorpheus format: ``{modname} on {residue}`` with optional terminal
    qualifiers like ``(Pep N-Term)`` or ``(Prot N-Term)``.

    Examples:
        ``Carbamidomethyl on C`` -> ``C[Carbamidomethyl]``
        ``Acetylation on X (Prot N-Term)`` -> ``Protein N-term[Acetylation]``
        ``Oxidation on M`` -> ``M[Oxidation]``
    """
    mod_str = mod_str.strip()
    if " on " not in mod_str:
        return mod_str
    name, residue_part = mod_str.split(" on ", 1)
    residue_part = residue_part.strip()
    if "(Prot N-Term)" in residue_part:
        return f"Protein N-term[{name}]"
    if "(Pep N-Term)" in residue_part:
        return f"N-term[{name}]"
    if "(Prot C-Term)" in residue_part:
        return f"Protein C-term[{name}]"
    if "(Pep C-Term)" in residue_part:
        return f"C-term[{name}]"
    return f"{residue_part}[{name}]"


def _parse_modifications(mods: str) -> str:
    """Convert MetaMorpheus tab-delimited mod blocks into a ``, ``-joined string."""
    parsed: list[str] = []
    for entry in mods.split("\t\t"):
        parts = entry.split("\t")
        if len(parts) > 1:
            parsed.append(_homogenize_mod(parts[1]))
    return ", ".join(parsed)


def _load_pair(file_a: _Source, file_b: _Source) -> tuple[str, dict]:
    """Identify which input is the version-text file and which is the TOML."""
    version_line: str | None = None
    settings: dict | None = None

    for source in (file_a, file_b):
        loaded = _try_load(source)
        if isinstance(loaded, dict):
            settings = loaded
        elif isinstance(loaded, str):
            version_line = loaded

    if version_line is None or settings is None:
        raise ValueError("expected one TOML file and one version-text file")
    return version_line, settings


def _try_load(source: _Source):
    """Return a parsed TOML mapping, or the first line of a version-text file."""
    text = read_text(source, errors="replace")
    try:
        return tomllib.loads(text)
    except tomllib.TOMLDecodeError:
        return text.splitlines()[0].strip()


def extract_params(file_a: _Source, file_b: _Source) -> Parameters:
    """Parse a MetaMorpheus TOML + version-text file pair (order-independent).

    Mirrors ``proteobench.io.params.metamorpheus.extract_params``.
    """
    version_line, settings = _load_pair(file_a, file_b)
    common = settings["CommonParameters"]
    search = settings["SearchParameters"]
    digestion = common["DigestionParams"]
    precursor = common["PrecursorDeconvolutionParameters"]

    return Parameters(
        software_name="MetaMorpheus",
        software_version=version_line.split()[2],
        search_engine="MetaMorpheus",
        enzyme=digestion["Protease"],
        allowed_miscleavages=digestion["MaxMissedCleavages"],
        fixed_mods=_parse_modifications(common["ListOfModsFixed"]),
        variable_mods=_parse_modifications(common["ListOfModsVariable"]),
        precursor_mass_tolerance=_format_tolerance(common["PrecursorMassTolerance"]),
        fragment_mass_tolerance=_format_tolerance(common["ProductMassTolerance"]),
        min_peptide_length=digestion["MinPeptideLength"],
        max_peptide_length=digestion["MaxPeptideLength"],
        max_mods=digestion["MaxModsForPeptide"],
        min_precursor_charge=precursor["MinAssumedChargeState"],
        max_precursor_charge=precursor["MaxAssumedChargeState"],
        enable_match_between_runs=bool(search["MatchBetweenRuns"]),
        quantification_method="FlashLFQ",
        protein_inference="Parsimony" if search.get("DoParsimony") else None,
        abundance_normalization_ions=bool(search.get("Normalize")),
        ident_fdr_psm=str(common["QValueThreshold"]),
    )
