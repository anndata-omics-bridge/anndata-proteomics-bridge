"""MetaMorpheus parameter-file parser (TOML + version text)."""

from __future__ import annotations

import tomllib
from io import BytesIO
from pathlib import Path
from typing import IO, Union

from anndata_proteomics.params.model import Parameters

_Source = Union[str, Path, IO]


def _format_tolerance(tolerance: str) -> str:
    """Format ``"±20.0000 PPM"`` → ``"[-20.00 PPM, 20.00 PPM]"``."""
    value, unit = tolerance.split()
    value = float(value.strip("±"))
    return f"[-{value:.2f} {unit}, {value:.2f} {unit}]"


def _parse_modifications(mods: str) -> str:
    """Convert MetaMorpheus tab-delimited mod blocks into a ``;``-joined string."""
    parsed: list[str] = []
    for entry in mods.split("\t\t"):
        parts = entry.split("\t")
        if len(parts) > 1:
            parsed.append(parts[1])
    return ";".join(parsed)


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
    if hasattr(source, "read"):
        return _try_load_filelike(source)
    path = Path(source)
    try:
        with open(path, "rb") as fh:
            return tomllib.load(fh)
    except tomllib.TOMLDecodeError:
        return path.read_text(encoding="utf-8").splitlines()[0].strip()


def _try_load_filelike(source: IO):
    try:
        source.seek(0)
    except Exception:
        pass
    try:
        content = source.read()
    finally:
        try:
            source.seek(0)
        except Exception:
            pass
    if isinstance(content, str):
        try:
            return tomllib.load(BytesIO(content.encode("utf-8")))
        except tomllib.TOMLDecodeError:
            return content.splitlines()[0].strip()
    if isinstance(content, bytes):
        try:
            return tomllib.load(BytesIO(content))
        except tomllib.TOMLDecodeError:
            return content.decode("utf-8", errors="replace").splitlines()[0].strip()
    raise TypeError(f"unsupported file-like content type: {type(content).__name__}")


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
