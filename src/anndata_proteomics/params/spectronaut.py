"""Spectronaut settings-text parameter parser."""

from __future__ import annotations

import re
from pathlib import Path
from typing import IO, Optional, Union

from anndata_proteomics.params.model import Parameters

_Source = Union[str, Path, IO]

_VENDOR_SYSTEM_MAP = {
    "Thermo": "Thermo Orbitrap",
    "Bruker": "TOF",
}
_MS1_STATIC = re.compile(r"MS1 Tolerance \(Th\):\s*(\d*)")
_MS2_STATIC = re.compile(r"MS2 Tolerance \(Th\):\s*(\d*)")
_MS1_RELATIVE = re.compile(r"MS1 Tolerance \(ppm\):\s*(\d*)")
_MS2_RELATIVE = re.compile(r"MS2 Tolerance \(ppm\):\s*(\d*)")
_MAIN_SEARCH = re.compile(r"Main Search:\s*(.*)")


def _clean(text: str) -> str:
    return re.sub(r"^[\s:,\t]+|[\s:,\t]+$", "", text)


def _value(lines: list[str], term: str) -> Optional[str]:
    for line in lines:
        if term in line:
            return _clean(line.split(term)[1])
    return None


def _value_regex(lines: list[str], pattern: str) -> Optional[str]:
    for line in lines:
        if re.search(pattern, line):
            return _clean(re.split(pattern, line)[1])
    return None


def _extract_tolerances(
    lines: list[str], system: str
) -> tuple[Optional[str], Optional[str]]:
    in_tolerance_block = False
    in_system_block = False
    calibration: Optional[str] = None
    ms1: Optional[str] = None
    ms2: Optional[str] = None

    for line in lines:
        if line.startswith("Pulsar Search\\Tolerances"):
            in_tolerance_block = True
            continue
        if not in_tolerance_block:
            continue
        if line.startswith(system):
            in_system_block = True
            continue
        if not in_system_block:
            continue
        if calibration is None:
            match = _MAIN_SEARCH.search(line)
            if match:
                calibration = match.group(1).strip()
        if calibration == "Dynamic":
            return "Dynamic", "Dynamic"
        if calibration in ("Static", "Relative"):
            unit = "Th" if calibration == "Static" else "ppm"
            ms1_pat, ms2_pat = (
                (_MS1_STATIC, _MS2_STATIC) if calibration == "Static" else (_MS1_RELATIVE, _MS2_RELATIVE)
            )
            if ms1 is None:
                hit = ms1_pat.search(line)
                if hit:
                    ms1 = hit.group(1)
            if ms2 is None:
                hit = ms2_pat.search(line)
                if hit:
                    ms2 = hit.group(1)
            if ms1 is not None and ms2 is not None:
                return f"[-{ms1} {unit}, {ms1} {unit}]", f"[-{ms2} {unit}, {ms2} {unit}]"
    return None, None


def _load_lines(source: _Source) -> list[str]:
    if hasattr(source, "read"):
        raw = source.read()
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8")
        return [line.strip() for line in raw.splitlines()]
    return [line.strip() for line in Path(source).read_text(encoding="utf-8").splitlines()]


def extract_params(source: _Source) -> Parameters:
    """Parse a Spectronaut settings-export text file into :class:`Parameters`.

    Mirrors ``proteobench.io.params.spectronaut.read_spectronaut_settings``.
    """
    lines = _load_lines(source)
    vendor = _value(lines, "Vendor:")
    if vendor not in _VENDOR_SYSTEM_MAP:
        raise ValueError(
            f"unknown Spectronaut vendor: {vendor!r}; expected one of "
            f"{sorted(_VENDOR_SYSTEM_MAP)}"
        )
    system = _VENDOR_SYSTEM_MAP[vendor]

    software_version = lines[0].split()[1]

    # Strip tree-drawing characters present in some Spectronaut exports.
    lines = [re.sub(r"^[\s│├─└]*", "", line).strip() for line in lines]

    precursor_tol, fragment_tol = _extract_tolerances(lines, system)

    psm_raw = _value(lines, "Precursor Qvalue Cutoff:")
    protein_raw = _value(lines, "Protein Qvalue Cutoff (Experiment):")
    ident_psm = float(psm_raw.replace(",", ".")) if psm_raw else None
    ident_protein = float(protein_raw.replace(",", ".")) if protein_raw else None

    charge_raw = _value(lines, "Peptide Charge:")
    if charge_raw is None or charge_raw == "False":
        min_z = max_z = None
    else:
        min_z = max_z = int(charge_raw)

    return Parameters(
        software_name="Spectronaut",
        software_version=software_version,
        search_engine="Spectronaut",
        search_engine_version=software_version,
        ident_fdr_psm=ident_psm,
        ident_fdr_protein=ident_protein,
        enable_match_between_runs=False,
        precursor_mass_tolerance=precursor_tol,
        fragment_mass_tolerance=fragment_tol,
        enzyme=_value(lines, "Enzymes / Cleavage Rules:"),
        semi_enzymatic=_value(lines, "Digest Type:") != "Specific",
        allowed_miscleavages=int(_value(lines, "Missed Cleavages:")),
        max_peptide_length=int(_value(lines, "Max Peptide Length:")),
        min_peptide_length=int(_value(lines, "Min Peptide Length:")),
        fixed_mods=_value(lines, "Fixed Modifications:"),
        variable_mods=_value_regex(lines, r"^Variable Modifications:"),
        max_mods=int(_value(lines, "Max Variable Modifications:")),
        min_precursor_charge=min_z,
        max_precursor_charge=max_z,
        scan_window=_value(lines, "XIC IM Extraction Window:"),
        quantification_method=_value(lines, "Quantity MS Level:"),
        protein_inference=_value(lines, "Inference Algorithm:"),
        abundance_normalization_ions=_value(lines, "Cross-Run Normalization:"),
    )
