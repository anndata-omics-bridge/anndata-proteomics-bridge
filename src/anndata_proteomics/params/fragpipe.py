"""FragPipe ``fragpipe.workflow`` parameter-file parser."""

from __future__ import annotations

import re
from collections import namedtuple
from io import BytesIO
from pathlib import Path
from typing import IO, Union

import pandas as pd

from anndata_proteomics.params._common import lookup_mass_mod, read_text
from anndata_proteomics.params.model import MassTolerance, Parameters

Parameter = namedtuple("Parameter", ["name", "value", "comment"])

_VERSION_NO_PATTERN = r"MSFragger-(.+)\.jar"

_DIANN_QUANT = {
    1: "Any LC (high accuracy)",
    2: "Any LC (high precision)",
    3: "Robust LC (high accuracy)",
    4: "Robust LC (high precision)",
}

# Common mass shifts mapped to modification names (ProForma notation).
_MASS_TO_MOD = {
    57.02146: "Carbamidomethyl",
    15.9949: "Oxidation",
    42.0106: "Acetyl",
    79.96633: "Phospho",
    114.04293: "GG",
    -17.0265: "Pyro-glu",
    -18.0106: "Pyro-glu",
    4.025107: "Label:2H(4)",
    6.020129: "Label:13C(6)",
    8.014199: "Label:13C(6)15N(2)",
    10.008269: "Label:13C(6)15N(4)",
}
_MASS_TOLERANCE = 0.001


def _lookup_mod_name(mass: float) -> str | None:
    """Look up a modification name by mass shift within tolerance."""
    return lookup_mass_mod(mass, _MASS_TO_MOD, tol=_MASS_TOLERANCE)


def _parse_fixed_mods(raw: str) -> str:
    """Parse MSFragger fixed modifications string into ProForma-like format.

    Input format: ``mass,residue_description,active,num_sites`` entries separated by ``; ``.
    Example: ``57.02146,C (cysteine),true,-1``
    """
    if not raw or not raw.strip():
        return ""
    results = []
    for entry in raw.split("; "):
        parts = entry.strip().split(",", 3)
        if len(parts) < 3:
            continue
        mass_str, residue_desc, active = parts[0].strip(), parts[1].strip(), parts[2].strip()
        if active != "true":
            continue
        mass = float(mass_str)
        if abs(mass) < _MASS_TOLERANCE:
            continue
        mod_name = _lookup_mod_name(mass) or mass_str.strip()
        residue_match = re.match(r"^([A-Z])\s*\(", residue_desc)
        if residue_match:
            residue = residue_match.group(1)
        elif "N-Term" in residue_desc:
            residue = "N-term"
        elif "C-Term" in residue_desc:
            residue = "C-term"
        else:
            residue = residue_desc
        results.append(f"{residue}[{mod_name}]")
    return ", ".join(results)


def _parse_variable_mods(raw: str) -> str:
    """Parse MSFragger variable modifications string into ProForma-like format.

    Input format: ``mass,residue,active,max_occurrences`` entries separated by ``; ``.
    Special residue notations: ``[^`` = protein N-term, ``nX`` = peptide N-term of residue X.
    """
    if not raw or not raw.strip():
        return ""
    results = []
    for entry in raw.split("; "):
        parts = entry.strip().split(",", 3)
        if len(parts) < 3:
            continue
        mass_str, residue_field, active = parts[0].strip(), parts[1].strip(), parts[2].strip()
        if active != "true":
            continue
        mass = float(mass_str)
        if abs(mass) < _MASS_TOLERANCE:
            continue
        mod_name = _lookup_mod_name(mass) or mass_str.strip()
        if residue_field == "[^":
            results.append(f"N-term[{mod_name}]")
        elif residue_field.startswith("n"):
            aa_residues = re.findall(r"n([A-Z])", residue_field)
            if aa_residues:
                for aa in aa_residues:
                    results.append(f"N-term {aa}[{mod_name}]")
            else:
                results.append(f"N-term[{mod_name}]")
        else:
            results.append(f"{residue_field}[{mod_name}]")
    return ", ".join(results)


def _parse_lines(lines: list[str], sep: str = "=") -> list[Parameter]:
    """Parse FragPipe ``key=value # comment`` style lines."""
    out: list[Parameter] = []
    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "#" in line:
            parts = line.split("#")
            if len(parts) == 1:
                out.append(Parameter(None, None, parts[0].strip()))
                continue
            param, comment = parts[0].strip(), parts[1].strip()
        else:
            param, comment = line, None
        kv = param.split(sep, maxsplit=1)
        if len(kv) == 1:
            out.append(Parameter(kv[0].strip(), None, comment))
            continue
        out.append(Parameter(kv[0].strip(), kv[1].strip(), comment))
    return out


def _parse_phi_report_filters(cmd: str) -> tuple[float, float, float]:
    """Read PSM/peptide/protein FDR triplet from a ``phi-report.filter`` value."""
    default = 0.01
    patterns = {
        "psm": r"--psm\s+(\d+\.\d+)",
        "peptide": r"--pep\s+(\d+\.\d+)",
        "protein": r"--prot\s+(\d+\.\d+)",
    }
    return tuple(
        float(m.group(1)) if (m := re.search(pat, cmd)) else default
        for pat in (patterns["psm"], patterns["peptide"], patterns["protein"])
    )


def _read_workflow(content: str) -> tuple[str, str | None, str | None, list[Parameter]]:
    lines = content.splitlines()
    header = lines[0][1:].strip()  # leading '#'
    msfragger_version = None
    fragpipe_version = None
    for line in lines[1:]:
        if line.startswith("# MSFragger version"):
            msfragger_version = line.split(" ")[-1].strip()
        elif line.startswith("fragpipe-config.bin-msfragger"):
            path = line.split("=")[-1].strip()
            filename = path.replace("\\", "/").rsplit("/", 1)[-1]
            match = re.search(_VERSION_NO_PATTERN, filename)
            if match:
                msfragger_version = match.group(1)
        if line.startswith("# FragPipe version"):
            fragpipe_version = line.split(" ")[-1].strip()
    return header, msfragger_version, fragpipe_version, _parse_lines(lines)


def extract_params(source: Union[str, Path, IO, BytesIO]) -> Parameters:
    """Parse a FragPipe ``.workflow`` file into :class:`Parameters`.

    Mirrors ``proteobench.io.params.fragger.extract_params``.
    """
    content = read_text(source)
    header, msfragger_version, fragpipe_version, records = _read_workflow(content)
    fp = pd.DataFrame.from_records(records, columns=Parameter._fields).set_index("name")["value"]

    if not fragpipe_version:
        match = re.match(r"FragPipe \((\d+\.\d+.*)\)", header)
        if match:
            fragpipe_version = match.group(1)

    enzyme = fp.loc["msfragger.search_enzyme_name_1"]
    second = fp.loc["msfragger.search_enzyme_name_2"]
    if second != "null":
        enzyme = f"{enzyme}|{second}"
    if enzyme == "stricttrypsin":
        enzyme = "Trypsin/P"
    elif enzyme == "trypsin":
        enzyme = "Trypsin"

    precursor_unit = "ppm" if int(fp.loc["msfragger.precursor_mass_units"]) else "Da"
    precursor_tol = (
        f'[{fp.loc["msfragger.precursor_mass_lower"]} {precursor_unit}, '
        f'{fp.loc["msfragger.precursor_mass_upper"]} {precursor_unit}]'
    )
    fragment_unit = "ppm" if int(fp.loc["msfragger.fragment_mass_units"]) else "Da"
    fragment_tol = MassTolerance(
        mode="absolute",
        value=float(fp.loc["msfragger.fragment_mass_tolerance"]),
        unit=fragment_unit,
    )

    if fp.loc["diann.run-dia-nn"] == "true":
        psm = pep = float(fp.loc["diann.q-value"])
        protein_fdr = float(fp.loc["diann.q-value"])
        peptide_fdr = None
        abundance_norm = None
    else:
        psm, pep, protein_fdr = _parse_phi_report_filters(fp.loc["phi-report.filter"])
        peptide_fdr = pep
        abundance_norm = None

    if fp.loc["msfragger.override_charge"] == "true":
        min_z = int(fp.loc["msfragger.misc.fragger.precursor-charge-lo"])
        max_z = int(fp.loc["msfragger.misc.fragger.precursor-charge-hi"])
    else:
        min_z, max_z = 1, None

    digest_lo = int(fp.loc["msfragger.misc.fragger.digest-mass-lo"])
    digest_hi = int(fp.loc["msfragger.misc.fragger.digest-mass-hi"])
    min_prec_mz = digest_lo / max_z if max_z else None
    max_prec_mz = digest_hi / min_z if min_z else None

    quantification_method = None
    enable_mbr: bool | None = None
    if fp.loc["quantitation.run-label-free-quant"] == "true":
        enable_mbr = bool(int(fp.loc["ionquant.mbr"]))
    elif fp.loc["diann.run-dia-nn"] == "true":
        enable_mbr = (
            ("diann.fragpipe.cmd-opts" in fp.index and "--reanalyse" in fp.loc["diann.fragpipe.cmd-opts"])
            or ("diann.cmd-opts" in fp.index and "--reanalyse" in fp.loc["diann.cmd-opts"])
        )
        quantification_method = _DIANN_QUANT[int(fp.loc["diann.quantification-strategy"])]

    protein_inference = None
    if fp.loc["protein-prophet.run-protein-prophet"] == "true":
        protein_inference = f"ProteinProphet: {fp.loc['protein-prophet.cmd-opts']}"

    return Parameters(
        software_name="FragPipe",
        software_version=fragpipe_version,
        search_engine="MSFragger",
        search_engine_version=msfragger_version,
        enzyme=enzyme,
        allowed_miscleavages=int(fp.loc["msfragger.allowed_missed_cleavage_1"]),
        semi_enzymatic=fp.loc["msfragger.num_enzyme_termini"] != "2",
        fixed_mods=_parse_fixed_mods(fp.loc["msfragger.table.fix-mods"]),
        variable_mods=_parse_variable_mods(fp.loc["msfragger.table.var-mods"]),
        max_mods=int(fp.loc["msfragger.max_variable_mods_per_peptide"]),
        min_peptide_length=int(fp.loc["msfragger.digest_min_length"]),
        max_peptide_length=int(fp.loc["msfragger.digest_max_length"]),
        precursor_mass_tolerance=precursor_tol,
        fragment_mass_tolerance=fragment_tol,
        ident_fdr_psm=psm,
        ident_fdr_peptide=peptide_fdr,
        ident_fdr_protein=protein_fdr,
        enable_match_between_runs=enable_mbr,
        quantification_method=quantification_method,
        protein_inference=protein_inference,
        min_precursor_charge=min_z,
        max_precursor_charge=max_z,
        min_precursor_mz=min_prec_mz,
        max_precursor_mz=max_prec_mz,
        abundance_normalization_ions=abundance_norm,
    )
