"""FragPipe ``fragpipe.workflow`` parameter-file parser."""

from __future__ import annotations

import re
from collections import namedtuple
from io import BytesIO
from pathlib import Path
from typing import IO, Union

import pandas as pd

from anndata_proteomics.params.model import Parameters

Parameter = namedtuple("Parameter", ["name", "value", "comment"])

_VERSION_NO_PATTERN = r"MSFragger-(.+)\.jar"

_DIANN_QUANT = {
    1: "An" "y LC (high accuracy)",
    2: "An" "y LC (high precision)",
    3: "Robust LC (high accuracy)",
    4: "Robust LC (high precision)",
}


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


def _load_text(source: Union[str, Path, IO]) -> str:
    if hasattr(source, "read"):
        try:
            source.seek(0)
        except Exception:
            pass
        raw = source.read()
        if isinstance(raw, bytes):
            return raw.decode("utf-8")
        return raw
    return Path(source).read_text(encoding="utf-8")


def extract_params(source: Union[str, Path, IO, BytesIO]) -> Parameters:
    """Parse a FragPipe ``.workflow`` file into :class:`Parameters`.

    Mirrors ``proteobench.io.params.fragger.extract_params``.
    """
    content = _load_text(source)
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
    fragment_tol_value = fp.loc["msfragger.fragment_mass_tolerance"]
    fragment_tol = f"[-{fragment_tol_value} {fragment_unit}, {fragment_tol_value} {fragment_unit}]"

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
        fixed_mods=fp.loc["msfragger.table.fix-mods"],
        variable_mods=fp.loc["msfragger.table.var-mods"],
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
