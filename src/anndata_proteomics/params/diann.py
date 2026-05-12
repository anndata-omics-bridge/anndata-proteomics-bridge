"""DIA-NN log/cfg parameter-file parser."""

from __future__ import annotations

import re
from pathlib import Path
from typing import IO, Any, Optional, Union

from packaging.version import Version

from anndata_proteomics.params.model import Parameters

_Source = Union[str, Path, IO]

_FRAGMENT_TOL = r"Optimised mass accuracy: (\d*\.?\d+) ppm"
_PRECURSOR_TOL = r"Recommended MS1 mass accuracy setting: (\d*\.?\d+) ppm"
_SOFTWARE_VERSION = r"DIA-NN\s(.*?)\s\(Data-Independent Acquisition by Neural Networks\)"
_SCAN_WINDOW = r"Scan window radius set to (\d+)"
_FDR = r"Output will be filtered at (\d+\.\d+) FDR"
_MIN_PEP_LEN = r"Min peptide length set to (\d+)"
_MAX_PEP_LEN = r"Max peptide length set to (\d+)"
_MIN_Z = r"Min precursor charge set to (\d+)"
_MAX_Z = r"Max precursor charge set to (\d+)"
_MIN_MZ_PREC = r"Min precursor m/z set to (\d+)"
_MAX_MZ_PREC = r"Max precursor m/z set to (\d+)"
_MIN_MZ_FRAG = r"Min fragment m/z set to (\d+)"
_MAX_MZ_FRAG = r"Max fragment m/z set to (\d+)"
_CLEAVAGE = r"In silico digest will involve cuts at (.*)"
_CLEAVAGE_EXC = r"But excluding cuts at (.*)"
_MISSED_CLEAVAGES = r"Maximum number of missed cleavages set to (\d+)"
_MAX_MODS = r"Maximum number of variable modifications set to (\d+)"
_FIXED_MODS_1 = r"(.*) enabled as a fixed modification"
_FIXED_MODS_2 = r"Modification (.*) with mass delta \d+\.*\d* at .+ will be considered as fixed"
_VAR_MODS = r"Modification (.*) with mass delta \d+\.*\d* at .+ will be considered as variable"
_QUANT_MODE = r"(.*?) quantification mode"
_PROTEIN_INFERENCE = r"Implicit protein grouping: (.*);"
_NORMALISATION_DISABLED = r"(Normalisation disabled)"
_MBR_FLAG = r"(MBR enabled)|(reanalyse them)"

_PARAM_CMD_DICT = {
    "ident_fdr_psm": "qvalue",
    "enable_match_between_runs": "reanalyse",
    "precursor_mass_tolerance": "mass-acc-ms1",
    "fragment_mass_tolerance": "mass-acc",
    "enzyme": "cut",
    "allowed_miscleavages": "missed-cleavages",
    "min_peptide_length": "min-pep-len",
    "max_peptide_length": "max-pep-len",
    "min_fragment_mz": "min-fr-mz",
    "max_fragment_mz": "max-fr-mz",
    "min_precursor_mz": "min-pr-mz",
    "max_precursor_mz": "max-pr-mz",
    "fixed_mods": "mod",
    "variable_mods": "var-mod",
    "max_mods": "var-mods",
    "min_precursor_charge": "min-pr-charge",
    "max_precursor_charge": "max-pr-charge",
    "scan_window": "window",
    "protein_inference": "pg-level",
}
_SETTINGS_PB_FLOAT = {
    "ident_fdr_psm",
    "ident_fdr_peptide",
    "ident_fdr_protein",
    "precursor_mass_tolerance",
    "fragment_mass_tolerance",
}
_SETTINGS_PB_INT = {
    "allowed_miscleavages",
    "min_peptide_length",
    "max_peptide_length",
    "max_mods",
    "min_precursor_charge",
    "max_precursor_charge",
    "scan_window",
}
_SETTINGS_PB_MOD = {"fixed_mods", "variable_mods"}

_PROT_INF_MAP = {"isoform IDs": "Isoforms", "protein names": "Protein_names", "genes": "Genes"}


def _find_cmdline(lines: list[str]) -> Optional[str]:
    for line in lines:
        if "diann" in line and "--" in line:
            return line.strip()
    return None


def _parse_cmdline(cmd: str, software_version: str) -> dict:
    settings: dict[str, Any] = {}
    var_mods: list[str] = []
    fixed_mods: list[str] = []
    below_1_8 = Version(software_version.split(" ")[0]) < Version("1.8")

    for parts in (s.split() for s in cmd.split(" --")):
        key, values = parts[0], parts[1:]
        if key.startswith("unimod"):
            if len(parts) != 1:
                raise ValueError(f"invalid `unimod` format: {parts}")
            if below_1_8:
                if key == "unimod4":
                    fixed_mods.append("Carbamidomethyl (C)")
                elif key == "unimod35":
                    var_mods.append("Oxidation (M)")
            else:
                fixed_mods.append(key)
        elif len(parts) == 1:
            settings[key] = True
        elif key == "var-mod":
            var_mods.append("".join(values).replace(",", "/"))
        else:
            settings[key] = values

    settings["var-mod"] = var_mods
    if "mod" not in settings:
        settings["mod"] = fixed_mods
    return settings


def _coerce(setting_name: str, values: list[str]):
    if setting_name in _SETTINGS_PB_FLOAT:
        return float(values[0])
    if setting_name in _SETTINGS_PB_INT:
        return int(values[0])
    if setting_name in _SETTINGS_PB_MOD:
        return ",".join(values)
    return "".join(values)


def _extract_with_regex(lines: list[str], regex: str, search_all: bool = False) -> Optional[str]:
    container: list[str] = []
    for line in lines:
        match = re.search(regex, line)
        if not match:
            continue
        if not search_all:
            return match.group(1)
        container.append(match.group(1))
    return container[-1] if container else None


def _extract_cfg(lines: list[str], regex: str, cast_type=str, default=None, search_all: bool = False):
    raw = _extract_with_regex(lines, regex, search_all=search_all)
    if raw is None:
        return default
    try:
        return cast_type(raw)
    except ValueError:
        return default


def _extract_modifications(lines: list[str], regexes: list[str]) -> Optional[str]:
    joined = "\n".join(lines)
    mods: list[str] = []
    for regex in regexes:
        for match in re.finditer(regex, joined):
            value = match.group(1)
            if not value.endswith("\n"):
                value = value + "\n"
            mods.append(value)
    return ",".join(mods).replace("\n", "") if mods else None


def _protein_inference(cmd_dict: dict) -> str:
    if "no-prot-inf" in cmd_dict:
        return "Disabled"
    if "pg-level" in cmd_dict:
        pg = cmd_dict["pg-level"][0]
        return {"0": "Isoforms", "1": "Protein_names", "2": "Genes"}.get(pg, "Genes")
    return "Genes"


def _quantification_strategy(cmd_dict: dict) -> str:
    if "direct-quant" in cmd_dict:
        return "Legacy"
    if "high-acc" in cmd_dict:
        return "QuantUMS high-accuracy"
    return "QuantUMS high-precision"


def _predictors_library(cmd_dict: dict):
    if "predictor" in cmd_dict:
        return {"RT": "DIANN", "IM": "DIANN", "MS2_int": "DIANN"}
    if "lib" in cmd_dict and not isinstance(cmd_dict["lib"], bool):
        return {"RT": "User defined speclib", "IM": "User defined speclib", "MS2_int": "User defined speclib"}
    return None


def _load_lines(source: _Source) -> list[str]:
    if hasattr(source, "read"):
        raw = source.read()
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8")
        return raw.splitlines()
    return Path(source).read_text(encoding="utf-8").splitlines()


def extract_params(source: _Source) -> Parameters:
    """Parse a DIA-NN log file into :class:`Parameters`.

    Mirrors ``proteobench.io.params.diann.extract_params``. Walks the log,
    finds the ``diann --...`` command line, applies command-line settings
    via ``_PARAM_CMD_DICT``, falls back to in-log regex extraction for
    fragment/precursor tolerances and the scan window, and finally
    re-reads from the ``--cfg`` free-text block when a config file was
    used.
    """
    lines = _load_lines(source)
    out: dict[str, Any] = {
        "software_name": "DIA-NN",
        "search_engine": "DIA-NN",
        "enable_match_between_runs": False,
        "quantification_method": "QuantUMS high-precision",
        "protein_inference": "Genes",
        "min_precursor_charge": 1,
        "max_precursor_charge": 4,
        "min_peptide_length": 7,
        "max_peptide_length": 30,
        "min_fragment_mz": 200,
        "max_fragment_mz": 1800,
        "min_precursor_mz": 300,
        "max_precursor_mz": 1800,
    }

    software_version = _extract_with_regex(lines, _SOFTWARE_VERSION)
    out["software_version"] = software_version
    out["search_engine_version"] = software_version

    cfg_used = False
    cmdline = _find_cmdline(lines) or ""
    if cmdline and "--cfg" in cmdline:
        cfg_used = True
    cmd_dict = _parse_cmdline(cmdline, software_version or "")

    out["quantification_method"] = _quantification_strategy(cmd_dict)
    out["protein_inference"] = _protein_inference(cmd_dict)
    out["predictors_library"] = _predictors_library(cmd_dict)

    for pb_name, cmd_name in _PARAM_CMD_DICT.items():
        if cmd_name not in cmd_dict:
            continue
        if isinstance(cmd_dict[cmd_name], bool):
            out[pb_name] = cmd_dict[cmd_name]
        else:
            out[pb_name] = _coerce(pb_name, cmd_dict[cmd_name])

    enzyme = out.get("enzyme")
    if enzyme is None:
        out["enzyme"] = "cut"
    elif enzyme == "K*,R*":
        out["enzyme"] = "Trypsin/P"
    elif enzyme == "K*,R*,!*P":
        out["enzyme"] = "Trypsin"

    if "fragment_mass_tolerance" not in out:
        frag = _extract_with_regex(lines, _FRAGMENT_TOL)
        out["fragment_mass_tolerance"] = f"[-{frag} ppm, {frag} ppm]"
    else:
        v = out["fragment_mass_tolerance"]
        out["fragment_mass_tolerance"] = f"[-{v} ppm, {v} ppm]"

    if "precursor_mass_tolerance" not in out:
        prec = _extract_with_regex(lines, _PRECURSOR_TOL)
        out["precursor_mass_tolerance"] = f"[-{prec} ppm, {prec} ppm]"
    else:
        v = out["precursor_mass_tolerance"]
        out["precursor_mass_tolerance"] = f"[-{v} ppm, {v} ppm]"

    scan_window = _extract_with_regex(lines, _SCAN_WINDOW)
    out["scan_window"] = int(scan_window) if scan_window is not None else None
    out["abundance_normalization_ions"] = "None" if "no-norm" in cmd_dict else "Cross-run normalization"

    if cfg_used:
        out.update(
            {
                "ident_fdr_psm": _extract_cfg(lines, _FDR, float),
                "ident_fdr_protein": None,
                "enable_match_between_runs": bool(re.search(_MBR_FLAG, "".join(lines))),
                "enzyme": (
                    f"{_extract_cfg(lines, _CLEAVAGE) or ''},"
                    f"!{_extract_cfg(lines, _CLEAVAGE_EXC) or ''}"
                ),
                "allowed_miscleavages": _extract_cfg(lines, _MISSED_CLEAVAGES, int),
                "min_peptide_length": _extract_cfg(lines, _MIN_PEP_LEN, int),
                "max_peptide_length": _extract_cfg(lines, _MAX_PEP_LEN, int),
                "min_precursor_charge": _extract_cfg(lines, _MIN_Z, int),
                "max_precursor_charge": _extract_cfg(lines, _MAX_Z, int),
                "max_mods": _extract_cfg(lines, _MAX_MODS, int),
                "quantification_method": _extract_cfg(
                    lines, _QUANT_MODE, str, "QuantUMS high-precision", search_all=True
                ),
                "fixed_mods": _extract_modifications(lines, [_FIXED_MODS_1, _FIXED_MODS_2]),
                "variable_mods": _extract_modifications(lines, [_VAR_MODS]),
                "min_fragment_mz": _extract_cfg(lines, _MIN_MZ_FRAG, int),
                "max_fragment_mz": _extract_cfg(lines, _MAX_MZ_FRAG, int),
                "min_precursor_mz": _extract_cfg(lines, _MIN_MZ_PREC, int),
                "max_precursor_mz": _extract_cfg(lines, _MAX_MZ_PREC, int),
            }
        )
        if re.search(_NORMALISATION_DISABLED, "".join(lines)):
            out["abundance_normalization_ions"] = "None"
        inference = _extract_cfg(lines, _PROTEIN_INFERENCE)
        out["protein_inference"] = _PROT_INF_MAP.get(inference, "Genes")

    return Parameters(**out)
