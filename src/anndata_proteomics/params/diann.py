"""DIA-NN log/cfg parameter-file parser."""

from __future__ import annotations

import re
from collections.abc import Callable
from pathlib import Path
from typing import IO, Optional, Union

from packaging.version import Version

from anndata_proteomics.params._common import read_lines
from anndata_proteomics.params.model import MassTolerance, Parameters

_Source = Union[str, Path, IO]

MODIFICATION_MAPPING = {
    # Command-line short forms
    "unimod4": "C[Carbamidomethyl]",
    # Descriptive forms
    "Carbamidomethyl (C)": "C[Carbamidomethyl]",
    "Cysteine carbamidomethylation": "C[Carbamidomethyl]",
    "Oxidation (M)": "M[Oxidation]",
    "Acetyl": "N-term[Acetyl]",
    # UniMod short forms (from cfg-extracted log text)
    "UniMod:4": "C[Carbamidomethyl]",
    "UniMod:35": "M[Oxidation]",
    "UniMod:1": "N-term[Acetyl]",
    "UniMod:21": "S[Phospho], T[Phospho], Y[Phospho]",
    "UniMod:121": "K[GG]",
    # UniMod full forms with slash separators (from command-line parsing)
    "UniMod:35/15.994915/M": "M[Oxidation]",
    "UniMod:1/42.010565/*n": "N-term[Acetyl]",
    "UniMod:21/79.966331/STY": "STY[Phospho]",
    "UniMod:121/114.042927/K": "K[GG]",
    # UniMod full forms with comma separators (alternative notation)
    "UniMod:1,42.010565,*n": "N-term[Acetyl]",
    "UniMod:21,79.966331,STY": "STY[Phospho]",
    "UniMod:121,114.042927,K": "K[GG]",
}

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
_SettingValue = bool | int | float | str | list[str] | dict[str, str]

# DIA-NN built-in defaults, reported when the log/cfg omits the corresponding
# setting. These mirror DIA-NN's own built-in defaults and are version-sensitive:
# re-verify against DIA-NN release notes when bumping supported versions.
_DIANN_IMPLICIT_DEFAULTS: dict[str, object] = {
    "min_precursor_charge": 1,
    "max_precursor_charge": 4,
    "min_peptide_length": 7,
    "max_peptide_length": 30,
    "min_fragment_mz": 200,
    "max_fragment_mz": 1800,
    "min_precursor_mz": 300,
    "max_precursor_mz": 1800,
}


def _find_cmdline(lines: list[str]) -> Optional[str]:
    for line in lines:
        if "diann" in line and "--" in line:
            return line.strip()
    return None


def _parse_cmdline(cmd: str, software_version: str) -> dict[str, _SettingValue]:
    settings: dict[str, _SettingValue] = {}
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


def _coerce(setting_name: str, values: list[str]) -> float | int | str:
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


def _extract_cfg(
    lines: list[str],
    regex: str,
    cast_type: Callable[[str], object] = str,
    default: object = None,
    search_all: bool = False,
) -> object:
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


def _protein_inference(cmd_dict: dict[str, _SettingValue]) -> str:
    if "no-prot-inf" in cmd_dict:
        return "Disabled"
    if "pg-level" in cmd_dict:
        pg = cmd_dict["pg-level"][0]
        return {"0": "Isoforms", "1": "Protein_names", "2": "Genes"}.get(pg, "Genes")
    return "Genes"


def _quantification_strategy(cmd_dict: dict[str, _SettingValue]) -> str:
    if "direct-quant" in cmd_dict:
        return "Legacy"
    if "high-acc" in cmd_dict:
        return "QuantUMS high-accuracy"
    return "QuantUMS high-precision"


def _predictors_library(cmd_dict: dict[str, _SettingValue]) -> str | None:
    if "predictor" in cmd_dict:
        return "{'RT': 'DIANN', 'IM': 'DIANN', 'MS2_int': 'DIANN'}"
    if "lib" in cmd_dict and not isinstance(cmd_dict["lib"], bool):
        return (
            "{'RT': 'User defined speclib', 'IM': 'User defined speclib', "
            "'MS2_int': 'User defined speclib'}"
        )
    return None


def _normalize_enzyme(enzyme_str: str) -> str:
    if enzyme_str == "K*,R*":
        return "Trypsin/P"
    if enzyme_str == "K*,R*,!P":
        return "Trypsin"
    return enzyme_str


def _defaults() -> dict[str, object]:
    """Static defaults seeded before any log/cmdline/cfg parsing."""
    return {
        "software_name": "DIA-NN",
        "search_engine": "DIA-NN",
        "enable_match_between_runs": False,
        "quantification_method": "QuantUMS high-precision",
        "protein_inference": "Genes",
        **_DIANN_IMPLICIT_DEFAULTS,
    }


def _from_cmdline(cmd_dict: dict[str, _SettingValue]) -> dict[str, object]:
    """Settings derived from the ``diann --...`` command line.

    Tolerance fields land here as raw numeric values; they are normalized to
    typed :class:`MassTolerance` once in :func:`extract_params`.
    """
    out: dict[str, object] = {
        "quantification_method": _quantification_strategy(cmd_dict),
        "protein_inference": _protein_inference(cmd_dict),
        "predictors_library": _predictors_library(cmd_dict),
    }
    for pb_name, cmd_name in _PARAM_CMD_DICT.items():
        if cmd_name not in cmd_dict:
            continue
        if isinstance(cmd_dict[cmd_name], bool):
            out[pb_name] = cmd_dict[cmd_name]
        else:
            out[pb_name] = _coerce(pb_name, cmd_dict[cmd_name])

    enzyme = out.get("enzyme")
    # Missing enzyme happens when running fragpipe-diann or kept as GUI default.
    out["enzyme"] = "Trypsin/P" if enzyme is None else _normalize_enzyme(enzyme)
    out["abundance_normalization_ions"] = (
        "None" if "no-norm" in cmd_dict else "Cross-run normalization"
    )
    return out


def _from_log_regex(lines: list[str], have: set[str]) -> dict[str, object]:
    """In-log regex fallbacks: tolerances gap-fill the command line, scan window overrides it."""
    out: dict[str, object] = {}
    if "fragment_mass_tolerance" not in have:
        out["fragment_mass_tolerance"] = _extract_with_regex(lines, _FRAGMENT_TOL)
    if "precursor_mass_tolerance" not in have:
        out["precursor_mass_tolerance"] = _extract_with_regex(lines, _PRECURSOR_TOL)
    scan_window = _extract_with_regex(lines, _SCAN_WINDOW)
    out["scan_window"] = int(scan_window) if scan_window is not None else None
    return out


def _from_cfg(lines: list[str]) -> dict[str, object]:
    """Settings re-read from the ``--cfg`` free-text block when a config file was used."""
    out: dict[str, object] = {
        "ident_fdr_psm": _extract_cfg(lines, _FDR, float),
        "ident_fdr_protein": None,
        "enable_match_between_runs": bool(re.search(_MBR_FLAG, "".join(lines))),
        "enzyme": _normalize_enzyme(
            f"{_extract_cfg(lines, _CLEAVAGE) or ''},"
            f"!{(_extract_cfg(lines, _CLEAVAGE_EXC) or '').strip('*')}"
        ),
        "allowed_miscleavages": _extract_cfg(lines, _MISSED_CLEAVAGES, int),
        "min_peptide_length": _extract_cfg(lines, _MIN_PEP_LEN, int),
        "max_peptide_length": _extract_cfg(lines, _MAX_PEP_LEN, int),
        "min_precursor_charge": _extract_cfg(lines, _MIN_Z, int),
        "max_precursor_charge": _extract_cfg(lines, _MAX_Z, int),
        "max_mods": _extract_cfg(lines, _MAX_MODS, int, default=0),
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
    if re.search(_NORMALISATION_DISABLED, "".join(lines)):
        out["abundance_normalization_ions"] = "None"
    inference = _extract_cfg(lines, _PROTEIN_INFERENCE)
    out["protein_inference"] = _PROT_INF_MAP.get(inference, "Genes")
    return out


def extract_params(source: _Source) -> Parameters:
    """Parse a DIA-NN log file into :class:`Parameters`.

    Mirrors ``proteobench.io.params.diann.extract_params``. Walks the log,
    finds the ``diann --...`` command line, applies command-line settings
    via ``_PARAM_CMD_DICT``, falls back to in-log regex extraction for
    fragment/precursor tolerances and the scan window, and finally
    re-reads from the ``--cfg`` free-text block when a config file was
    used. Stages merge with explicit precedence:
    defaults < command line < log regex (gap-fill) < cfg block.
    """
    lines = read_lines(source)
    software_version = _extract_with_regex(lines, _SOFTWARE_VERSION)
    cmdline = _find_cmdline(lines) or ""
    cfg_used = bool(cmdline) and "--cfg" in cmdline
    cmd_dict = _parse_cmdline(cmdline, software_version or "")

    out = _defaults()
    out["software_version"] = software_version
    out["search_engine_version"] = software_version
    out.update(_from_cmdline(cmd_dict))
    out.update(_from_log_regex(lines, have=set(out)))
    if cfg_used:
        out.update(_from_cfg(lines))

    # Normalize tolerances to typed MassTolerance once. DIA-NN tolerances are
    # always a symmetric ppm half-width; the value comes from either the command
    # line (numeric) or the in-log regex (string), so coerce via float().
    for key in ("fragment_mass_tolerance", "precursor_mass_tolerance"):
        value = out.get(key)
        if value in (None, ""):
            continue
        out[key] = MassTolerance(mode="absolute", value=float(value), unit="ppm")

    # Map modification strings to ProForma-like notation.
    for mod_key in ("fixed_mods", "variable_mods"):
        raw = out.get(mod_key)
        if not isinstance(raw, str) or not raw:
            continue
        mapped = [MODIFICATION_MAPPING.get(mod.strip(), mod.strip()) for mod in raw.split(",")]
        out[mod_key] = ", ".join(mapped)

    return Parameters(**out)
