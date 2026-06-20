"""MaxQuant ``mqpar.xml`` parameter-file parser."""

from __future__ import annotations

import collections.abc
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import IO, Union

import pandas as pd

from anndata_proteomics.params._common import homogenize_paren_mods
from anndata_proteomics.params.model import MassTolerance, Parameters

XmlValue = str | dict[str, "XmlValue"] | list["XmlValue"] | None
FlatValue = str | None

# Fallback mapping for modifications without parenthesized residue specifiers.
_MODIFICATION_MAPPING = {
    "Cys-Cys": "C[Disulfide]",
    "Cysteinyl": "C[Cysteinyl]",
    "Cysteinyl - carbamidomethyl": "C[Cysteinyl + Carbamidomethyl]",
}


def _homogenize_mods(raw_mods: str, sep: str = ",") -> str:
    """Parse and homogenize a separator-delimited ``{name} ({residues})`` string."""
    if not raw_mods or not raw_mods.strip():
        return ""
    return ", ".join(
        homogenize_paren_mods(mod, _MODIFICATION_MAPPING)
        for mod in raw_mods.split(sep)
        if mod.strip()
    )


def _add_record(data: dict[str, XmlValue], tag: str, record: XmlValue) -> dict[str, XmlValue]:
    if tag in data:
        if isinstance(data[tag], list):
            data[tag].append(record)
        else:
            data[tag] = [data[tag], record]
    else:
        data[tag] = record
    return data


def _read_element(element: ET.Element) -> XmlValue:
    data: dict[str, XmlValue] = {}
    if element.attrib:
        data.update(element.attrib)
    for child in element:
        if len(child) > 1 and child.tag:
            # Each list item wraps grandchild as {grandchild.tag: parsed-value}.
            data[child.tag] = [
                _add_record(
                    {},
                    tag=grand.tag,
                    record=(
                        grand.text.strip()
                        if (grand.text and grand.text.strip())
                        else _read_element(grand)
                    ),
                )
                for grand in child
            ]
        elif child.text and child.text.strip():
            _add_record(data, child.tag, child.text.strip())
        else:
            _add_record(data, child.tag, _read_element(child))
    return data or None


def _read_xml(source: Union[str, Path, IO[bytes], IO[str]]) -> dict[str, XmlValue]:
    tree = ET.parse(source)
    parsed = _read_element(tree.getroot())
    if not isinstance(parsed, dict):
        raise ValueError("mqpar root did not parse to a mapping")
    return parsed


def _extend(t: tuple, target_length: int) -> tuple:
    if len(t) > target_length:
        raise ValueError(f"tuple too long for index width {target_length}: {t!r}")
    return t + (None,) * (target_length - len(t))


def _flatten(
    d: dict[str, XmlValue], parent_key: tuple = ()
) -> list[tuple[tuple, FlatValue]]:
    items: list[tuple[tuple, FlatValue]] = []
    for key, value in d.items():
        new_key = parent_key + (key,)
        if isinstance(value, collections.abc.MutableMapping):
            items.extend(_flatten(value, parent_key=new_key))
        elif isinstance(value, list):
            for item in value:
                if isinstance(item, collections.abc.MutableMapping):
                    items.extend(_flatten(item, parent_key=new_key))
                elif isinstance(item, str) or item is None:
                    items.append((new_key, item))
        else:
            items.append((new_key, value))
    return items


def _build_series(record: dict, index_length: int = 4) -> pd.Series:
    items = _flatten(record)
    idx = pd.MultiIndex.from_tuples(_extend(k, index_length) for (k, _) in items)
    return pd.Series((v for (_, v) in items), index=idx)


def _tolerance_pair(series: pd.Series) -> tuple[MassTolerance, MassTolerance]:
    """Build precursor (ppm) and fragment (ppm/Da) tolerances from the mqpar series."""
    prec_value = float(series.loc[pd.IndexSlice["parameterGroups", "parameterGroup", "mainSearchTol", :]].squeeze())
    precursor = MassTolerance(mode="absolute", value=prec_value, unit="ppm")
    frag_value = float(series.loc[pd.IndexSlice["msmsParamsArray", "msmsParams", "MatchTolerance", :]].squeeze())
    in_ppm = bool(series.loc[pd.IndexSlice["msmsParamsArray", "msmsParams", "MatchToleranceInPpm", :]].squeeze())
    fragment = MassTolerance(mode="absolute", value=frag_value, unit="ppm" if in_ppm else "Da")
    return precursor, fragment


def _min_peptide_length(series: pd.Series) -> int:
    """Read the minimum peptide length, tolerating the pre/post-rename key."""
    try:
        return int(series.loc["minPepLen"].squeeze())
    except KeyError:
        return int(series.loc["minPeptideLength"].squeeze())


def _mods_for_version(series: pd.Series, version: str) -> tuple[str, str]:
    """Homogenize fixed/variable modifications, handling the 1.6.0.0 path change."""
    if version > "1.6.0.0":
        fixed_path = pd.IndexSlice["parameterGroups", "parameterGroup", "fixedModifications", :]
    else:
        fixed_path = pd.IndexSlice["fixedModifications", :]
    fixed_mods = series.loc[fixed_path].squeeze()
    if not isinstance(fixed_mods, str):
        fixed_mods = ",".join(fixed_mods)

    variable_mods = series.loc[pd.IndexSlice["parameterGroups", "parameterGroup", "variableModifications", :]].squeeze()
    if not isinstance(variable_mods, str):
        variable_mods = ",".join(variable_mods)

    return _homogenize_mods(fixed_mods), _homogenize_mods(variable_mods)


def extract_params(
    source: Union[str, Path, IO[bytes], IO[str]],
    ms2frac: str = "FTMS",
) -> Parameters:
    """Parse a MaxQuant ``mqpar.xml`` into :class:`Parameters`.

    Mirrors ``proteobench.io.params.maxquant.extract_params``: MS2
    fragmentation method must be selected explicitly (``"FTMS"`` by
    default) because mqpar.xml carries one entry per fragmentation
    method.
    """
    record = _read_xml(source)
    record["msmsParamsArray"] = [
        d for d in record["msmsParamsArray"] if d["msmsParams"]["Name"] == ms2frac
    ]
    series = _build_series(record, 4).sort_index()

    version = str(series.loc["maxQuantVersion"].squeeze())
    precursor_tolerance, fragment_tolerance = _tolerance_pair(series)
    enzyme_mode = int(series.loc[("parameterGroups", "parameterGroup", "enzymeMode")].squeeze())
    fixed_mods, variable_mods = _mods_for_version(series, version)

    return Parameters(
        software_name="MaxQuant",
        software_version=version,
        search_engine="Andromeda",
        ident_fdr_psm=float(series.loc["peptideFdr"].squeeze()),
        ident_fdr_peptide=None,
        ident_fdr_protein=float(series.loc["proteinFdr"].squeeze()),
        enable_match_between_runs=series.loc["matchBetweenRuns"].squeeze().lower() == "true",
        precursor_mass_tolerance=precursor_tolerance,
        fragment_mass_tolerance=fragment_tolerance,
        enzyme=series.loc[("parameterGroups", "parameterGroup", "enzymes", "string")].squeeze(),
        semi_enzymatic=enzyme_mode != 0,
        allowed_miscleavages=int(series.loc[pd.IndexSlice["parameterGroups", "parameterGroup", "maxMissedCleavages", :]].squeeze()),
        min_peptide_length=_min_peptide_length(series),
        max_peptide_length=None,
        fixed_mods=fixed_mods,
        variable_mods=variable_mods,
        max_mods=int(series.loc[("parameterGroups", "parameterGroup", "maxNmods")].squeeze()),
        min_precursor_charge=None,
        max_precursor_charge=int(series.loc[pd.IndexSlice["parameterGroups", "parameterGroup", "maxCharge", :]].squeeze()),
    )
