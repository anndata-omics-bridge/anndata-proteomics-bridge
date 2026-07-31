"""Cached-corpus regressions for AlphaDIA ion coverage across three export shapes.

Pins the contract established in TODO_alphadia_coverage.md: all 11 cached AlphaDIA
submissions convert to the exact six-run observation axis their ProteoBench module
annotation declares, across three shape-distinct rule documents:

* ``v1_10`` — wide TSV, one bare run column per sample, features repeated verbatim
* ``v1_12`` — long TSV, ``run`` + ``intensity``
* ``v2``    — long parquet, dotted namespaces, ``raw.name`` + ``precursor.intensity``
"""

from __future__ import annotations

import tomllib
from dataclasses import dataclass
from pathlib import Path

import pytest
from conversion_support import convert_parameterized_level_to_anndata

from anndata_proteomics.converters.pipeline import PresentRuleVersion, resolve_parameters
from anndata_proteomics.readers.dispatch import read_table
from anndata_proteomics.rules.loader import resolve_rule_locator_for_version
from anndata_proteomics.rules.registry import RuleLocator

REPO_ROOT = Path(__file__).resolve().parents[1]
ROOT = REPO_ROOT / "test_data_download" / "json_dir"
ANNOTATIONS = REPO_ROOT / "test_data_download" / "annotations"


@dataclass(frozen=True)
class AlphaDiaCase:
    """One cached AlphaDIA submission and the contract it pins."""

    relative_input: str
    module: str
    expected_version: str
    expected_document: str

    @property
    def input_path(self) -> Path:
        return ROOT / self.relative_input

    @property
    def params_path(self) -> Path:
        return self.input_path.parent / "param_0..txt"


CASES = [
    # Shape A — wide TSV.
    AlphaDiaCase(
        "Results_quant_ion_DIA_Astral/11d1e3cad24d35b31110e82274bd0fef65435a34/input_file.tsv",
        "dia_astral",
        "1.10.3",
        "v1_10",
    ),
    AlphaDiaCase(
        "Results_quant_ion_DIA_diaPASEF/40e7248950fb9c487f3f7903b5a53068280acfee/input_file.tsv",
        "dia_diapasef",
        "1.10.3",
        "v1_10",
    ),
    AlphaDiaCase(
        "Results_quant_ion_DIA_Astral/059a69e4eccf183869f956c4c9494d31945f6849/input_file.tsv",
        "dia_astral",
        "1.10.4-dev0",
        "v1_10",
    ),
    AlphaDiaCase(
        "Results_quant_ion_DIA_diaPASEF/ff31788bc5a1c8a8da1907f916acc76c4bef4cef/input_file.tsv",
        "dia_diapasef",
        "1.10.4-dev0",
        "v1_10",
    ),
    # Shape B — long TSV.
    AlphaDiaCase(
        "Results_quant_ion_DIA_Astral/1bd69ca0cc5060987715afa2933b54f61fd1522e/input_file.tsv",
        "dia_astral",
        "1.12.1",
        "v1_12",
    ),
    AlphaDiaCase(
        "Results_quant_ion_DIA_diaPASEF/61e5aa4f0b50dc3a6f3fe8fafdd0aa131905a3a3/input_file.tsv",
        "dia_diapasef",
        "1.12.1",
        "v1_12",
    ),
    AlphaDiaCase(
        "Results_quant_ion_DIA_ZenoTOF/e465ee503ca0a029ce05047b514e20658238dc99/input_file.tsv",
        "dia_zenotof",
        "1.12.1",
        "v1_12",
    ),
    AlphaDiaCase(
        "Results_quant_ion_DIA_ZenoTOF/0d6d0437f5779b8be0790596b546af7ecc033b69/input_file.tsv",
        "dia_zenotof",
        "1.12.2",
        "v1_12",
    ),
    # Shape C — long parquet.
    AlphaDiaCase(
        "Results_quant_ion_DIA_Astral/a17457871aa99ec961cbeef60e81f2400f95fa4e/input_file.parquet",
        "dia_astral",
        "2.1.0",
        "v2",
    ),
    AlphaDiaCase(
        "Results_quant_ion_DIA_singlecell/1d6460eccc4bcd4a60b2c96f94ef35f5c9674599/input_file.parquet",
        "dia_singlecell",
        "2.1.0",
        "v2",
    ),
    AlphaDiaCase(
        "Results_quant_ion_DIA_Astral/39dad6944c982c86f0566488fc2f0737ba1ebb84/input_file.parquet",
        "dia_astral",
        "2.1.1",
        "v2",
    ),
]


def _expected_runs(module: str) -> set[str]:
    annotation = ANNOTATIONS / f"{module}.toml"
    with annotation.open("rb") as handle:
        return {sample["raw_file"] for sample in tomllib.load(handle)["samples"]}


@pytest.mark.parametrize(
    "case",
    CASES,
    ids=lambda case: f"{case.expected_version}-{case.module}",
)
def test_cached_alphadia_inputs_convert_with_the_annotated_run_axis(case: AlphaDiaCase) -> None:
    if not case.input_path.exists():
        pytest.skip(f"cached ProteoBench input is absent: {case.input_path}")

    resolution = resolve_parameters(case.params_path, "alphadia")
    assert resolution.version == PresentRuleVersion(case.expected_version)

    locator = resolve_rule_locator_for_version("alphadia", "ion", case.expected_version)
    assert isinstance(locator, RuleLocator)
    assert locator.path.parent.name == case.expected_document

    result = convert_parameterized_level_to_anndata(
        read_table(case.input_path),
        "alphadia",
        "ion",
        resolution,
    )

    assert set(result.obs_names) == _expected_runs(case.module)
    assert result.n_vars > 0
    assert "ProForma_peptidoform" in result.var
    assert "ProForma_peptide" in result.var


@pytest.mark.parametrize(
    "case",
    [c for c in CASES if c.expected_document == "v1_10"],
    ids=lambda case: f"{case.expected_version}-{case.module}",
)
def test_wide_shape_collapses_verbatim_repeated_features(case: AlphaDiaCase) -> None:
    """AlphaDIA 1.10 repeats a precursor's whole row, intensities included.

    ``keep_first`` is therefore correct and ``aggregate`` would multiply the
    intensity by the repeat count.
    """
    if not case.input_path.exists():
        pytest.skip(f"cached ProteoBench input is absent: {case.input_path}")

    frame = read_table(case.input_path)

    # The precursor identity is (sequence, mods, mod_sites, charge) — not the vendor's
    # mod_seq_charge_hash, which is a uint64 that pandas reads from TSV as float64 and
    # silently collapses (81 949 distinct values become 75 292). That is why the hash is
    # not carried as a var column. Measured before converting: conversion normalizes
    # modifications and renames selected columns on the frame in place.
    identity = ["sequence", "mods", "mod_sites", "charge"]
    keys = frame[identity].fillna("")
    expected_features = len(keys.drop_duplicates())
    repeated_rows_are_exact_copies = frame[keys.duplicated(keep=False)].duplicated(keep=False).all()

    resolution = resolve_parameters(case.params_path, "alphadia")
    result = convert_parameterized_level_to_anndata(
        frame,
        "alphadia",
        "ion",
        resolution,
    )

    assert result.n_vars == expected_features
    assert result.var_names.is_unique
    # Every repeat is an exact copy of the whole row, intensities included, so
    # keep_first loses nothing and aggregate would multiply the intensity.
    assert repeated_rows_are_exact_copies


def test_peptidoforms_of_one_sequence_stay_distinct_features() -> None:
    """A sequence carrying different modifications must not collapse into one feature."""
    case = next(c for c in CASES if c.expected_document == "v1_12")
    if not case.input_path.exists():
        pytest.skip(f"cached ProteoBench input is absent: {case.input_path}")

    resolution = resolve_parameters(case.params_path, "alphadia")
    result = convert_parameterized_level_to_anndata(
        read_table(case.input_path),
        "alphadia",
        "ion",
        resolution,
    )

    var = result.var
    assert var["ProForma_peptidoform"].nunique() > var["ProForma_peptide"].nunique()
    modified = var["ProForma_peptidoform"].str.contains("[UNIMOD:", regex=False)
    assert modified.any()
    assert result.var_names.is_unique
