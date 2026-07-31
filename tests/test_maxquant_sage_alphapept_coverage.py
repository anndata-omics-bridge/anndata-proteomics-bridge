"""Cached-corpus regressions for MaxQuant, Sage, and AlphaPept ion coverage.

Pins the contract established in TODO_maxquant_sage_alphapept_dda_coverage.md: every cached
MaxQuant submission converts across four modules despite configuration-dependent columns,
and the two new vendor documents convert with the exact six-run observation axis the
ProteoBench module annotations expect.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pytest
from conversion_support import convert_parameterized_level_to_anndata

from anndata_proteomics.converters.pipeline import PresentRuleVersion, resolve_parameters
from anndata_proteomics.params.model import Parameters
from anndata_proteomics.readers.dispatch import read_table
from anndata_proteomics.rules.loader import (
    RuleLocatorUnavailable,
    load_packaged_rule_for_version,
    resolve_parameterized_rule_locator_for_version,
    resolve_parameterized_rule_locator_without_version,
    resolve_rule_locator_for_version,
)
from anndata_proteomics.rules.registry import RuleLocator

ROOT = Path(__file__).resolve().parents[1] / "test_data_download" / "json_dir"


@dataclass(frozen=True)
class CorpusCase:
    """One cached vendor input whose conversion contract is pinned here."""

    slug: str
    relative_input: str
    expected_version: str

    @property
    def input_path(self) -> Path:
        return ROOT / self.relative_input


MAXQUANT_CASES = [
    CorpusCase(
        "maxquant",
        "Results_quant_ion_DDA/8f4fa9a7dd1f44ac4ae7a7e7fb9b9606660f4578/input_file.txt",
        "1.5.2.8",
    ),
    CorpusCase(
        "maxquant",
        "Results_quant_ion_DDA/f4f23a743baef55ba419a1cf0e8dd67a4cb5b7ac/input_file.txt",
        "1.5.3.30",
    ),
    CorpusCase(
        "maxquant",
        "Results_quant_ion_DDA/254d6c77ce656888918e738772ad5f5f6f1543e4/input_file.txt",
        "1.5.8.2",
    ),
    CorpusCase(
        "maxquant",
        "Results_quant_ion_DDA/7912158e0522a315917e20fe434966ca9a17192e/input_file.txt",
        "1.6.3.3",
    ),
    CorpusCase(
        "maxquant",
        "Results_quant_ion_DDA/a1140a31b414d7b3110ee9b9c0456cc4f1709782/input_file.txt",
        "2.1.3.0",
    ),
    CorpusCase(
        "maxquant",
        "Results_quant_ion_DDA/36b7b01b380f641722b3b34633bb53d72348eb80/input_file.txt",
        "2.1.4.0",
    ),
    CorpusCase(
        "maxquant",
        "Results_quant_ion_DDA/0280a06fabdbe84746419d0810deae56e7ab2406/input_file.txt",
        "2.3.1.0",
    ),
    CorpusCase(
        "maxquant",
        "Results_quant_ion_DDA/a3d801fcb75c46b2e76fa7078ae0a004360ebe44/input_file.txt",
        "2.4.13.0",
    ),
    CorpusCase(
        "maxquant",
        "Results_quant_ion_DDA/00e2f863939301a2a71178652972dad895b27520/input_file.txt",
        "2.5.1.0",
    ),
]

NEW_VENDOR_CASES = [
    CorpusCase(
        "sage",
        "Results_quant_ion_DDA/05906d50cf20634dc684bda319e3954b8cc2bf5b/input_file.tsv",
        "0.15.0-beta.2",
    ),
    CorpusCase(
        "alphapept",
        "Results_quant_ion_DDA/e8e80290fb48ff02de5ee54eb6b0114ff661bace/input_file.txt",
        "0.5.0",
    ),
    CorpusCase(
        "alphapept",
        "Results_quant_ion_DDA/94c4c0b7d00761d24fdde05276053b087cb99ea1/input_file.txt",
        "0.5.3",
    ),
]


def _convert(case: CorpusCase):
    if not case.input_path.exists():
        pytest.skip(f"cached ProteoBench input is absent: {case.input_path}")
    param_paths = sorted(case.input_path.parent.glob("param_0.*"))
    if not param_paths:
        pytest.skip(f"cached parameter file is absent beside {case.input_path}")
    param_path = param_paths[0]
    resolution = resolve_parameters(param_path, case.slug)
    assert resolution.version == PresentRuleVersion(case.expected_version)
    return convert_parameterized_level_to_anndata(
        read_table(case.input_path),
        case.slug,
        "ion",
        resolution,
    )


@pytest.mark.parametrize(
    "case",
    MAXQUANT_CASES + NEW_VENDOR_CASES,
    ids=lambda case: f"{case.slug}-{case.expected_version}",
)
def test_cached_dda_inputs_convert_with_exact_run_axis(case: CorpusCase) -> None:
    result = _convert(case)

    assert result.n_obs == 6
    assert result.n_vars > 0
    assert all(name.startswith("LFQ_Orbitrap_DDA_Condition_") for name in result.obs_names)
    # The module annotation joins on bare run names; a surviving extension breaks it.
    assert not any(name.endswith((".raw", ".mzML", ".hdf")) for name in result.obs_names)
    metadata = result.uns["anndata_proteomics"]
    assert metadata["rule_selection_method"] == "software_version"


@pytest.mark.parametrize(
    "case",
    MAXQUANT_CASES,
    ids=lambda case: case.expected_version,
)
def test_maxquant_converts_without_fraction_column(case: CorpusCase) -> None:
    """No cached DDA `evidence.txt` carries `Fraction`; it must not gate conversion."""
    frame = read_table(case.input_path) if case.input_path.exists() else None
    if frame is None:
        pytest.skip(f"cached ProteoBench input is absent: {case.input_path}")
    assert "Fraction" not in frame.columns

    result = _convert(case)

    assert "Fraction" not in result.obs.columns
    assert result.var["Proteins"].notna().any()


def test_maxquant_1_5_2_8_uses_title_case_leading_protein_columns() -> None:
    """1.5.2.8 spells the leading-protein columns in title case and has no `Experiment`."""
    case = MAXQUANT_CASES[0]
    if not case.input_path.exists():
        pytest.skip(f"cached ProteoBench input is absent: {case.input_path}")
    columns = set(read_table(case.input_path).columns)
    assert "Leading Proteins" in columns
    assert "Leading proteins" not in columns
    assert "Experiment" not in columns

    result = _convert(case)

    assert "Experiment" not in result.obs.columns
    assert "Leading_Proteins_Legacy" in result.var.columns
    assert "Leading_Proteins" not in result.var.columns


def test_maxquant_optional_columns_stay_declared() -> None:
    """The volatile columns move to `optional_select`; only `Raw file` stays required."""
    rule = load_packaged_rule_for_version("maxquant", "ion", "2.6.7.0")

    assert rule.columns.obs.select == {"Raw_File": "Raw file"}
    assert rule.columns.obs.optional_select == {
        "Experiment": "Experiment",
        "Fraction": "Fraction",
    }


@pytest.mark.parametrize(
    ("relative_input", "version"),
    [
        (
            "Results_quant_ion_DDA_Astral/e5a709af0d04b71d0572029df66753a06195a0f6/input_file.txt",
            "2.6.7.0",
        ),
        (
            "Results_quant_ion_DIA_Astral/be39a8defd1e85b69bcad3896a82e721e1bdfb0c/input_file.txt",
            "2.6.3.0",
        ),
        (
            "Results_quant_ion_DIA_singlecell/601cf4e2ba80aa84103611eec97353f0cabc7604/input_file.txt",
            "2.7.5.0",
        ),
    ],
    ids=["dda_astral-2.6.7.0", "dia_astral-2.6.3.0", "dia_singlecell-2.7.5.0"],
)
def test_maxquant_other_modules_do_not_regress(relative_input: str, version: str) -> None:
    """2.6.7.0 was the only supported version before; the broadened regex must keep all three."""
    case = CorpusCase("maxquant", relative_input, version)

    result = _convert(case)

    assert result.n_obs > 0
    assert result.n_vars > 0
    assert "Experiment" in result.obs.columns or version == "2.6.3.0"


SAGE_COLLAPSED = ROOT / "Results_quant_ion_DDA_Astral/fc1bf3c26e8323b1e724c1b6756b3f1a80941ab9"


def test_sage_level_follows_combine_charge_states_not_the_version() -> None:
    """Sage's level is decided by `lfq_settings.combine_charge_states`, not by version.

    Sage's DOCS.md defaults it to `true` ("Combine all charge states for quantification"),
    and a combined row is written with `charge = -1`
    (`crates/sage-cli/src/runner.rs`: `charge.unwrap_or(-1)`), so the same `lfq.tsv` schema
    is ion-level or peptidoform-level depending only on that setting. Neither the version
    regex nor the headers can separate the two.
    """
    collapsed = Parameters(software_name="Sage", combine_charge_states=True)
    resolved = Parameters(software_name="Sage", combine_charge_states=False)

    for version in ("0.14.6", "0.15.0", "0.15.0-beta.2", "0.16.0"):
        assert isinstance(
            resolve_parameterized_rule_locator_for_version("sage", "ion", version, collapsed),
            RuleLocatorUnavailable,
        )
        assert isinstance(
            resolve_parameterized_rule_locator_for_version(
                "sage", "peptidoform", version, collapsed
            ),
            RuleLocator,
        )
        assert isinstance(
            resolve_parameterized_rule_locator_for_version("sage", "ion", version, resolved),
            RuleLocator,
        )
        assert isinstance(
            resolve_parameterized_rule_locator_for_version(
                "sage", "peptidoform", version, resolved
            ),
            RuleLocatorUnavailable,
        )

    assert isinstance(
        resolve_parameterized_rule_locator_without_version("sage", "ion", collapsed),
        RuleLocatorUnavailable,
    )
    assert isinstance(
        resolve_parameterized_rule_locator_without_version("sage", "peptidoform", collapsed),
        RuleLocator,
    )
    assert isinstance(
        resolve_parameterized_rule_locator_without_version("sage", "ion", resolved),
        RuleLocator,
    )
    assert isinstance(
        resolve_parameterized_rule_locator_without_version("sage", "peptidoform", resolved),
        RuleLocatorUnavailable,
    )


def test_sage_gated_levels_are_unavailable_without_parameters() -> None:
    """No parameters means no charge-state decision: offer nothing rather than guess wrong."""
    assert isinstance(
        resolve_rule_locator_for_version("sage", "ion", "0.15.0"),
        RuleLocatorUnavailable,
    )
    assert isinstance(
        resolve_rule_locator_for_version("sage", "peptidoform", "0.15.0"),
        RuleLocatorUnavailable,
    )


def test_sage_charge_collapsed_export_converts_as_peptidoform() -> None:
    """The 0.14.6 Astral submission is peptidoform-level and now converts as one."""
    case = CorpusCase("sage", str(SAGE_COLLAPSED.relative_to(ROOT) / "input_file.tsv"), "0.14.6")
    if not case.input_path.exists():
        pytest.skip(f"cached ProteoBench input is absent: {case.input_path}")
    param_path = sorted(case.input_path.parent.glob("param_0.*"))[0]
    resolution = resolve_parameters(param_path, "sage")
    assert resolution.parameters is not None
    assert resolution.parameters.combine_charge_states is True

    result = convert_parameterized_level_to_anndata(
        read_table(case.input_path),
        "sage",
        "peptidoform",
        resolution,
    )

    assert result.n_obs == 6
    assert result.n_vars == 31200
    # The -1 sentinel must not reach the output: no charge column, no charge in the feature key.
    assert "Charge" not in result.var.columns
    assert not any("/" in name for name in result.var_names)


def test_sage_charge_resolved_export_is_not_reachable_as_peptidoform() -> None:
    """The charge-resolved 0.15.0-beta.2 submission stays ion-level only."""
    fixture = ROOT / "Results_quant_ion_DDA/05906d50cf20634dc684bda319e3954b8cc2bf5b"
    param_paths = sorted(fixture.glob("param_0.*"))
    if not param_paths:
        pytest.skip(f"cached parameter file is absent under {fixture}")
    resolution = resolve_parameters(param_paths[0], "sage")
    assert resolution.parameters is not None

    assert resolution.parameters.combine_charge_states is False
    assert isinstance(resolution.version, PresentRuleVersion)
    assert isinstance(
        resolve_parameterized_rule_locator_for_version(
            "sage",
            "peptidoform",
            resolution.version.value,
            resolution.parameters,
        ),
        RuleLocatorUnavailable,
    )


def test_alphapept_intensity_layer_is_the_proteobench_column() -> None:
    """`ms1_int_sum_apex_dn` reproduces ProteoBench's per-run intensities; the others do not."""
    rule = load_packaged_rule_for_version("alphapept", "ion", "0.5.3")
    x_layer = next(layer for layer in rule.layers if layer.name == rule.axis.x_layer)

    assert x_layer.source == "ms1_int_sum_apex_dn"
    assert rule.axis.duplicates.mode == "keep_first"


def test_alphapept_txt_is_detected_as_comma_delimited() -> None:
    case = NEW_VENDOR_CASES[1]
    if not case.input_path.exists():
        pytest.skip(f"cached ProteoBench input is absent: {case.input_path}")

    assert len(read_table(case.input_path).columns) == 85
