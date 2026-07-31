"""Focused edge-path coverage for parameter models and shared parsers."""

from __future__ import annotations

from io import BytesIO
from pathlib import Path
from typing import Any, cast

import anndata as ad
import numpy as np
import pytest
from pydantic import ValidationError

from anndata_proteomics.adapters.anndata.params import (
    has_search_parameters,
    write_search_parameters,
)
from anndata_proteomics.modifications.model import SearchedModification
from anndata_proteomics.params import model
from anndata_proteomics.params.model import MassTolerance, Parameters, Probability
from anndata_proteomics.params.parsers import _common, alphapept, msaid, wombat


class _Unseekable:
    def seekable(self) -> bool:
        return False

    def read(self) -> str:
        return "content"


def test_probability_and_tolerance_edge_inputs() -> None:
    assert Probability.model_validate({"value": "50%"}).value == 0.5
    assert Probability.model_validate({"value": "0.25"}).value == 0.25
    assert Probability.parse(object()) is None

    with pytest.raises(ValidationError, match="requires value"):
        MassTolerance(mode="absolute", unit="ppm")
    with pytest.raises(ValidationError, match="cannot define unit"):
        MassTolerance(mode="automatic", unit="ppm")
    with pytest.raises(ValidationError, match="cannot define numeric"):
        MassTolerance(mode="automatic", value=1)
    assert MassTolerance.parse(0) == MassTolerance(
        mode="automatic",
        label="Automatic calibration",
    )
    with pytest.raises(ValueError, match="explicit unit"):
        MassTolerance.parse(1)
    with pytest.raises(TypeError, match="unsupported"):
        MassTolerance.parse(object())
    with pytest.raises(ValueError, match="could not parse"):
        MassTolerance.parse("not a tolerance")


def test_parameter_coercion_edge_inputs() -> None:
    with pytest.raises(ValidationError):
        Parameters.model_validate({"enzyme": object()})
    assert Parameters(enable_match_between_runs=None).enable_match_between_runs is None
    assert (
        Parameters.model_validate({"enable_match_between_runs": 1}).enable_match_between_runs
        is True
    )
    with pytest.raises(ValidationError, match="cannot coerce boolean"):
        Parameters.model_validate({"enable_match_between_runs": "maybe"})
    with pytest.raises(ValidationError, match="cannot coerce boolean"):
        Parameters.model_validate({"enable_match_between_runs": object()})

    assert Parameters.model_validate({"scan_window": 3.0}).scan_window == 3
    assert Parameters.model_validate({"scan_window": 3.5}).scan_window == "3.5"
    with pytest.raises(ValidationError):
        Parameters.model_validate({"scan_window": object()})

    modification = SearchedModification(name="Oxidation")
    assert Parameters.model_validate({"fixed_mods": modification}).fixed_mods == [modification]
    dict_mods = Parameters.model_validate({"fixed_mods": {"C": "57.0"}}).fixed_mods
    assert dict_mods[0].target == "C"
    with pytest.raises(TypeError, match="unsupported modification"):
        Parameters.model_validate({"fixed_mods": 1})
    assert Parameters.model_validate({"unparsed_parameters": None}).unparsed_parameters == []


def test_parameter_low_level_serialization_helpers() -> None:
    with pytest.raises(ValueError, match="requires unit"):
        model._normalize_unit("")
    assert model._coerce_float(None) is None
    assert model._coerce_float(object()) is None
    assert model._split_mod_string(" ") == []

    modification = SearchedModification(name="Oxidation")
    assert model._modification_from_item(modification) is modification
    marker = object()
    assert model._to_scalar(marker) == str(marker)


def test_shared_parser_io_and_tolerance_errors(tmp_path: Path) -> None:
    path = tmp_path / "text.txt"
    path.write_text("path", encoding="utf-8")
    assert _common.read_text(path) == "path"
    assert _common.read_text(BytesIO(b"bytes")) == "bytes"
    assert _common.read_text(cast(Any, _Unseekable())) == "content"

    with pytest.raises(ValueError, match="must be dict"):
        _common.format_tolerance_range(cast(Any, []))
    with pytest.raises(KeyError, match="unsupported"):
        _common.format_tolerance_range({"invalid": [1, 2]})
    assert _common.format_tolerance_range({"invalid": [0], "ppm": [-1, 1]}) == "[-1 ppm, 1 ppm]"
    assert _common.homogenize_paren_mods("known", {"known": "mapped"}) == "mapped"
    assert _common.lookup_mass_mod(1.0, {1.0: "known"}) == (_common.MassModificationMatch("known"))
    assert _common.lookup_mass_mod(2.0, {1.0: "known"}) == (
        _common.UnrecognizedModificationMass(2.0)
    )


def test_vendor_specific_small_helpers() -> None:
    assert alphapept._map_modifications(["cC", "oxM"]) == ("C[Carbamidomethyl], M[Oxidation]")
    assert msaid._homogenize_mods(" ") == ""
    assert wombat._homogenize_mod_xtandem("plain") == "plain"
    assert wombat._homogenize_mod_xtandem("Acetyl of N-term") == ("N-term[Acetyl]")
    assert wombat._homogenize_mod_xtandem("Acetyl of Protein C-term") == ("Protein C-term[Acetyl]")
    assert wombat._homogenize_mod_xtandem("Acetyl of C-term") == ("C-term[Acetyl]")


def test_parameter_storage_without_payload_or_source_path() -> None:
    target = ad.AnnData(np.ones((1, 1)))
    target.uns["anndata_proteomics"] = {"other": "metadata"}
    assert not has_search_parameters(target)
    write_search_parameters(target, Parameters())
    assert "search_parameters_path" not in target.uns["anndata_proteomics"]
