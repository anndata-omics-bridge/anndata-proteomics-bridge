"""Validate every packaged parse_*.toml against the ParseRule schema."""

from __future__ import annotations

import tomllib
from pathlib import Path

import pytest

from anndata_proteomics.rules.schema import ParseRule


PARSING_RULES_ROOT = (
    Path(__file__).resolve().parent.parent
    / "src"
    / "anndata_proteomics"
    / "parsing_rules"
)


def _packaged_rule_files() -> list[Path]:
    return sorted(PARSING_RULES_ROOT.glob("*/parse_*.toml"))


@pytest.mark.parametrize(
    "toml_path",
    _packaged_rule_files(),
    ids=lambda p: f"{p.parent.name}/{p.name}",
)
def test_packaged_rule_validates(toml_path: Path) -> None:
    data = tomllib.loads(toml_path.read_text())
    ParseRule.model_validate(data)


def test_at_least_one_long_and_one_wide_rule() -> None:
    files = _packaged_rule_files()
    shapes = {
        f.parent.name: ParseRule.model_validate(tomllib.loads(f.read_text())).input_shape
        for f in files
    }
    assert "long" in shapes.values(), f"no long rule found: {shapes}"
    assert "wide" in shapes.values(), f"no wide rule found: {shapes}"


@pytest.mark.parametrize(
    "toml_path",
    _packaged_rule_files(),
    ids=lambda p: f"{p.parent.name}/{p.name}",
)
def test_filename_quant_level_matches_toml(toml_path: Path) -> None:
    # Filename convention: parse_<software_tokens...>_<level>_<file_version>.toml
    # The level is always the second-to-last underscore-separated token of the stem.
    parts = toml_path.stem.split("_")
    filename_level = parts[-2]
    rule = ParseRule.model_validate(tomllib.loads(toml_path.read_text()))
    assert rule.quantification_level == filename_level, (
        f"{toml_path.name}: filename says level={filename_level!r} "
        f"but TOML has quantification_level={rule.quantification_level!r}"
    )
