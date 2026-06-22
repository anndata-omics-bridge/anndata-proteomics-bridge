"""Verify every packaged parse_*.toml validates and matches its filename."""

from __future__ import annotations

import tomllib
from pathlib import Path

import pytest

from anndata_proteomics.rules.registry import iter_packaged_rules
from anndata_proteomics.rules.schema import ParseRule
from anndata_proteomics.rules.validate import validate_all_packaged


def test_all_packaged_rules_validate() -> None:
    results = validate_all_packaged()
    failed = [r for r in results if not r.ok]
    assert not failed, "\n".join(f"{r.path}: {r.error}" for r in failed)


def test_at_least_one_long_and_one_wide_rule() -> None:
    rules = [
        ParseRule.model_validate(tomllib.loads(p.read_text()))
        for p in iter_packaged_rules()
    ]
    shapes = {r.input_shape for r in rules}
    assert "long" in shapes, f"no long rule found: {shapes}"
    assert "wide" in shapes, f"no wide rule found: {shapes}"


@pytest.mark.parametrize(
    "toml_path",
    list(iter_packaged_rules()),
    ids=lambda p: f"{p.parent.name}/{p.name}",
)
def test_filename_quant_level_matches_toml(toml_path: Path) -> None:
    # parse_<software>_<level>.toml (version-foldered) or parse_<software>_<level>_<n>.toml (flat).
    parts = toml_path.stem.split("_")
    if parts[-1].isdigit():  # drop a trailing flat-file version token
        parts = parts[:-1]
    filename_level = parts[-1]
    rule = ParseRule.model_validate(tomllib.loads(toml_path.read_text()))
    assert rule.quantification_level == filename_level, (
        f"{toml_path.name}: filename says level={filename_level!r} "
        f"but TOML has quantification_level={rule.quantification_level!r}"
    )
