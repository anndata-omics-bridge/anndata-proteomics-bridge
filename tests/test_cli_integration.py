"""End-to-end CLI tests via subprocess.

These run the installed `apb` binary, exercising the cyclopts
dispatch layer. Complements tests/test_cli.py which calls subcommand functions
directly (unit-level).

CLI output goes through loguru → stderr; tests assert against `r.stderr`.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from anndata_proteomics.rules.registry import find_rule
from anndata_proteomics.test_data import find_param_file, find_test_data

# The console script lives next to the python that's running pytest.
# Resolving via sys.executable means tests work regardless of whether the venv
# is activated in the parent shell.
_CLI = str(Path(sys.executable).parent / "apb")


def _run(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [_CLI, *args],
        capture_output=True,
        text=True,
        check=False,
    )


def test_cli_validate_no_args_returns_zero() -> None:
    r = _run("validate")
    assert r.returncode == 0, r.stderr
    assert "PASS" in r.stderr
    assert "0 failed" in r.stderr


def test_cli_validate_path_happy() -> None:
    p = find_rule("diann", "ion")
    r = _run("validate", str(p))
    assert r.returncode == 0, r.stderr
    assert "PASS" in r.stderr


def test_cli_validate_path_bad(tmp_path: Path) -> None:
    bad = tmp_path / "bad.toml"
    bad.write_text("not = valid [[")
    r = _run("validate", str(bad))
    assert r.returncode == 1
    assert "FAIL" in r.stderr
    assert "1 failed" in r.stderr


def test_cli_list_outputs_eleven_rules() -> None:
    r = _run("list")
    assert r.returncode == 0, r.stderr
    lines = [line for line in r.stderr.splitlines() if line.strip()]
    assert len(lines) == 11
    assert "diann" in r.stderr
    assert "wombat" in r.stderr


def _require(software: str):
    """Return (data_file, param_file) for a tool, or skip when either is unavailable."""
    import pytest

    data_file = find_test_data(software)
    param_file = find_param_file(software)
    if data_file is None or not data_file.exists():
        pytest.skip(f"no {software} test data available")
    if param_file is None or not param_file.exists():
        pytest.skip(f"no {software} param fixture available")
    return data_file, param_file


def test_cli_convert_with_rule_toml_writes_h5ad(tmp_path: Path) -> None:
    """The --rule-toml override (single level, version-agnostic) writes a .h5ad."""
    import pytest

    data_file = find_test_data("WOMBAT")
    if data_file is None or not data_file.exists():
        pytest.skip("no WOMBAT test data available")
    rule = find_rule("wombat", "peptidoform")
    out = tmp_path / "wombat.h5ad"
    r = _run("convert", str(data_file), "--rule-toml", str(rule), "--output", str(out))
    assert r.returncode == 0, r.stderr
    assert out.exists()
    assert "wrote" in r.stderr


def test_cli_convert_default_writes_h5mu(tmp_path: Path) -> None:
    """A multi-level vendor (DIA-NN) with no level argument writes a multi-modality .h5mu."""
    import mudata

    data_file, param_file = _require("DIA-NN")
    out = tmp_path / "diann.h5mu"
    r = _run("convert", str(data_file), "--params", str(param_file), "--output", str(out))
    assert r.returncode == 0, r.stderr
    assert out.exists()
    md = mudata.read_h5mu(out)
    assert len(md.mod) >= 2


def test_cli_convert_explicit_level_writes_h5ad(tmp_path: Path) -> None:
    """A single level argument (DIA-NN protein) writes a .h5ad."""
    data_file, param_file = _require("DIA-NN")
    out = tmp_path / "diann_protein.h5ad"
    r = _run(
        "convert", str(data_file), "protein", "--params", str(param_file), "--output", str(out)
    )
    assert r.returncode == 0, r.stderr
    assert out.exists()
    assert "wrote" in r.stderr
