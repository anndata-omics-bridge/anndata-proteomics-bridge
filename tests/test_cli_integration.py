"""End-to-end CLI tests via subprocess.

These run the installed `anndata-proteomics` binary, exercising the cyclopts
dispatch layer. Complements tests/test_cli.py which calls subcommand functions
directly (unit-level).

CLI output goes through loguru → stderr; tests assert against `r.stderr`.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from anndata_proteomics.rules.registry import find_rule
from anndata_proteomics.test_data import find_test_data

# The console script lives next to the python that's running pytest.
# Resolving via sys.executable means tests work regardless of whether the venv
# is activated in the parent shell.
_CLI = str(Path(sys.executable).parent / "anndata-proteomics")


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


def test_cli_convert_writes_h5ad(tmp_path: Path) -> None:
    """Pick a small vendor file; convert via the binary; assert .h5ad exists."""
    import pytest

    data_file = find_test_data("WOMBAT")
    if data_file is None or not data_file.exists():
        pytest.skip("no WOMBAT test data available")

    out = tmp_path / "wombat.h5ad"
    r = _run("convert", str(data_file), "--output", str(out))
    assert r.returncode == 0, r.stderr
    assert out.exists()
    assert "wrote" in r.stderr
