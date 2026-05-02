"""End-to-end CLI tests via subprocess.

These run the installed `anndata-proteomics` binary, exercising the cyclopts
dispatch layer. Complements tests/test_cli.py which calls subcommand functions
directly (unit-level).
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from anndata_proteomics.rules.registry import find_rule

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
    assert "PASS" in r.stdout
    assert "0 failed" in r.stdout


def test_cli_validate_path_happy() -> None:
    p = find_rule("diann", "ion")
    r = _run("validate", str(p))
    assert r.returncode == 0, r.stderr
    assert "PASS" in r.stdout


def test_cli_validate_path_bad(tmp_path: Path) -> None:
    bad = tmp_path / "bad.toml"
    bad.write_text("not = valid [[")
    r = _run("validate", str(bad))
    assert r.returncode == 1
    assert "FAIL" in r.stdout
    assert "1 failed" in r.stdout


def test_cli_list_outputs_six_rules() -> None:
    r = _run("list")
    assert r.returncode == 0, r.stderr
    lines = [line for line in r.stdout.splitlines() if line.strip()]
    assert len(lines) == 6
    assert "diann" in r.stdout
    assert "wombat" in r.stdout


def test_cli_convert_stub_returns_two() -> None:
    r = _run("convert", "dummy.tsv", "rule.toml")
    assert r.returncode == 2
    assert "not yet implemented" in (r.stdout + r.stderr).lower()
