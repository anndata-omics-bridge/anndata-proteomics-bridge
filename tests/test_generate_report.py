"""Smoke tests for tools/generate_report.py.

Skips gracefully when Rscript / annProtSum / the test_data_download cache are
absent (the orchestrator pulls inputs from the cache and shells out to
render_report.R).
"""

from __future__ import annotations

import importlib.util
import json
import shutil
import sys
from pathlib import Path

import pytest

from anndata_proteomics.test_data import find_test_data

PROJECT_ROOT = Path(__file__).resolve().parent.parent
TOOL_PATH = PROJECT_ROOT / "tools" / "generate_report.py"


def _load_tool():
    """Import tools/generate_report.py as a module."""
    spec = importlib.util.spec_from_file_location("generate_report", TOOL_PATH)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules["generate_report"] = mod
    spec.loader.exec_module(mod)
    return mod


def _have_render_script() -> bool:
    """Quick check that Rscript is installed and the render script is reachable."""
    if shutil.which("Rscript") is None:
        return False
    mod = _load_tool()
    try:
        mod._resolve_render_script()
        return True
    except RuntimeError:
        return False


def test_rebuild_index_with_no_meta_files_writes_empty_table(tmp_path: Path) -> None:
    mod = _load_tool()
    index = mod.rebuild_index(tmp_path)
    assert index.exists()
    body = index.read_text()
    assert "input size" in body
    assert ".h5ad size" in body
    assert "no conversions yet" in body


def test_rebuild_index_renders_file_sizes(tmp_path: Path) -> None:
    mod = _load_tool()
    meta = {
        "status": "ok",
        "software": "ExampleSoft",
        "stem": "examplesoft_12345678",
        "input_path": "/tmp/input.tsv",
        "h5ad_path": "examplesoft_12345678.h5ad",
        "html_path": "examplesoft_12345678.html",
        "log_path": "examplesoft_12345678.log",
        "input_size_bytes": 1536,
        "h5ad_size_bytes": 2 * 1024 * 1024,
        "n_obs": 2,
        "n_var": 3,
        "layers": [{"name": "abundance", "n_obs": 2, "n_var": 3}],
        "error": None,
    }
    (tmp_path / "examplesoft_12345678.meta.json").write_text(json.dumps(meta) + "\n")

    index = mod.rebuild_index(tmp_path)
    body = index.read_text()

    assert "1.5 KiB" in body
    assert "2.0 MiB" in body


def test_main_runs_one_converter_end_to_end(tmp_path: Path) -> None:
    """Run main() against the WOMBAT rule; check artifacts + meta + index."""
    if find_test_data("WOMBAT") is None:
        pytest.skip("test_data_download cache not present")
    if not _have_render_script():
        pytest.skip("annProtSum render_report.R not reachable")

    mod = _load_tool()
    rc = mod.main(["--rule", "WOMBAT", "--output-dir", str(tmp_path)])
    assert rc == 0

    metas = list(tmp_path.glob("*.meta.json"))
    assert len(metas) == 1
    meta = json.loads(metas[0].read_text())
    assert meta["status"] == "ok"
    assert meta["software"] == "WOMBAT"
    assert meta["input_size_bytes"] is not None
    assert meta["input_size_bytes"] > 0
    assert meta["h5ad_size_bytes"] is not None
    assert meta["h5ad_size_bytes"] > 0
    assert (tmp_path / meta["h5ad_path"]).exists()
    assert (tmp_path / meta["html_path"]).exists()
    assert (tmp_path / meta["log_path"]).exists()

    index = tmp_path / "index.html"
    assert index.exists()
    body = index.read_text()
    assert meta["h5ad_path"] in body
    assert meta["html_path"] in body
    assert meta["log_path"] in body
    assert "WOMBAT" in body


def test_main_emits_skipped_row_when_cache_missing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """If find_test_data returns None for every software, every row is `skipped`."""
    mod = _load_tool()
    monkeypatch.setattr(mod, "find_test_data", lambda _name: None)
    rc = mod.main(["--rule", "WOMBAT", "--output-dir", str(tmp_path)])
    # rc == 0 when there are no failures; skipped is not a failure.
    assert rc == 0

    metas = list(tmp_path.glob("*.meta.json"))
    assert len(metas) == 1
    meta = json.loads(metas[0].read_text())
    assert meta["status"] == "skipped"
    assert meta["h5ad_path"] is None
    assert meta["html_path"] is None
    assert meta["input_size_bytes"] is None
    assert meta["h5ad_size_bytes"] is None
    assert (tmp_path / meta["log_path"]).exists()
