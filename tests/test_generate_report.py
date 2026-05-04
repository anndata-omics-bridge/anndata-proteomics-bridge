"""Smoke tests for tools/generate_report.py.

Skips gracefully when Rscript or annProtSum isn't installed (the orchestrator
shells out to render_report.R).
"""

from __future__ import annotations

import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

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
    assert "no conversions yet" in index.read_text()


def test_convert_one_writes_artifacts(tmp_path: Path) -> None:
    if not _have_render_script():
        pytest.skip("annProtSum render_report.R not reachable (Rscript or package missing)")

    mod = _load_tool()

    # Synthesise a tiny long DIA-NN-shaped DataFrame so the packaged DIA-NN
    # rule can be auto-recognized.
    data_path = tmp_path / "tiny.tsv"
    pd.DataFrame(
        {
            "Run": ["S1", "S1", "S2", "S2"],
            "Modified.Sequence": ["P1", "P2", "P1", "P2"],
            "Stripped.Sequence": ["P1", "P2", "P1", "P2"],
            "Precursor.Charge": [2, 2, 2, 2],
            "Precursor.Id": ["P1_2", "P2_2", "P1_2", "P2_2"],
            "Protein.Group": ["A", "B", "A", "B"],
            "Protein.Ids": ["A", "B", "A", "B"],
            "Protein.Names": ["A", "B", "A", "B"],
            "Genes": ["g1", "g2", "g1", "g2"],
            "Precursor.Normalised": [10.0, 20.0, 11.0, 21.0],
            "Precursor.Quantity": [100.0, 200.0, 110.0, 210.0],
            "Ms1.Area": [1000.0, 2000.0, 1100.0, 2100.0],
            "Q.Value": [0.01, 0.02, 0.01, 0.02],
            "RT": [10.0, 20.0, 10.5, 20.5],
        }
    ).to_csv(data_path, sep="\t", index=False)

    out_dir = tmp_path / "out"
    conv = mod.convert_one(data_path, rule_toml=None, output_dir=out_dir)

    assert conv.h5ad_path.exists()
    assert conv.html_path.exists()
    assert conv.meta_path.exists()
    assert conv.html_path.stat().st_size > 1000  # non-trivial HTML

    meta = json.loads(conv.meta_path.read_text())
    assert meta["software"] == "DIA-NN"
    assert meta["n_obs"] == 2
    assert meta["n_var"] == 2
    assert {l["name"] for l in meta["layers"]} == {
        "Precursor_Normalised", "Precursor_Quantity", "Ms1_Area", "Q_Value", "RT"
    }

    # Rebuild index — single row referencing this conversion.
    index = mod.rebuild_index(out_dir)
    body = index.read_text()
    assert conv.h5ad_path.name in body
    assert conv.html_path.name in body
    assert "DIA-NN" in body
