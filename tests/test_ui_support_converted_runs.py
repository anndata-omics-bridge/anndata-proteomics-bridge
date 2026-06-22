from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from anndata_proteomics.params.anndata_io import write_search_parameters
from anndata_proteomics.params.model import MassTolerance, Parameters
from anndata_proteomics.scripts import _ui_support as ui


def _params() -> Parameters:
    return Parameters(
        software_name="DIA-NN",
        software_version="2.3.0 Academia ",
        search_engine_version="2.3.0 Academia ",
        enzyme="Trypsin/P",
        precursor_mass_tolerance=MassTolerance(value=15.0, unit="ppm", mode="absolute"),
        fragment_mass_tolerance=MassTolerance(value=20.0, unit="ppm", mode="absolute"),
        enable_match_between_runs=False,
    )


def test_list_converted_runs_empty_dir(tmp_path) -> None:
    runs = ui.list_converted_runs(tmp_path / "missing")

    assert list(runs.columns) == ui.CONVERTED_COLUMNS
    assert runs.empty


def test_list_converted_runs_discovers_h5ad_result(tmp_path) -> None:
    import anndata as ad

    run_dir = tmp_path / "20260622T120000_diann_ion"
    run_dir.mkdir()
    result = run_dir / "result.h5ad"
    ad.AnnData(
        X=np.array([[1.0, 2.0]]),
        obs=pd.DataFrame(index=["run1"]),
        var=pd.DataFrame(index=["feature1", "feature2"]),
    ).write_h5ad(result)
    (run_dir / "console.log").write_text("DONE\n")

    runs = ui.list_converted_runs(tmp_path)

    assert len(runs) == 1
    row = runs.iloc[0]
    assert row["run_name"] == "20260622T120000_diann_ion"
    assert row["timestamp"] == "20260622T120000"
    assert row["slug"] == "diann"
    assert row["target"] == "ion"
    assert row["status"] == "finished"
    assert row["result_type"] == "h5ad"
    assert row["result_path"] == str(result)
    loaded = ui.load_converted_result(row["result_path"])
    assert loaded.shape == (1, 2)


def test_list_converted_runs_joins_catalog_metadata(tmp_path, monkeypatch) -> None:
    import anndata as ad

    rel_input = "Results_quant_ion_DIA_AIF/example/input_file.txt"
    param_path = "/tmp/param_0..txt"
    catalog = pd.DataFrame(
        [
            {
                "software_name": "DIA-NN",
                "software_version": "1.9.2",
                "nr_prec": 123,
                "size_mb": 45.6,
                "input_file_path": rel_input,
                "param_path": param_path,
            }
        ]
    )
    monkeypatch.setattr(ui, "load_catalog", lambda: catalog)
    run_dir = tmp_path / "20260622T120003_diann_ion"
    run_dir.mkdir()
    ad.AnnData(
        X=np.array([[1.0]]),
        obs=pd.DataFrame(index=["run1"]),
        var=pd.DataFrame(index=["feature1"]),
    ).write_h5ad(run_dir / "result.h5ad")
    (run_dir / "console.log").write_text(
        "$ python -m anndata_proteomics.scripts.convert_one "
        f"--input {rel_input} --slug diann --target ion --params {param_path} "
        f"--outdir {run_dir}\nDONE\n"
    )

    runs = ui.list_converted_runs(tmp_path)
    row = runs.iloc[0]
    assert row["software_name"] == "DIA-NN"
    assert row["software_version"] == "1.9.2"
    assert row["nr_prec"] == "123"
    assert row["size_mb"] == "45.6"
    assert row["input_file_path"] == rel_input
    assert row["param_path"] == param_path

    display = ui.converted_runs_table(runs)
    assert list(display.columns) == [
        "run_name",
        "software_name",
        "software_version",
        "target",
        "status",
        "result_type",
        "nr_prec",
        "size_mb",
    ]
    assert "result_path" not in display.columns
    assert display.iloc[0]["software_name"] == "DIA-NN"


def test_finished_run_uses_artifact_parameters_without_catalog(tmp_path, monkeypatch) -> None:
    import anndata as ad

    monkeypatch.setattr(ui, "load_catalog", lambda: pd.DataFrame())
    run_dir = tmp_path / "20260622T120004_diann_ion"
    run_dir.mkdir()
    adata = ad.AnnData(
        X=np.array([[1.0]]),
        obs=pd.DataFrame(index=["run1"]),
        var=pd.DataFrame(index=["feature1"]),
    )
    write_search_parameters(adata, _params(), source_path="/tmp/diann.log.txt")
    adata.write_h5ad(run_dir / "result.h5ad")

    runs = ui.list_converted_runs(tmp_path)
    row = runs.iloc[0]

    assert row["software_name"] == "DIA-NN"
    assert row["software_version"] == "2.3.0 Academia "
    assert row["param_path"] == "/tmp/diann.log.txt"


def test_list_converted_runs_discovers_h5mu_result(tmp_path) -> None:
    from anndata import AnnData
    from mudata import MuData

    run_dir = tmp_path / "20260622T120001_diann_mudata"
    run_dir.mkdir()
    result = run_dir / "result.h5mu"
    MuData(
        {
            "ion": AnnData(
                X=np.array([[1.0]]),
                obs=pd.DataFrame(index=["run1"]),
                var=pd.DataFrame(index=["ion1"]),
            )
        },
        axis=0,
    ).write(result)

    runs = ui.list_converted_runs(tmp_path)

    assert len(runs) == 1
    row = runs.iloc[0]
    assert row["target"] == "mudata"
    assert row["status"] == "finished"
    assert row["result_type"] == "h5mu"
    loaded = ui.load_converted_result(row["result_path"])
    assert "ion" in loaded.mod


def test_summarize_shows_full_search_parameters() -> None:
    from anndata import AnnData

    adata = AnnData(
        X=np.array([[1.0]]),
        obs=pd.DataFrame(index=["run1"]),
        var=pd.DataFrame(index=["feature1"]),
    )
    write_search_parameters(adata, _params(), source_path="/tmp/diann.log.txt")

    summary = ui.summarize(adata)
    params = summary["search_parameters"]

    assert params["software_name"] == "DIA-NN"
    assert params["software_version"] == "2.3.0 Academia "
    assert params["enzyme"] == "Trypsin/P"
    assert params["precursor_mass_tolerance"] == {
        "value": 15.0,
        "unit": "ppm",
        "mode": "absolute",
    }
    assert params["fragment_mass_tolerance"] == {
        "value": 20.0,
        "unit": "ppm",
        "mode": "absolute",
    }
    assert params["search_parameters_path"] == "/tmp/diann.log.txt"


def test_summarize_mudata_shows_root_search_parameters() -> None:
    from anndata import AnnData
    from mudata import MuData

    mdata = MuData(
        {
            "ion": AnnData(
                X=np.array([[1.0]]),
                obs=pd.DataFrame(index=["run1"]),
                var=pd.DataFrame(index=["ion1"]),
            )
        },
        axis=0,
    )
    write_search_parameters(mdata, _params(), source_path="/tmp/diann.log.txt")

    summary = ui.summarize(mdata)

    assert summary["search_parameters"]["software_name"] == "DIA-NN"
    assert summary["search_parameters"]["software_version"] == "2.3.0 Academia "
    assert summary["search_parameters"]["search_parameters_path"] == "/tmp/diann.log.txt"
    assert "search_parameters" not in summary["modalities"]["ion"]


def test_summarize_mudata_falls_back_to_modality_search_parameters() -> None:
    from anndata import AnnData
    from mudata import MuData

    adata = AnnData(
        X=np.array([[1.0]]),
        obs=pd.DataFrame(index=["run1"]),
        var=pd.DataFrame(index=["ion1"]),
    )
    write_search_parameters(adata, _params(), source_path="/tmp/diann.log.txt")
    mdata = MuData({"ion": adata}, axis=0)

    summary = ui.summarize(mdata)

    assert summary["search_parameters"]["software_name"] == "DIA-NN"
    assert summary["search_parameters"]["source"] == "modalities"
    assert summary["search_parameters"]["modalities"] == ["ion"]
    assert "search_parameters" not in summary["modalities"]["ion"]


def test_list_converted_runs_keeps_log_only_run(tmp_path) -> None:
    run_dir = tmp_path / "20260622T120002_diann_mudata"
    run_dir.mkdir()
    log = run_dir / "console.log"
    log.write_text("Traceback\n")

    runs = ui.list_converted_runs(tmp_path)

    assert len(runs) == 1
    row = runs.iloc[0]
    assert row["status"] == "incomplete"
    assert row["result_path"] == ""
    assert row["log_path"] == str(log)


def test_load_converted_result_rejects_unsupported_suffix(tmp_path) -> None:
    result = tmp_path / "result.txt"
    result.write_text("not an AnnData result")

    with pytest.raises(ValueError, match="unsupported converted result type"):
        ui.load_converted_result(result)


def test_convert_level_passes_params_path(monkeypatch) -> None:
    from anndata import AnnData
    from anndata_proteomics.converters import assemble

    captured = {}

    def fake_convert(df, rule, *, params_path=None):
        captured["params_path"] = params_path
        return AnnData(
            X=np.array([[1.0]]),
            obs=pd.DataFrame(index=["run1"]),
            var=pd.DataFrame(index=["feature1"]),
        )

    monkeypatch.setattr(ui, "select_rule", lambda slug, level, version, headers: object())
    monkeypatch.setattr(assemble, "convert", fake_convert)

    adata = ui._convert_level(
        pd.DataFrame({"x": [1]}),
        "diann",
        "ion",
        "1.9.2",
        params_path="/tmp/param_0..txt",
    )

    assert adata.shape == (1, 1)
    assert captured["params_path"] == "/tmp/param_0..txt"
