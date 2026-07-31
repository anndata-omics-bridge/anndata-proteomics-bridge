"""Tests for table-driven observation annotation."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import mudata
import numpy as np
import pandas as pd
import pytest
from mudata import MuData

from anndata_proteomics.annotation.apply import annotate_obs
from anndata_proteomics.annotation.loader import AnnotationTable, load_annotation
from anndata_proteomics.readers.summary import describe

RUNS = ["runA1", "runA2", "runB1", "runB2"]

_BASIC_TOML = """
schema_version = "0.1"

[obs]
match_on = "index"
key_field = "raw_file"

[[obs.samples]]
raw_file = "runA1"
sample_name = "A_rep1"
condition = "A"

[[obs.samples]]
raw_file = "runA2"
sample_name = "A_rep2"
condition = "A"

[[obs.samples]]
raw_file = "runB1"
sample_name = "B_rep1"
condition = "B"

[[obs.samples]]
raw_file = "runB2"
sample_name = "B_rep2"
condition = "B"
"""


def _adata(
    var_prefix: str = "ion:",
    n_var: int = 3,
    runs: list[str] = RUNS,
) -> ad.AnnData:
    var_names = [f"{var_prefix}{index}" for index in range(n_var)]
    return ad.AnnData(
        X=np.arange(len(runs) * n_var, dtype="float64").reshape(len(runs), n_var),
        obs=pd.DataFrame(index=pd.Index(list(runs), name="R_FileName")),
        var=pd.DataFrame(index=pd.Index(var_names, name="ProForma_ion")),
    )


def _mudata() -> MuData:
    mods = {"ion": _adata("ion:", 3), "protein": _adata("prt:", 2)}
    with mudata.set_options(pull_on_update=False):
        return MuData(mods, axis=0)


def _annotation_from(tmp_path: Path, text: str = _BASIC_TOML) -> AnnotationTable:
    path = tmp_path / "annotation.toml"
    path.write_text(text)
    return load_annotation(path)


def test_obs_join_by_index(tmp_path: Path) -> None:
    adata = _adata()
    annotate_obs(adata, _annotation_from(tmp_path))
    assert list(adata.obs["condition"]) == ["A", "A", "B", "B"]
    assert list(adata.obs["sample_name"]) == ["A_rep1", "A_rep2", "B_rep1", "B_rep2"]


def test_obs_join_falls_back_to_exact_raw_file_alias(tmp_path: Path) -> None:
    text = _BASIC_TOML.replace(
        'raw_file = "runA1"',
        'raw_file = "runA1"\nraw_file_alias = "aliasA1"',
    )
    text = text.replace(
        'raw_file = "runA2"',
        'raw_file = "runA2"\nraw_file_alias = "aliasA2"',
    )
    text = text.replace(
        'raw_file = "runB1"',
        'raw_file = "runB1"\nraw_file_alias = "aliasB1"',
    )
    text = text.replace(
        'raw_file = "runB2"',
        'raw_file = "runB2"\nraw_file_alias = "aliasB2"',
    )
    adata = _adata(runs=["aliasA1", "aliasA2", "aliasB1", "aliasB2"])

    annotate_obs(adata, _annotation_from(tmp_path, text))

    assert list(adata.obs["condition"]) == ["A", "A", "B", "B"]
    assert list(adata.obs["sample_name"]) == ["A_rep1", "A_rep2", "B_rep1", "B_rep2"]
    provenance = adata.uns["anndata_proteomics"]["obs_annotations_json"]
    assert '"key_field": "raw_file_alias"' in provenance


def test_obs_join_supports_multiple_exact_raw_file_aliases(tmp_path: Path) -> None:
    text = _BASIC_TOML.replace(
        'raw_file = "runA1"',
        'raw_file = "runA1"\nraw_file_aliases = ["aliasA1", "secondA1"]',
    )
    text = text.replace(
        'raw_file = "runA2"',
        'raw_file = "runA2"\nraw_file_aliases = ["aliasA2", "secondA2"]',
    )
    text = text.replace(
        'raw_file = "runB1"',
        'raw_file = "runB1"\nraw_file_aliases = ["aliasB1", "secondB1"]',
    )
    text = text.replace(
        'raw_file = "runB2"',
        'raw_file = "runB2"\nraw_file_aliases = ["aliasB2", "secondB2"]',
    )
    adata = _adata(runs=["secondA1", "secondA2", "secondB1", "secondB2"])

    annotate_obs(adata, _annotation_from(tmp_path, text))

    assert list(adata.obs["condition"]) == ["A", "A", "B", "B"]
    assert list(adata.obs["sample_name"]) == ["A_rep1", "A_rep2", "B_rep1", "B_rep2"]
    provenance = adata.uns["anndata_proteomics"]["obs_annotations_json"]
    assert '"key_field": "raw_file_aliases"' in provenance


def test_duplicate_raw_file_aliases_raise(tmp_path: Path) -> None:
    text = _BASIC_TOML.replace(
        'raw_file = "runA1"',
        'raw_file = "runA1"\nraw_file_aliases = ["shared", "aliasA1"]',
    )
    text = text.replace(
        'raw_file = "runA2"',
        'raw_file = "runA2"\nraw_file_aliases = ["shared", "aliasA2"]',
    )

    with pytest.raises(ValueError, match="duplicate 'raw_file_aliases'"):
        annotate_obs(_adata(runs=["shared"]), _annotation_from(tmp_path, text))


def test_obs_join_does_not_fall_back_to_sample_name(tmp_path: Path) -> None:
    adata = _adata(runs=["A_rep1", "A_rep2", "B_rep1", "B_rep2"])

    with pytest.raises(ValueError, match="no obs rows matched"):
        annotate_obs(adata, _annotation_from(tmp_path))


def test_annotation_preserves_quantification_view_and_adds_provenance(tmp_path: Path) -> None:
    adata = _adata()
    matrix = adata.X
    assert isinstance(matrix, np.ndarray)
    adata.layers["intensity"] = matrix.copy()
    adata.uns["anndata_proteomics"] = {
        "quantification_level": "ion",
        "software_name": "Synthetic",
    }
    before = describe(adata)["quantification"]

    annotate_obs(adata, _annotation_from(tmp_path))

    result = describe(adata)
    assert result["quantification"] == before
    assert result["annotations"]["obs"] == [
        {
            "source": str(tmp_path / "annotation.toml"),
            "source_format": "toml",
            "match_on": "index",
            "key_field": "raw_file",
            "obs_columns_added": ["sample_name", "condition"],
            "n_obs_matched": 4,
        }
    ]


def test_join_respects_obs_order(tmp_path: Path) -> None:
    adata = _adata(runs=["runB2", "runA1", "runB1", "runA2"])
    annotate_obs(adata, _annotation_from(tmp_path))
    assert list(adata.obs["condition"]) == ["B", "A", "B", "A"]


def test_match_on_named_column(tmp_path: Path) -> None:
    adata = _adata()
    adata.obs_names = ["x0", "x1", "x2", "x3"]
    adata.obs["Run"] = RUNS
    annotation = _annotation_from(
        tmp_path,
        _BASIC_TOML.replace('match_on = "index"', 'match_on = "Run"'),
    )
    annotate_obs(adata, annotation)
    assert list(adata.obs["condition"]) == ["A", "A", "B", "B"]


def test_freeform_extra_columns(tmp_path: Path) -> None:
    text = _BASIC_TOML.replace('sample_name = "A_rep1"', 'batch = 1\ngenotype = "wt"')
    text = text.replace('sample_name = "A_rep2"', 'batch = 2\ngenotype = "ko"')
    text = text.replace('sample_name = "B_rep1"', 'batch = 1\ngenotype = "wt"')
    text = text.replace('sample_name = "B_rep2"', 'batch = 2\ngenotype = "ko"')
    adata = _adata()
    annotate_obs(adata, _annotation_from(tmp_path, text))
    assert list(adata.obs["batch"]) == [1, 2, 1, 2]
    assert list(adata.obs["genotype"]) == ["wt", "ko", "wt", "ko"]


def test_mudata_annotates_global_and_modalities(tmp_path: Path) -> None:
    md = _mudata()
    annotate_obs(md, _annotation_from(tmp_path))
    for frame in (md.obs, md.mod["ion"].obs, md.mod["protein"].obs):
        assert list(frame["condition"]) == ["A", "A", "B", "B"]


def test_mudata_roundtrip(tmp_path: Path) -> None:
    md = _mudata()
    annotate_obs(md, _annotation_from(tmp_path))
    output = tmp_path / "md.annotated.h5mu"
    md.write_h5mu(output)
    with mudata.set_options(pull_on_update=False):
        roundtrip = mudata.read_h5mu(output)
    assert list(roundtrip.obs["condition"]) == ["A", "A", "B", "B"]
    assert list(roundtrip.mod["ion"].obs["condition"]) == ["A", "A", "B", "B"]


def test_anndata_roundtrip_records_provenance(tmp_path: Path) -> None:
    adata = _adata()
    annotate_obs(adata, _annotation_from(tmp_path))
    output = tmp_path / "a.annotated.h5ad"
    adata.write_h5ad(output)
    roundtrip = ad.read_h5ad(output)
    assert list(roundtrip.obs["condition"]) == ["A", "A", "B", "B"]
    assert "obs_annotations_json" in roundtrip.uns["anndata_proteomics"]
    provenance = describe(roundtrip)["annotations"]["obs"]
    assert provenance[0]["obs_columns_added"] == ["sample_name", "condition"]
    assert provenance[0]["n_obs_matched"] == 4


def test_obs_column_names_sanitised(tmp_path: Path) -> None:
    text = _BASIC_TOML.replace("sample_name", '"Sample Name"')
    adata = _adata()
    annotate_obs(adata, _annotation_from(tmp_path, text))
    assert "Sample_Name" in adata.obs.columns
    assert "Sample Name" not in adata.obs.columns


def test_sanitisation_collision_raises(tmp_path: Path) -> None:
    text = _BASIC_TOML.replace(
        'sample_name = "A_rep1"',
        '"Sample Name" = "x"\n"Sample-Name" = "y"',
    )
    adata = _adata()
    with pytest.raises(ValueError, match="collision after sanitisation"):
        annotate_obs(adata, _annotation_from(tmp_path, text))


def test_no_match_raises(tmp_path: Path) -> None:
    adata = _adata(runs=["nope1", "nope2", "nope3", "nope4"])
    with pytest.raises(ValueError, match="no obs rows matched"):
        annotate_obs(adata, _annotation_from(tmp_path))


def test_partial_match_warns(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    adata = _adata(runs=["runA1", "runA2", "runB1", "extra_run"])
    annotate_obs(adata, _annotation_from(tmp_path))
    error = capsys.readouterr().err
    assert "1/4 obs rows had no matching" in error
    assert "annotation record(s) matched no obs row" in error
    assert list(adata.obs["condition"])[:3] == ["A", "A", "B"]
    assert pd.isna(adata.obs["condition"].iloc[3])


def test_collision_with_existing_obs_column_raises(tmp_path: Path) -> None:
    adata = _adata()
    adata.obs["condition"] = ["pre", "pre", "pre", "pre"]
    with pytest.raises(ValueError, match="already present in obs"):
        annotate_obs(adata, _annotation_from(tmp_path))


def test_duplicate_key_field_raises(tmp_path: Path) -> None:
    duplicate = """
[[obs.samples]]
raw_file = "runA1"
sample_name = "duplicate"
condition = "A"
"""
    with pytest.raises(ValueError, match="duplicate 'raw_file'"):
        annotate_obs(_adata(), _annotation_from(tmp_path, _BASIC_TOML + duplicate))


def test_unknown_match_on_column_raises(tmp_path: Path) -> None:
    text = _BASIC_TOML.replace('match_on = "index"', 'match_on = "NoSuchColumn"')
    with pytest.raises(ValueError, match="match_on column 'NoSuchColumn' not found"):
        annotate_obs(_adata(), _annotation_from(tmp_path, text))


def test_loads_proteobench_module_settings_without_modelling_extra_sections(
    tmp_path: Path,
) -> None:
    path = tmp_path / "module_settings.toml"
    path.write_text(
        """
[general]
level = "ion"

[species_expected_ratio.HUMAN]
A_vs_B = 1.0

[[samples]]
raw_file = "runA1"
sample_name = "A_rep1"
condition = "A"
"""
    )

    annotation = load_annotation(path)

    assert annotation.match_on == "index"
    assert annotation.key_field == "raw_file"
    assert annotation.samples.to_dict(orient="records") == [
        {"raw_file": "runA1", "sample_name": "A_rep1", "condition": "A"}
    ]


@pytest.mark.parametrize(("suffix", "separator"), [(".csv", ","), (".tsv", "\t")])
def test_loads_delimited_annotation_tables(
    tmp_path: Path,
    suffix: str,
    separator: str,
) -> None:
    path = tmp_path / f"samples{suffix}"
    path.write_text(f"raw_file{separator}condition\nrunA1{separator}A\n")

    annotation = load_annotation(path)

    assert annotation.samples.to_dict(orient="records") == [{"raw_file": "runA1", "condition": "A"}]


def test_loader_rejects_missing_key_field(tmp_path: Path) -> None:
    path = tmp_path / "bad.tsv"
    path.write_text("sample_name\tcondition\nA_rep1\tA\n")
    with pytest.raises(ValueError, match="missing key field 'raw_file'"):
        load_annotation(path)


def test_loader_rejects_json(tmp_path: Path) -> None:
    path = tmp_path / "annotation.json"
    path.write_text("{}")
    with pytest.raises(ValueError, match="Unsupported annotation format"):
        load_annotation(path)


def test_loader_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        load_annotation(tmp_path / "does_not_exist.toml")


def test_committed_aif_fixture_applies() -> None:
    fixture = Path(__file__).parent / "data" / "annotation" / "obs_samples_AIF.toml"
    annotation = load_annotation(fixture)
    runs = [
        "LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_01",
        "LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_02",
        "LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_03",
        "LFQ_Orbitrap_AIF_Condition_B_Sample_Alpha_01",
        "LFQ_Orbitrap_AIF_Condition_B_Sample_Alpha_02",
        "LFQ_Orbitrap_AIF_Condition_B_Sample_Alpha_03",
    ]
    adata = _adata(runs=runs)
    annotate_obs(adata, annotation)
    assert list(adata.obs["condition"]) == ["A", "A", "A", "B", "B", "B"]
