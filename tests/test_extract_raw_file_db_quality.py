"""Behavior and failure-path coverage for the test-data cache CLI."""

from __future__ import annotations

import io
import json
import zipfile
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from bs4 import BeautifulSoup

from anndata_proteomics.scripts import extract_raw_file_db as rawdb


class _Response:
    def __init__(self, *, content: bytes = b"", text: str = "") -> None:
        self.content = content
        self.text = text
        self.status_checked = False

    def raise_for_status(self) -> None:
        self.status_checked = True

    def iter_content(self, _size: int) -> tuple[bytes, ...]:
        return (self.content,)


def _zip_bytes(files: dict[str, bytes]) -> bytes:
    stream = io.BytesIO()
    with zipfile.ZipFile(stream, "w") as archive:
        for name, content in files.items():
            archive.writestr(name, content)
    return stream.getvalue()


def test_feature_count_and_safe_zip_extraction(tmp_path: Path) -> None:
    assert rawdb._feature_count({"nr_feature": 3, "nr_prec": 2}) == 3
    assert rawdb._feature_count({"nr_feature": None, "nr_prec": 2}) == 2

    with zipfile.ZipFile(io.BytesIO(_zip_bytes({"root/file.txt": b"x"}))) as archive:
        rawdb._extract_zip(archive, tmp_path)
    assert (tmp_path / "root/file.txt").read_text(encoding="utf-8") == "x"

    unsafe = _zip_bytes({"../escape.txt": b"x"})
    with zipfile.ZipFile(io.BytesIO(unsafe)) as archive:
        with pytest.raises(RuntimeError, match="Unsafe ZIP member"):
            rawdb._extract_zip(archive, tmp_path)


def test_get_merged_json_discovers_redirected_root(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    response = _Response(
        content=_zip_bytes({"Renamed-main/a.json": b"{}", "Renamed-main/b.json": b"{}"})
    )
    monkeypatch.setattr(rawdb.requests, "get", lambda _url: response)
    monkeypatch.chdir(tmp_path)
    extracted = rawdb.get_merged_json("https://github.com/org/Original/archive/refs/heads/main.zip")
    assert extracted == Path("Original") / "Renamed-main"
    assert (tmp_path / extracted / "a.json").exists()
    assert response.status_checked

    response.content = _zip_bytes({"one/a": b"", "two/b": b""})
    with pytest.raises(RuntimeError, match="one root folder"):
        rawdb.get_merged_json("https://github.com/org/Original/archive/refs/heads/main.zip")


def test_href_filter_and_raw_download_idempotency(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    soup = BeautifulSoup(
        '<a href="hash/">hash</a><a href="archive.zip">zip</a><a>none</a>',
        "html.parser",
    )
    assert rawdb._hrefs_ending_with(soup, "/") == ["hash/"]
    assert rawdb._hrefs_ending_with(soup, ".zip") == ["archive.zip"]

    output = tmp_path / "downloads"
    existing = output / "present"
    existing.mkdir(parents=True)
    (existing / "input_file.tsv").write_text("old", encoding="utf-8")
    archive = _zip_bytes({"input_file.tsv": b"new"})
    responses = {
        "https://server/": _Response(
            text='<a href="present/">present</a><a href="fresh/">fresh</a>'
        ),
        "https://server/fresh/": _Response(text='<a href="data.zip">data</a>'),
        "https://server/fresh/data.zip": _Response(content=archive),
    }
    monkeypatch.setattr(
        rawdb.requests,
        "get",
        lambda url, **_kwargs: responses[url],
    )
    monkeypatch.chdir(tmp_path)
    found = rawdb.get_raw_data(
        pd.DataFrame({"intermediate_hash": ["present", "fresh"]}),
        base_url="https://server/",
        output_directory=str(output),
    )
    assert found["present"] == str(existing)
    assert (Path(found["fresh"]) / "input_file.tsv").read_text(encoding="utf-8") == "new"
    assert not (tmp_path / "data.zip").exists()


def test_metadata_refresh_preserves_and_restores_previous_snapshot(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cache = tmp_path / "cache"
    target = cache / "Repo" / "Repo-main"
    target.mkdir(parents=True)
    (target / "old.json").write_text("{}", encoding="utf-8")

    def extract(repo_url: str) -> Path:
        assert repo_url
        root = Path("Repo-main")
        root.mkdir()
        (root / "new.json").write_text("{}", encoding="utf-8")
        return root

    monkeypatch.setattr(rawdb, "get_merged_json", extract)
    refreshed = rawdb._download_module_jsons(
        "https://github.com/org/Repo/archive/refs/heads/main.zip",
        cache,
    )
    assert refreshed == target
    assert (target / "new.json").exists()
    assert not (target / "old.json").exists()

    def missing_extract(repo_url: str) -> Path:
        assert repo_url
        return Path("missing")

    monkeypatch.setattr(rawdb, "get_merged_json", missing_extract)
    with pytest.raises(RuntimeError, match="Expected extracted folder"):
        rawdb._download_module_jsons(
            "https://github.com/org/Repo/archive/refs/heads/main.zip",
            cache,
        )


def test_metadata_refresh_rolls_back_failed_replace(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cache = tmp_path / "cache"
    target = cache / "Repo" / "Repo-main"
    target.mkdir(parents=True)
    (target / "old.json").write_text("{}", encoding="utf-8")

    def extract(repo_url: str) -> Path:
        assert repo_url
        root = Path("Repo-main")
        root.mkdir()
        (root / "new.json").write_text("{}", encoding="utf-8")
        return root

    monkeypatch.setattr(rawdb, "get_merged_json", extract)
    original_replace = Path.replace

    def fail_new_snapshot(path: Path, destination: Path) -> Path:
        if path.name == "Repo-main" and "metadata-" in str(path.parent):
            raise OSError("replace failed")
        return original_replace(path, destination)

    monkeypatch.setattr(Path, "replace", fail_new_snapshot)
    with pytest.raises(OSError, match="replace failed"):
        rawdb._download_module_jsons(
            "https://github.com/org/Repo/archive/refs/heads/main.zip",
            cache,
        )
    assert (target / "old.json").exists()


def test_metadata_refresh_failure_without_previous_and_final_guard(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cache = tmp_path / "cache"

    def extract(repo_url: str) -> Path:
        assert repo_url
        root = Path("Repo-main")
        root.mkdir()
        return root

    monkeypatch.setattr(rawdb, "get_merged_json", extract)
    original_replace = Path.replace

    def fail_new_snapshot(path: Path, destination: Path) -> Path:
        if path.name == "Repo-main" and "metadata-" in str(path.parent):
            raise OSError("replace failed")
        return original_replace(path, destination)

    monkeypatch.setattr(Path, "replace", fail_new_snapshot)
    with pytest.raises(OSError, match="replace failed"):
        rawdb._download_module_jsons(
            "https://github.com/org/Repo/archive/refs/heads/main.zip",
            cache,
        )

    monkeypatch.setattr(Path, "replace", original_replace)
    original_is_dir = Path.is_dir
    target = cache / "Repo" / "Repo-main"

    def hide_final_target(path: Path) -> bool:
        return False if path == target else original_is_dir(path)

    monkeypatch.setattr(Path, "is_dir", hide_final_target)
    with pytest.raises(RuntimeError, match="Expected extracted folder"):
        rawdb._download_module_jsons(
            "https://github.com/org/Repo/archive/refs/heads/main.zip",
            cache,
        )


def test_catalog_and_selection_commands(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    metadata = tmp_path / "metadata"
    metadata.mkdir()
    (metadata / "good.json").write_text(
        json.dumps(
            {
                "intermediate_hash": "hash-b",
                "software_name": "Tool",
                "software_version": "1",
                "nr_prec": 5,
            }
        ),
        encoding="utf-8",
    )
    (metadata / "bad.json").write_text("{", encoding="utf-8")
    monkeypatch.setattr(
        rawdb,
        "CONFIGS",
        {"dda_qexactive": {"repo_url": "https://github.com/org/Repo/archive/refs/heads/main.zip"}},
    )
    monkeypatch.setattr(rawdb, "_download_module_jsons", lambda *_args: metadata)
    catalog_csv = tmp_path / "catalog.csv"
    rawdb.catalog(catalog_csv=catalog_csv, cache_dir=tmp_path / "cache")
    catalog = pd.read_csv(catalog_csv)
    assert catalog["intermediate_hash"].tolist() == ["hash-b"]
    assert catalog["nr_feature"].tolist() == [5]

    extra = pd.concat(
        [
            catalog,
            catalog.assign(intermediate_hash="hash-a", nr_feature=5),
            catalog.assign(intermediate_hash="missing", nr_feature=None),
        ],
        ignore_index=True,
    )
    extra.to_csv(catalog_csv, index=False)
    selection_csv = tmp_path / "selection.csv"
    rawdb.select(
        catalog_csv=catalog_csv,
        selection_csv=selection_csv,
        strategy="smallest-per-software-version",
        module="dda_qexactive",
    )
    selected = pd.read_csv(selection_csv)
    assert selected["intermediate_hash"].tolist() == ["hash-a"]
    rawdb.select(
        catalog_csv=catalog_csv,
        selection_csv=selection_csv,
        strategy="all",
    )
    assert len(pd.read_csv(selection_csv)) == 3


@pytest.mark.parametrize(
    ("strategy", "expected"),
    [
        ("all", 4),
        ("smallest-per-software-version", 2),
        ("smallest-per-software", 1),
        ("smallest-per-module", 1),
    ],
)
def test_selection_strategies(strategy: Any, expected: int) -> None:
    frame = pd.DataFrame(
        {
            "module": ["m"] * 4,
            "software_name": ["A", "A", "A", "B"],
            "software_version": ["1", "1", "2", "1"],
            "nr_feature": [2, 1, 3, None],
            "intermediate_hash": ["b", "a", "c", "d"],
        }
    )
    selected, missing, _rule = rawdb._select_rows(frame, strategy)
    assert len(selected) == expected
    assert missing == (0 if strategy == "all" else 1)


def test_existing_dataset_detection(tmp_path: Path) -> None:
    output = tmp_path / "repo"
    output.mkdir()
    (output / "present").mkdir()
    (output / "unrelated").mkdir()
    frame = pd.DataFrame({"intermediate_hash": ["present", "missing"]})
    remaining, present = rawdb.get_datasets_to_download(frame, output)
    assert remaining["intermediate_hash"].tolist() == ["missing"]
    assert present == {"present": str(output / "present")}
    untouched, empty = rawdb.get_datasets_to_download(frame, tmp_path / "absent")
    assert untouched is frame
    assert empty == {}


def test_download_manifest_statuses_and_validation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    invalid = tmp_path / "invalid.csv"
    pd.DataFrame({"module": ["m"]}).to_csv(invalid, index=False)
    with pytest.raises(SystemExit, match="missing required columns"):
        rawdb.download(
            selection_csv=invalid,
            cache_dir=tmp_path / "cache",
            manifest_csv=tmp_path / "manifest.csv",
        )

    selection = tmp_path / "selection.csv"
    rows = pd.DataFrame(
        {
            "module": ["m", "m", "m"],
            "repo_name": ["Repo", "Repo", "Repo"],
            "intermediate_hash": ["ok", "empty", "absent"],
        }
    )
    rows.to_csv(selection, index=False)
    cache = tmp_path / "cache"
    ok = cache / "Repo" / "ok"
    empty = cache / "Repo" / "empty"
    ok.mkdir(parents=True)
    empty.mkdir(parents=True)
    (ok / "input_file.tsv").write_text("x", encoding="utf-8")
    monkeypatch.setattr(rawdb, "get_raw_data", lambda *_args, **_kwargs: {})
    manifest = tmp_path / "manifest.csv"
    rawdb.download(selection_csv=selection, cache_dir=cache, manifest_csv=manifest)
    result = pd.read_csv(manifest).set_index("intermediate_hash")
    assert result.loc["ok", "status"] == "ok"
    assert result.loc["empty", "status"] == "input_file_missing"
    assert result.loc["absent", "status"] == "not_on_server"

    numeric = tmp_path / "numeric.csv"
    rows.assign(repo_name=1).to_csv(numeric, index=False)
    with pytest.raises(TypeError, match="repo_name"):
        rawdb.download(
            selection_csv=numeric,
            cache_dir=cache,
            manifest_csv=manifest,
        )

    present_only = tmp_path / "present.csv"
    rows.iloc[[0]].to_csv(present_only, index=False)
    monkeypatch.setattr(
        rawdb,
        "get_raw_data",
        lambda *_args, **_kwargs: pytest.fail("download should be skipped"),
    )
    rawdb.download(
        selection_csv=present_only,
        cache_dir=cache,
        manifest_csv=manifest,
    )


def test_fasta_annotations_and_cleanup(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    archive = _zip_bytes(
        {
            "reference.fasta": b">P1\nAAAA\n",
            "__MACOSX/metadata": b"x",
        }
    )
    monkeypatch.setattr(rawdb, "FASTA_URLS", ("https://server/fasta.zip",))
    monkeypatch.setattr(rawdb.requests, "get", lambda _url: _Response(content=archive))
    fasta_dir = tmp_path / "fasta"
    rawdb.fasta(fasta_dir=fasta_dir)
    assert (fasta_dir / "reference.fasta").exists()
    assert not (fasta_dir / "__MACOSX").exists()

    monkeypatch.setattr(rawdb, "ANNOTATION_URLS", {"dda": "https://server/module"})
    monkeypatch.setattr(
        rawdb.requests,
        "get",
        lambda _url: _Response(content=b"[general]\nlevel='ion'\n"),
    )
    monkeypatch.setattr(rawdb, "load_annotation", lambda _path: pd.DataFrame())
    monkeypatch.setattr(rawdb, "load_module_settings", lambda _path: object())
    annotation_dir = tmp_path / "annotations"
    rawdb.annotations(annotation_dir=annotation_dir)
    assert (annotation_dir / "dda.toml").exists()
    assert not list(annotation_dir.glob(".*.download.toml"))

    data_dir = tmp_path / "generated"
    for name in ("json_dir", "fasta", "annotations"):
        (data_dir / name).mkdir(parents=True, exist_ok=True)
    for name in rawdb.GENERATED_CSV_NAMES:
        (data_dir / name).write_text("x", encoding="utf-8")
    removed = rawdb.clean_generated_data(data_dir)
    assert len(removed) == 6
    rawdb.clean(data_dir=data_dir)
    with pytest.raises(ValueError, match="filesystem"):
        rawdb.clean_generated_data(Path("/"))


def test_annotations_remove_invalid_temporary_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(rawdb, "ANNOTATION_URLS", {"dda": "https://server/module"})
    monkeypatch.setattr(
        rawdb.requests,
        "get",
        lambda _url: _Response(content=b"invalid"),
    )
    monkeypatch.setattr(
        rawdb,
        "load_annotation",
        lambda _path: (_ for _ in ()).throw(ValueError("invalid")),
    )
    annotation_dir = tmp_path / "annotations"
    with pytest.raises(ValueError, match="invalid"):
        rawdb.annotations(annotation_dir=annotation_dir)
    assert not list(annotation_dir.glob(".*.download.toml"))


def test_build_database_and_main(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    (tmp_path / "plain.txt").write_text("skip", encoding="utf-8")
    without_main = tmp_path / "without_main"
    without_main.mkdir()
    module = tmp_path / "module"
    json_dir = module / "module-main"
    json_dir.mkdir(parents=True)
    missing_json = module / "missing-json"
    missing_json.mkdir()
    missing_input = module / "missing-input"
    missing_input.mkdir()
    ok = module / "ok"
    ok.mkdir()
    (json_dir / "missing-input.json").write_text(
        '{"software_name":"Tool","software_version":"1"}',
        encoding="utf-8",
    )
    (json_dir / "ok.json").write_text(
        '{"software_name":"Tool","software_version":"2"}',
        encoding="utf-8",
    )
    (ok / "input_file.txt").write_text("data", encoding="utf-8")
    result = rawdb.build_database(tmp_path).set_index("intermediate_hash")
    assert result.loc["missing-json", "status"] == "json_missing"
    assert result.loc["missing-input", "status"] == "input_file_missing"
    assert result.loc["ok", "status"] == "ok"

    called: list[bool] = []
    monkeypatch.setattr(rawdb, "app", lambda: called.append(True))
    rawdb.main()
    assert called == [True]
