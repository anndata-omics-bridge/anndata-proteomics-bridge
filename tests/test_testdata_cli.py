import io
import zipfile
from pathlib import Path

import pytest

from anndata_proteomics.scripts.extract_raw_file_db import (
    ANNOTATION_URLS,
    GENERATED_CSV_NAMES,
    PROTEOBENCH_SETTINGS_REVISION,
    _download_module_jsons,
    _feature_count,
    _SubmissionMetadata,
    annotations,
    clean,
    clean_generated_data,
    fasta,
)


def test_feature_count_uses_current_proteobench_field() -> None:
    assert _feature_count(_SubmissionMetadata(nr_feature=42)) == 42


def test_feature_count_supports_legacy_proteobench_field() -> None:
    assert _feature_count(_SubmissionMetadata(nr_prec=17)) == 17


def test_clean_generated_data_removes_only_known_artifacts(tmp_path: Path) -> None:
    for directory in ("json_dir", "fasta", "annotations"):
        path = tmp_path / directory
        path.mkdir()
        (path / "generated").touch()
    for name in GENERATED_CSV_NAMES:
        (tmp_path / name).touch()
    keep = tmp_path / "keep.txt"
    keep.touch()

    removed = clean_generated_data(tmp_path)

    assert len(removed) == 6
    assert keep.exists()
    assert not any((tmp_path / name).exists() for name in GENERATED_CSV_NAMES)
    assert not (tmp_path / "json_dir").exists()
    assert not (tmp_path / "fasta").exists()
    assert not (tmp_path / "annotations").exists()


def test_clean_generated_data_is_idempotent(tmp_path: Path) -> None:
    assert clean_generated_data(tmp_path) == []


def test_clean_generated_data_rejects_filesystem_root() -> None:
    filesystem_root = Path(Path.cwd().anchor)

    with pytest.raises(ValueError, match="Refusing to clean"):
        clean_generated_data(filesystem_root)


def test_clean_command_uses_explicit_data_directory(tmp_path: Path) -> None:
    custom_data = tmp_path / "custom-test-data"
    custom_data.mkdir()
    (custom_data / GENERATED_CSV_NAMES[0]).touch()

    clean(data_dir=custom_data)

    assert not (custom_data / GENERATED_CSV_NAMES[0]).exists()


def test_download_module_jsons_normalizes_renamed_github_archive(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    archive = io.BytesIO()
    with zipfile.ZipFile(archive, "w") as zip_file:
        zip_file.writestr("renamed-repository-main/abc.json", "{}")

    class Response:
        content = archive.getvalue()

        @staticmethod
        def raise_for_status() -> None:
            return None

    monkeypatch.setattr(
        "anndata_proteomics.scripts.extract_raw_file_db.requests.get",
        lambda _url: Response(),
    )

    result = _download_module_jsons(
        "https://github.com/Proteobench/historical-name/archive/refs/heads/main.zip",
        tmp_path,
    )

    assert result == tmp_path / "historical-name" / "historical-name-main"
    assert (result / "abc.json").exists()


def test_download_module_jsons_replaces_only_metadata_snapshot(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    repo_dir = tmp_path / "historical-name"
    previous_metadata = repo_dir / "historical-name-main"
    downloaded_fixture = repo_dir / "fixture-hash"
    previous_metadata.mkdir(parents=True)
    downloaded_fixture.mkdir()
    (previous_metadata / "stale.json").write_text("{}")
    (downloaded_fixture / "input_file.tsv").write_text("quantity\n1\n")

    archive = io.BytesIO()
    with zipfile.ZipFile(archive, "w") as zip_file:
        zip_file.writestr("renamed-repository-main/current.json", "{}")

    class Response:
        content = archive.getvalue()

        @staticmethod
        def raise_for_status() -> None:
            return None

    monkeypatch.setattr(
        "anndata_proteomics.scripts.extract_raw_file_db.requests.get",
        lambda _url: Response(),
    )

    result = _download_module_jsons(
        "https://github.com/Proteobench/historical-name/archive/refs/heads/main.zip",
        tmp_path,
    )

    assert result == previous_metadata
    assert (result / "current.json").exists()
    assert not (result / "stale.json").exists()
    assert (downloaded_fixture / "input_file.tsv").read_text() == "quantity\n1\n"


def test_fasta_download_extracts_archives_without_zip_files(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    archive = io.BytesIO()
    with zipfile.ZipFile(archive, "w") as zip_file:
        zip_file.writestr("reference.fasta", ">protein\nPEPTIDE\n")

    class Response:
        content = archive.getvalue()

        @staticmethod
        def raise_for_status() -> None:
            return None

    monkeypatch.setattr(
        "anndata_proteomics.scripts.extract_raw_file_db.requests.get",
        lambda _url: Response(),
    )

    fasta(fasta_dir=tmp_path)

    assert (tmp_path / "reference.fasta").read_text() == ">protein\nPEPTIDE\n"
    assert not list(tmp_path.glob("*.zip"))


def test_annotation_download_fetches_and_validates_each_module_toml(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    annotation_content = b"""
[species_expected_ratio.HUMAN]
A_vs_B = 1.0
[species_mapper]
_HUMAN = "HUMAN"
[general]
min_count_multispec = 1
level = "ion"
[[samples]]
raw_file = "run1"
sample_name = "A1"
condition = "A"
[[samples]]
raw_file = "run2"
sample_name = "B1"
condition = "B"
"""

    class Response:
        def __init__(self, content: bytes) -> None:
            self.content = content

        @staticmethod
        def raise_for_status() -> None:
            return None

    requested: list[str] = []

    def get(url: str) -> Response:
        requested.append(url)
        return Response(annotation_content)

    monkeypatch.setattr(
        "anndata_proteomics.scripts.extract_raw_file_db.requests.get",
        get,
    )

    annotations(annotation_dir=tmp_path)

    assert requested == [*ANNOTATION_URLS.values()]
    assert {path.name for path in tmp_path.glob("*.toml")} == {
        f"{module}.toml" for module in ANNOTATION_URLS
    }
    assert not list(tmp_path.glob(".*.download.toml"))


def test_settings_urls_use_intermediate_format_contract_revision() -> None:
    assert PROTEOBENCH_SETTINGS_REVISION == "2738c47f8d621f0ee1fa4a6d3d358846f2bfa261"
    assert all(PROTEOBENCH_SETTINGS_REVISION in url for url in ANNOTATION_URLS.values())
