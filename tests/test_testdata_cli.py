import io
import zipfile
from pathlib import Path

import pytest

from anndata_proteomics.scripts.extract_raw_file_db import (
    GENERATED_CSV_NAMES,
    _download_module_jsons,
    _feature_count,
    clean,
    clean_generated_data,
    fasta,
)


def test_feature_count_uses_current_proteobench_field() -> None:
    assert _feature_count({"nr_feature": 42}) == 42


def test_feature_count_supports_legacy_proteobench_field() -> None:
    assert _feature_count({"nr_prec": 17}) == 17


def test_clean_generated_data_removes_only_known_artifacts(tmp_path: Path) -> None:
    for directory in ("json_dir", "fasta"):
        path = tmp_path / directory
        path.mkdir()
        (path / "generated").touch()
    for name in GENERATED_CSV_NAMES:
        (tmp_path / name).touch()
    keep = tmp_path / "keep.txt"
    keep.touch()

    removed = clean_generated_data(tmp_path)

    assert len(removed) == 5
    assert keep.exists()
    assert not any((tmp_path / name).exists() for name in GENERATED_CSV_NAMES)
    assert not (tmp_path / "json_dir").exists()
    assert not (tmp_path / "fasta").exists()


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
    tmp_path: Path, monkeypatch
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


def test_fasta_download_extracts_archives_without_zip_files(tmp_path: Path, monkeypatch) -> None:
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
