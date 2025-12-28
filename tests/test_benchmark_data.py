"""Tests for BenchmarkDataLoader class."""

import json
import zipfile
from io import BytesIO
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from anndata_proteomics.test_data import BENCHMARK_REPOS, BenchmarkDataLoader


class TestBenchmarkDataLoaderInit:
    """Test class initialization."""

    def test_valid_benchmark(self):
        loader = BenchmarkDataLoader("dda_qexactive")
        assert loader.benchmark == "dda_qexactive"
        assert loader.config == BENCHMARK_REPOS["dda_qexactive"]

    def test_invalid_benchmark_raises(self):
        with pytest.raises(ValueError, match="Unknown benchmark"):
            BenchmarkDataLoader("invalid_benchmark")

    def test_output_dir_created(self, tmp_path):
        output_dir = tmp_path / "test_output"
        loader = BenchmarkDataLoader("dda_qexactive", output_dir=str(output_dir))
        assert output_dir.exists()

    def test_properties(self):
        loader = BenchmarkDataLoader("dda_qexactive")
        assert "github.com" in loader.repo_url
        assert "proteobench" in loader.data_url
        assert loader.repo_name == "Results_quant_ion_DDA"


class TestBenchmarkDataLoaderMetadata:
    """Test metadata loading."""

    @pytest.fixture
    def mock_zip_response(self):
        """Create a mock zip file with JSON metadata."""
        # Create in-memory zip with JSON files
        zip_buffer = BytesIO()
        with zipfile.ZipFile(zip_buffer, "w") as zf:
            # Add JSON files matching expected structure
            json_data = {
                "software_name": "MaxQuant",
                "software_version": "2.1.3.0",
                "intermediate_hash": "abc123def456",
                "search_engine": "Andromeda",
            }
            zf.writestr(
                "Results_quant_ion_DDA-main/abc123def456.json",
                json.dumps(json_data),
            )
            json_data2 = {
                "software_name": "FragPipe",
                "software_version": "22.0",
                "intermediate_hash": "xyz789",
                "search_engine": "MSFragger",
            }
            zf.writestr(
                "Results_quant_ion_DDA-main/xyz789.json",
                json.dumps(json_data2),
            )
        zip_buffer.seek(0)
        return zip_buffer.getvalue()

    def test_load_metadata(self, tmp_path, mock_zip_response):
        loader = BenchmarkDataLoader("dda_qexactive", output_dir=str(tmp_path))

        with patch("requests.get") as mock_get:
            mock_response = MagicMock()
            mock_response.content = mock_zip_response
            mock_response.raise_for_status = MagicMock()
            mock_get.return_value = mock_response

            df = loader.load_metadata()

        assert len(df) == 2
        assert "software_name" in df.columns
        assert set(df["software_name"]) == {"MaxQuant", "FragPipe"}

    def test_metadata_property_lazy_loads(self, tmp_path, mock_zip_response):
        loader = BenchmarkDataLoader("dda_qexactive", output_dir=str(tmp_path))
        assert loader._metadata is None

        with patch("requests.get") as mock_get:
            mock_response = MagicMock()
            mock_response.content = mock_zip_response
            mock_response.raise_for_status = MagicMock()
            mock_get.return_value = mock_response

            # Access property triggers load
            _ = loader.metadata

        assert loader._metadata is not None

    def test_list_software(self, tmp_path, mock_zip_response):
        loader = BenchmarkDataLoader("dda_qexactive", output_dir=str(tmp_path))

        with patch("requests.get") as mock_get:
            mock_response = MagicMock()
            mock_response.content = mock_zip_response
            mock_response.raise_for_status = MagicMock()
            mock_get.return_value = mock_response

            software = loader.list_software()

        assert "MaxQuant" in software
        assert "FragPipe" in software
        assert "2.1.3.0" in software["MaxQuant"]


class TestBenchmarkDataLoaderDownload:
    """Test dataset downloading."""

    @pytest.fixture
    def loader_with_metadata(self, tmp_path):
        """Create loader with pre-loaded metadata."""
        loader = BenchmarkDataLoader("dda_qexactive", output_dir=str(tmp_path))
        loader._metadata = pd.DataFrame([
            {"software_name": "MaxQuant", "software_version": "2.1.3.0", "intermediate_hash": "hash123"},
            {"software_name": "FragPipe", "software_version": "22.0", "intermediate_hash": "hash456"},
        ])
        return loader

    def test_fetch_available_folders(self, loader_with_metadata):
        html_content = """
        <html><body>
        <a href="hash123/">hash123/</a>
        <a href="hash456/">hash456/</a>
        <a href="other/">other/</a>
        </body></html>
        """
        with patch("requests.get") as mock_get:
            mock_response = MagicMock()
            mock_response.text = html_content
            mock_response.raise_for_status = MagicMock()
            mock_get.return_value = mock_response

            folders = loader_with_metadata._fetch_available_folders()

        assert "hash123" in folders
        assert "hash456" in folders

    def test_find_zip_files(self, loader_with_metadata):
        html_content = """
        <html><body>
        <a href="data.zip">data.zip</a>
        <a href="params.zip">params.zip</a>
        <a href="readme.txt">readme.txt</a>
        </body></html>
        """
        with patch("requests.get") as mock_get:
            mock_response = MagicMock()
            mock_response.text = html_content
            mock_response.raise_for_status = MagicMock()
            mock_get.return_value = mock_response

            zips = loader_with_metadata._find_zip_files("http://example.com/folder/")

        assert len(zips) == 2
        assert any("data.zip" in z for z in zips)

    def test_download_and_extract_zip(self, loader_with_metadata, tmp_path):
        # Create mock zip with a test file
        zip_buffer = BytesIO()
        with zipfile.ZipFile(zip_buffer, "w") as zf:
            zf.writestr("input_file.txt", "test content")
        zip_buffer.seek(0)

        with patch("requests.get") as mock_get:
            mock_response = MagicMock()
            mock_response.content = zip_buffer.getvalue()
            mock_response.raise_for_status = MagicMock()
            mock_get.return_value = mock_response

            extract_dir = tmp_path / "extracted"
            loader_with_metadata._download_and_extract_zip(
                "http://example.com/data.zip", extract_dir
            )

        assert extract_dir.exists()
        assert (extract_dir / "input_file.txt").exists()

    def test_get_input_file_path_exists(self, tmp_path):
        # Create directory with input file
        data_dir = tmp_path / "hash123"
        data_dir.mkdir()
        (data_dir / "input_file.txt").write_text("test")

        result = BenchmarkDataLoader._get_input_file_path(str(data_dir))
        assert result is not None
        assert "input_file.txt" in result

    def test_get_input_file_path_not_exists(self, tmp_path):
        data_dir = tmp_path / "empty"
        data_dir.mkdir()

        result = BenchmarkDataLoader._get_input_file_path(str(data_dir))
        assert result is None

    def test_get_input_file_path_empty_string(self):
        result = BenchmarkDataLoader._get_input_file_path("")
        assert result is None


class TestBenchmarkDataLoaderIntegration:
    """Integration tests (skipped by default, run with --run-integration)."""

    @pytest.fixture
    def loader(self, tmp_path):
        return BenchmarkDataLoader("dda_qexactive", output_dir=str(tmp_path))

    @pytest.mark.integration
    def test_load_real_metadata(self, loader):
        """Test loading real metadata from GitHub."""
        df = loader.load_metadata()
        assert len(df) > 0
        assert "software_name" in df.columns
        assert "intermediate_hash" in df.columns

    @pytest.mark.integration
    def test_list_real_software(self, loader):
        """Test listing real software from benchmark."""
        software = loader.list_software()
        assert len(software) > 0
        # Known software in DDA benchmark
        assert any(sw in software for sw in ["MaxQuant", "FragPipe", "AlphaPept"])


