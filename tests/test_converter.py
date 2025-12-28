"""Tests for Converter and ConverterBuilder."""

from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pandas as pd
import pytest

from anndata_proteomics.core import Converter
from anndata_proteomics.builder import ConverterBuilder
from anndata_proteomics.strategies.maxquant import MaxQuantStrategy
from anndata_proteomics.strategies.diann import DIANNStrategy


class SimpleTestStrategy:
    """Simple strategy for testing - mimics the strategy interface."""

    name = "TestSoftware"
    obs_id = "Run"
    var_id = "precursor_id"
    VAR_COLUMNS = ["Sequence", "Charge"]
    LAYER_COLUMNS = ["Intensity"]

    def __init__(self, layer_columns=None):
        if layer_columns:
            self.LAYER_COLUMNS = layer_columns

    def load(self, path):
        df = pd.read_csv(path, sep="\t")
        df["precursor_id"] = df["Sequence"].astype(str) + "/" + df["Charge"].astype(str)
        return df

    def get_obs(self, df):
        return pd.DataFrame(index=df[self.obs_id].unique())

    def get_var(self, df):
        available = [c for c in self.VAR_COLUMNS if c in df.columns]
        var_df = df.groupby(self.var_id).first()[available]
        return var_df

    def get_layers(self, df):
        cols = [self.obs_id, self.var_id] + [c for c in self.LAYER_COLUMNS if c in df.columns]
        return df[cols].copy()


class TestConverter:
    """Test Converter class."""

    @pytest.fixture
    def mock_strategy(self):
        """Create a mock strategy."""
        strategy = MagicMock()
        strategy.name = "MockSoftware"
        return strategy

    def test_init(self, mock_strategy):
        converter = Converter(mock_strategy)
        assert converter.strategy == mock_strategy

    def test_repr(self, mock_strategy):
        converter = Converter(mock_strategy)
        assert "Converter" in repr(converter)

    def test_convert_with_real_strategy(self, tmp_path):
        """Test conversion with a real strategy and simple data."""
        data_file = tmp_path / "data.tsv"
        data_file.write_text(
            "Sequence\tCharge\tIntensity\tRun\n"
            "PEPTIDE\t2\t1000\tSample1\n"
            "SEQUENCE\t2\t2000\tSample1\n"
            "PEPTIDE\t2\t1500\tSample2\n"
            "SEQUENCE\t2\t2500\tSample2\n"
        )

        ann_file = tmp_path / "annotation.csv"
        ann_file.write_text("sample,condition\nSample1,A\nSample2,B\n")

        strategy = SimpleTestStrategy()
        converter = Converter(strategy)
        adata = converter.convert(data_file, ann_file, log2_transform=False)

        assert adata.shape == (2, 2)  # 2 samples × 2 precursors
        assert set(adata.obs_names) == {"Sample1", "Sample2"}
        assert "Intensity" in adata.layers

    def test_convert_with_log2(self, tmp_path):
        """Test log2 transformation."""
        data_file = tmp_path / "data.tsv"
        data_file.write_text(
            "Sequence\tCharge\tIntensity\tRun\n"
            "PEP\t2\t1000\tS1\n"
            "PEP\t2\t2000\tS2\n"
        )

        ann_file = tmp_path / "annotation.csv"
        ann_file.write_text("sample,condition\nS1,A\nS2,B\n")

        strategy = SimpleTestStrategy()
        converter = Converter(strategy)
        adata = converter.convert(data_file, ann_file, log2_transform=True)

        # Check log2 transformation applied
        raw = adata.layers["Intensity"]
        transformed = adata.X
        np.testing.assert_array_almost_equal(transformed, np.log2(raw + 1))

    def test_x_layer_selection(self, tmp_path):
        """Test selecting which layer to use for X."""
        data_file = tmp_path / "data.tsv"
        data_file.write_text(
            "Sequence\tCharge\tIntensity\tScore\tRun\n"
            "PEP\t2\t1000\t0.95\tS1\n"
            "PEP\t2\t2000\t0.99\tS2\n"
        )

        ann_file = tmp_path / "annotation.csv"
        ann_file.write_text("sample,condition\nS1,A\nS2,B\n")

        strategy = SimpleTestStrategy(layer_columns=["Intensity", "Score"])
        converter = Converter(strategy)

        # Default: first layer (Intensity)
        adata1 = converter.convert(data_file, ann_file, log2_transform=False)
        assert adata1.uns["exploreDE"]["primary_quantification"] == "Intensity"

        # Explicit: use Score as X
        adata2 = converter.convert(data_file, ann_file, x_layer="Score", log2_transform=False)
        assert adata2.uns["exploreDE"]["primary_quantification"] == "Score"


class TestConverterBuilder:
    """Test ConverterBuilder class."""

    def test_list_supported(self):
        supported = ConverterBuilder.list_supported()
        assert "diann" in supported
        assert "maxquant" in supported
        assert "spectronaut" in supported

    def test_for_software_valid(self):
        converter = ConverterBuilder.for_software("MaxQuant")
        assert isinstance(converter, Converter)
        assert isinstance(converter.strategy, MaxQuantStrategy)

    def test_for_software_invalid(self):
        with pytest.raises(ValueError, match="Unknown software"):
            ConverterBuilder.for_software("NotARealSoftware")

    def test_for_software_case_insensitive(self):
        converter1 = ConverterBuilder.for_software("MaxQuant")
        converter2 = ConverterBuilder.for_software("maxquant")
        assert type(converter1.strategy) == type(converter2.strategy)

    def test_from_file_diann(self, tmp_path):
        # Create a DIA-NN-like file
        data_file = tmp_path / "report.tsv"
        data_file.write_text(
            "Modified.Sequence\tRun\tPrecursor.Quantity\tPrecursor.Id\n"
            "PEPTIDEK\tSample1\t1000\tPEPTIDEK/2\n"
        )

        converter = ConverterBuilder.from_file(data_file)
        assert isinstance(converter.strategy, DIANNStrategy)

    def test_from_file_maxquant(self, tmp_path):
        # Create a MaxQuant-like file
        data_file = tmp_path / "evidence.txt"
        data_file.write_text(
            "Modified sequence\tRaw file\tIntensity\tProteins\n"
            "_PEPTIDEK_\tSample1\t1000\tP12345\n"
        )

        converter = ConverterBuilder.from_file(data_file)
        assert isinstance(converter.strategy, MaxQuantStrategy)

    def test_from_file_unknown(self, tmp_path):
        data_file = tmp_path / "unknown.txt"
        data_file.write_text("RandomColumn\tAnotherColumn\nval1\tval2\n")

        with pytest.raises(ValueError, match="Could not detect"):
            ConverterBuilder.from_file(data_file)


class TestConverterIntegration:
    """Integration tests using real test data (skipped by default)."""

    @pytest.fixture
    def diann_file(self):
        path = Path("/Users/wolski/projects/ProteoBench/test/data/quant/quant_lfq_ion_DIA_AIF/DIANN_1.9_beta_sample_report.tsv")
        if not path.exists():
            pytest.skip("DIA-NN test file not found")
        return path

    @pytest.fixture
    def diann_annotation(self):
        path = Path("/Users/wolski/projects/anndata_proteomics_bridge/examples/diann_annotation.csv")
        if not path.exists():
            pytest.skip("DIA-NN annotation file not found")
        return path

    @pytest.mark.integration
    def test_diann_full_conversion(self, diann_file, diann_annotation):
        converter = ConverterBuilder.from_file(diann_file)
        adata = converter.convert(diann_file, diann_annotation)

        assert adata.shape[0] == 6  # 6 samples
        assert adata.shape[1] > 0  # Some precursors
        # Layer names use original column names (Precursor.Quantity for DIA-NN)
        assert "Precursor.Quantity" in adata.layers
        assert adata.uns["exploreDE"]["software"] == "DIA-NN"
