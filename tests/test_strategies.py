"""Tests for converter strategies."""

from pathlib import Path

import pandas as pd
import pytest

from anndata_proteomics.strategies.maxquant import MaxQuantStrategy
from anndata_proteomics.strategies.diann import DIANNStrategy
from anndata_proteomics.strategies.spectronaut import SpectronautStrategy


class TestMaxQuantStrategy:
    """Test MaxQuantStrategy."""

    def test_obs_var_ids(self):
        strategy = MaxQuantStrategy()
        assert strategy.obs_id == "Raw file"
        assert strategy.var_id == "precursor_id"

    def test_detection_columns(self):
        strategy = MaxQuantStrategy()
        assert "Modified sequence" in strategy.DETECTION_COLUMNS
        assert "Raw file" in strategy.DETECTION_COLUMNS

    def test_detect_positive(self, tmp_path):
        data_file = tmp_path / "evidence.txt"
        data_file.write_text(
            "Modified sequence\tRaw file\tIntensity\tProteins\tCharge\n"
            "_PEPTIDEK_\tSample1\t1000\tP12345\t2\n"
        )

        strategy = MaxQuantStrategy()
        assert strategy.detect(data_file)

    def test_detect_negative(self, tmp_path):
        data_file = tmp_path / "other.tsv"
        data_file.write_text("Column1\tColumn2\nval1\tval2\n")

        strategy = MaxQuantStrategy()
        assert not strategy.detect(data_file)

    def test_load_creates_precursor_id(self, tmp_path):
        data_file = tmp_path / "evidence.txt"
        data_file.write_text(
            "Modified sequence\tRaw file\tIntensity\tProteins\tCharge\n"
            "_PEPTIDEK_\tSample1\t1000\tP12345\t2\n"
        )

        strategy = MaxQuantStrategy()
        df = strategy.load(data_file)
        assert "precursor_id" in df.columns
        assert df["precursor_id"].iloc[0] == "_PEPTIDEK_/2"

    def test_var_columns(self):
        strategy = MaxQuantStrategy()
        assert "Modified sequence" in strategy.VAR_COLUMNS
        assert "Charge" in strategy.VAR_COLUMNS

    def test_layer_columns(self):
        strategy = MaxQuantStrategy()
        assert "Intensity" in strategy.LAYER_COLUMNS
        assert strategy.LAYER_COLUMNS[0] == "Intensity"  # First is default X


class TestDIANNStrategy:
    """Test DIANNStrategy."""

    def test_obs_var_ids(self):
        strategy = DIANNStrategy()
        assert strategy.obs_id == "Run"
        assert strategy.var_id == "Precursor.Id"

    def test_detection_columns(self):
        strategy = DIANNStrategy()
        assert "Modified.Sequence" in strategy.DETECTION_COLUMNS
        assert "Precursor.Quantity" in strategy.DETECTION_COLUMNS

    def test_detect_tsv(self, tmp_path):
        data_file = tmp_path / "report.tsv"
        data_file.write_text(
            "Modified.Sequence\tRun\tPrecursor.Quantity\tPrecursor.Id\n"
            "PEPTIDEK\tSample1\t1000\tPEPTIDEK/2\n"
        )

        strategy = DIANNStrategy()
        assert strategy.detect(data_file)

    def test_detect_negative(self, tmp_path):
        data_file = tmp_path / "other.tsv"
        data_file.write_text("Column1\tColumn2\nval1\tval2\n")

        strategy = DIANNStrategy()
        assert not strategy.detect(data_file)

    def test_var_columns(self):
        strategy = DIANNStrategy()
        assert "Modified.Sequence" in strategy.VAR_COLUMNS
        assert "Precursor.Charge" in strategy.VAR_COLUMNS
        assert "Protein.Names" in strategy.VAR_COLUMNS

    def test_layer_columns(self):
        strategy = DIANNStrategy()
        assert "Precursor.Quantity" in strategy.LAYER_COLUMNS
        assert strategy.LAYER_COLUMNS[0] == "Precursor.Quantity"  # First is default X
        assert "Q.Value" in strategy.LAYER_COLUMNS  # Q-values are layers too

    def test_get_var_with_real_data(self, tmp_path):
        data_file = tmp_path / "report.tsv"
        data_file.write_text(
            "Modified.Sequence\tRun\tPrecursor.Quantity\tPrecursor.Id\tPrecursor.Charge\tProtein.Names\n"
            "PEPTIDEK\tSample1\t1000\tPEPTIDEK/2\t2\tMyProtein\n"
            "PEPTIDEK\tSample2\t1500\tPEPTIDEK/2\t2\tMyProtein\n"
        )

        strategy = DIANNStrategy()
        df = strategy.load(data_file)
        var_df = strategy.get_var(df)

        assert len(var_df) == 1  # One unique precursor
        assert "Modified.Sequence" in var_df.columns
        assert "proforma" in var_df.columns


class TestSpectronautStrategy:
    """Test SpectronautStrategy."""

    def test_obs_var_ids(self):
        strategy = SpectronautStrategy()
        assert strategy.obs_id == "R.FileName"
        assert strategy.var_id == "precursor_id"

    def test_detection_columns(self):
        strategy = SpectronautStrategy()
        assert "FG.LabeledSequence" in strategy.DETECTION_COLUMNS
        assert "FG.Quantity" in strategy.DETECTION_COLUMNS

    def test_detect_positive(self, tmp_path):
        data_file = tmp_path / "spectronaut.tsv"
        data_file.write_text(
            "FG.LabeledSequence\tR.FileName\tFG.Quantity\tFG.Charge\n"
            "_PEPTIDEK_\tSample1\t1000\t2\n"
        )

        strategy = SpectronautStrategy()
        assert strategy.detect(data_file)

    def test_detect_negative(self, tmp_path):
        data_file = tmp_path / "other.tsv"
        data_file.write_text("Column1\tColumn2\nval1\tval2\n")

        strategy = SpectronautStrategy()
        assert not strategy.detect(data_file)

    def test_var_columns(self):
        strategy = SpectronautStrategy()
        assert "FG.LabeledSequence" in strategy.VAR_COLUMNS
        assert "FG.Charge" in strategy.VAR_COLUMNS

    def test_layer_columns(self):
        strategy = SpectronautStrategy()
        assert "FG.Quantity" in strategy.LAYER_COLUMNS
        assert "EG.Qvalue" in strategy.LAYER_COLUMNS

    def test_load_strips_underscores(self, tmp_path):
        data_file = tmp_path / "spectronaut.tsv"
        data_file.write_text(
            "FG.LabeledSequence\tR.FileName\tFG.Quantity\tFG.Charge\n"
            "_PEPTIDEK_\tSample1\t1000\t2\n"
        )

        strategy = SpectronautStrategy()
        df = strategy.load(data_file)
        assert df["FG.LabeledSequence"].iloc[0] == "PEPTIDEK"

    def test_load_creates_precursor_id(self, tmp_path):
        data_file = tmp_path / "spectronaut.tsv"
        data_file.write_text(
            "FG.LabeledSequence\tR.FileName\tFG.Quantity\tFG.Charge\n"
            "_PEPTIDEK_\tSample1\t1000\t2\n"
        )

        strategy = SpectronautStrategy()
        df = strategy.load(data_file)
        assert "precursor_id" in df.columns
        assert df["precursor_id"].iloc[0] == "PEPTIDEK/2"
