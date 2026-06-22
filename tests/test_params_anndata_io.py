"""Tests for typed read/write of search parameters in ``AnnData.uns``."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import anndata as ad

from anndata_proteomics.params.anndata_io import (
    get_search_parameters_path,
    read_search_parameters,
    write_search_parameters,
)
from anndata_proteomics.params.model import MassTolerance, Parameters, Probability


def _empty_adata() -> ad.AnnData:
    return ad.AnnData(
        X=np.zeros((1, 1)),
        obs=pd.DataFrame(index=["o0"]),
        var=pd.DataFrame(index=["v0"]),
    )


def test_read_none_when_no_uns_entry():
    adata = _empty_adata()
    assert read_search_parameters(adata) is None
    assert get_search_parameters_path(adata) is None


def test_write_then_read_roundtrip():
    adata = _empty_adata()
    params = Parameters(
        software_name="Sage",
        software_version="0.14.6",
        enzyme="KR",
        allowed_miscleavages=1,
        fixed_mods="{'C': 57.02146}",
    )
    write_search_parameters(adata, params, source_path="/tmp/fake.json")

    recovered = read_search_parameters(adata)
    assert isinstance(recovered, Parameters)
    assert recovered.software_name == "Sage"
    assert recovered.software_version == "0.14.6"
    assert recovered.fixed_mods[0].source == "{'C': 57.02146}"
    assert recovered.to_series()["fixed_mods"] == "{'C': 57.02146}"
    assert get_search_parameters_path(adata) == "/tmp/fake.json"


def test_write_then_read_typed_value_roundtrip():
    adata = _empty_adata()
    params = Parameters(
        software_name="DIA-NN",
        software_version="2.3.0 Academia ",
        ident_fdr_psm=Probability(value=0.01),
        precursor_mass_tolerance=MassTolerance(
            mode="absolute", value=15.0, unit="ppm"
        ),
        fragment_mass_tolerance=MassTolerance(
            mode="absolute", value=20.0, unit="ppm"
        ),
    )
    write_search_parameters(adata, params, source_path="/tmp/diann.log.txt")

    recovered = read_search_parameters(adata)

    assert recovered.ident_fdr_psm == Probability(value=0.01)
    assert recovered.precursor_mass_tolerance == MassTolerance(
        mode="absolute", value=15.0, unit="ppm"
    )
    assert recovered.fragment_mass_tolerance == MassTolerance(
        mode="absolute", value=20.0, unit="ppm"
    )
    assert get_search_parameters_path(adata) == "/tmp/diann.log.txt"


def test_read_rejects_extra_fields():
    adata = _empty_adata()
    adata.uns["anndata_proteomics"] = {
        "search_parameters": '{"software_name": "Sage", "vendor_specific": 42}',
    }
    with pytest.raises(Exception):
        read_search_parameters(adata)


def test_read_validates_against_current_schema():
    adata = _empty_adata()
    adata.uns["anndata_proteomics"] = {
        "search_parameters": '{"software_name": "Sage"}',
    }
    recovered = read_search_parameters(adata)
    assert recovered.software_name == "Sage"


def test_read_raises_on_corrupt_payload():
    adata = _empty_adata()
    adata.uns["anndata_proteomics"] = {
        "search_parameters": "not-valid-json",
    }
    with pytest.raises(Exception):
        read_search_parameters(adata)
