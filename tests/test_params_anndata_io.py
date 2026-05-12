"""Tests for typed read/write of search parameters in ``AnnData.uns``."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import anndata as ad
from pydantic import ValidationError

from anndata_proteomics.params.anndata_io import (
    get_search_parameters_path,
    read_search_parameters,
    write_search_parameters,
)
from anndata_proteomics.params.model import Parameters


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
    assert recovered.fixed_mods == "{'C': 57.02146}"
    assert get_search_parameters_path(adata) == "/tmp/fake.json"


def test_write_preserves_extra_fields():
    adata = _empty_adata()
    params = Parameters(software_name="Sage", vendor_specific=42)
    write_search_parameters(adata, params)
    recovered = read_search_parameters(adata)
    assert recovered.model_dump()["vendor_specific"] == 42


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
