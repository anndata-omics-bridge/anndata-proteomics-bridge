"""Tests for inferred and persisted FASTA identifier configuration."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import mudata
import numpy as np
import pytest
from pydantic import ValidationError

from anndata_proteomics.adapters.anndata.fasta import (
    require_fasta_config,
    write_anndata_fasta_config,
    write_mudata_fasta_config,
)
from anndata_proteomics.fasta.config import (
    ExplicitPatterns,
    FastaConfig,
    resolve_fasta_config,
)


def test_inference_records_matching_candidates_and_counts() -> None:
    resolved = resolve_fasta_config(
        ("sp|P1|A", "REV_sp|P1|A", "rev_P2", "zz|C1|CONT"),
        FastaConfig(),
    )
    assert resolved.decoy.patterns == (r"^REV_", r"^rev_")
    assert resolved.decoy.match_counts[r"^REV_"] == 1
    assert resolved.decoy.match_counts[r"^rev_"] == 1
    assert resolved.contaminant.patterns == (r"^zz(?:\||_)",)
    assert resolved.n_fasta_ids == 4


def test_explicit_patterns_override_candidates() -> None:
    resolved = resolve_fasta_config(
        ("REV_P1", "CUSTOM_P2"),
        FastaConfig(decoy=ExplicitPatterns(patterns=(r"^CUSTOM_",))),
    )
    assert resolved.decoy.source == "explicit"
    assert resolved.decoy.patterns == (r"^CUSTOM_",)
    assert resolved.decoy.match_counts == {r"^CUSTOM_": 1}


def test_explicit_empty_tuple_disables_classification() -> None:
    resolved = resolve_fasta_config(
        ("REV_P1", "zz|C1|CONT"),
        FastaConfig(
            decoy=ExplicitPatterns(patterns=()),
            contaminant=ExplicitPatterns(patterns=()),
        ),
    )
    assert resolved.decoy.source == "explicit"
    assert resolved.decoy.patterns == ()
    assert resolved.contaminant.patterns == ()


def test_no_candidate_hit_resolves_to_none() -> None:
    resolved = resolve_fasta_config(("sp|P1|A",), FastaConfig())
    assert resolved.decoy.source == "none"
    assert resolved.decoy.patterns == ()
    assert resolved.contaminant.source == "none"


def test_invalid_regex_is_rejected() -> None:
    with pytest.raises(ValidationError, match="invalid FASTA identifier regex"):
        FastaConfig(decoy=ExplicitPatterns(patterns=("[",)))


def test_config_round_trips_through_h5ad_and_h5mu(tmp_path: Path) -> None:
    resolved = resolve_fasta_config(("REV_P1",), FastaConfig())
    adata = ad.AnnData(np.zeros((1, 1)))
    write_anndata_fasta_config(adata, resolved)
    h5ad_path = tmp_path / "config.h5ad"
    adata.write_h5ad(h5ad_path)
    assert require_fasta_config(ad.read_h5ad(h5ad_path)) == resolved

    with mudata.set_options(pull_on_update=False):
        mdata = mudata.MuData({"ion": adata}, axis=0)
    write_mudata_fasta_config(mdata, resolved)
    h5mu_path = tmp_path / "config.h5mu"
    mdata.write_h5mu(h5mu_path)
    with mudata.set_options(pull_on_update=False):
        restored = mudata.read_h5mu(h5mu_path)
    assert require_fasta_config(restored) == resolved
