"""PEAKS rule regressions for run-level wide-column recognition."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd
from conversion_support import convert_to_anndata

from anndata_proteomics.adapters.anndata.matrix import layer_names
from anndata_proteomics.rules.loader import load_packaged_rule


def _frame(rows: int = 1, **overrides: Sequence[object]) -> pd.DataFrame:
    """A minimal PEAKS ion export; `overrides` replace single columns."""
    base: dict[str, Sequence[object]] = {
        "Peptide": ["PEPTIDE", "PEPTIDEK", "PEPTIDER"][:rows],
        "Quality": [99.0] * rows,
        "Significance": [1.0] * rows,
        "m/z": [400.0] * rows,
        "RT range": ["1 - 2"] * rows,
        "z": [2] * rows,
        "Avg. Area": [10.0] * rows,
        "Accession": ["P1"] * rows,
        "LFQ_Run_1 Normalized Area": [100.0] * rows,
        "LFQ_Run_1.raw m/z": [400.0] * rows,
        "LFQ_Run_1_raw RT mean": [1.5] * rows,
        "LFQ_Run_1.raw AScore": ["S3:Phosphorylation:1000.00"] * rows,
        "Group 1 Normalized Area": [999.0] * rows,
        "Condition A Normalized Area": [999.0] * rows,
        "Best AScore": [999.0] * rows,
    }
    return pd.DataFrame(base | overrides)


def test_peaks_rule_excludes_summary_columns_and_normalizes_raw_suffixes() -> None:
    result = convert_to_anndata(_frame(), load_packaged_rule("peaks", "ion"))

    assert result.obs_names.tolist() == ["LFQ_Run_1"]
    assert set(layer_names(result)) == {
        "Normalized_Area",
        "Sample_Mz",
        "Sample_RT_Mean",
        "AScore",
    }


def test_peaks_ascore_extracts_the_score_from_the_vendor_site_string() -> None:
    """PEAKS writes `site:modification:score`, not a bare number.

    Without `value_pattern` the whole column coerces to NaN and the layer is silently
    empty -- the defect this test pins.
    """
    frame = _frame(
        rows=3,
        **{
            "LFQ_Run_1.raw AScore": [
                "S3:Phosphorylation:76.54",
                # Multi-site cells are ';'-separated; the first score wins.
                "C1:Carbamidomethylation:1000.00;C1:Acetylation (Protein N-term):12.30",
                # Unmodified peptides carry no AScore at all.
                None,
            ]
        },
    )

    result = convert_to_anndata(frame, load_packaged_rule("peaks", "ion"))

    ascore = np.asarray(result.layers["AScore"], dtype="float64")
    assert ascore[0, 0] == 76.54
    assert ascore[0, 1] == 1000.00
    assert np.isnan(ascore[0, 2])


def test_peaks_zero_normalized_area_is_the_not_detected_sentinel() -> None:
    """PEAKS writes 0 for an ion it did not quantify; m/z carries '-' in the same cell."""
    frame = _frame(
        rows=2,
        **{
            "LFQ_Run_1 Normalized Area": [100.0, 0.0],
            "LFQ_Run_1.raw m/z": [400.0, "-"],
        },
    )

    result = convert_to_anndata(frame, load_packaged_rule("peaks", "ion"))

    area = np.asarray(result.layers["Normalized_Area"], dtype="float64")
    mz = np.asarray(result.layers["Sample_Mz"], dtype="float64")
    assert area[0, 0] == 100.0
    assert np.isnan(area[0, 1])
    assert np.isnan(mz[0, 1])
