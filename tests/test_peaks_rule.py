"""PEAKS rule regressions for run-level wide-column recognition."""

from __future__ import annotations

import pandas as pd

from anndata_proteomics.converters.assemble import convert
from anndata_proteomics.rules.loader import load_packaged_rule


def test_peaks_rule_excludes_summary_columns_and_normalizes_raw_suffixes() -> None:
    frame = pd.DataFrame(
        {
            "Peptide": ["PEPTIDE"],
            "Quality": [99.0],
            "Significance": [1.0],
            "m/z": [400.0],
            "RT range": ["1 - 2"],
            "z": [2],
            "Avg. Area": [10.0],
            "Accession": ["P1"],
            "LFQ_Run_1 Normalized Area": [100.0],
            "LFQ_Run_1.raw m/z": [400.0],
            "LFQ_Run_1_raw RT mean": [1.5],
            "LFQ_Run_1.raw AScore": [42.0],
            "Group 1 Normalized Area": [999.0],
            "Condition A Normalized Area": [999.0],
            "Best AScore": [999.0],
        }
    )

    result = convert(frame, load_packaged_rule("peaks", "ion"))

    assert result.obs_names.tolist() == ["LFQ_Run_1"]
    assert set(result.layers) == {
        "Normalized_Area",
        "Sample_Mz",
        "Sample_RT_Mean",
        "AScore",
    }
