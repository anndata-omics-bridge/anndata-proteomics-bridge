"""FragPipe rule regressions for auxiliary-layer sample normalization."""

from __future__ import annotations

import pandas as pd

from anndata_proteomics.converters.assemble import convert
from anndata_proteomics.rules.loader import load_packaged_rule


def test_fragpipe_rule_aligns_auxiliary_group_suffixes_to_intensity_runs() -> None:
    frame = pd.DataFrame(
        {
            "Peptide Sequence": ["PEPTIDE"],
            "Modified Sequence": ["PEPTIDE"],
            "Charge": [2],
            "M/Z": [400.0],
            "Protein": ["P1"],
            "Mapped Proteins": [""],
            "Protein ID": ["P1"],
            "Gene": ["G1"],
            "Assigned Modifications": [""],
            "LFQ_Run_01 Intensity": [100.0],
            "LFQ_Run_01_1 Spectral Count": [2],
            "LFQ_Run_01_1 Apex Retention Time": [1.5],
            "LFQ_Run_01_1 Match Type": ["MS/MS"],
        }
    )

    result = convert(frame, load_packaged_rule("fragpipe", "ion"))

    assert result.obs_names.tolist() == ["LFQ_Run_01"]
    assert set(result.layers) == {
        "Intensity",
        "Spectral_Count",
        "Apex_Retention_Time",
        "Match_Type",
    }
