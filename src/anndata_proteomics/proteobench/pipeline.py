"""Container-independent orchestration for ProteoBench scoring."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from anndata_proteomics.proteobench.config import ModuleSettings
from anndata_proteomics.proteobench.contracts import QuantMatrix
from anndata_proteomics.proteobench.intermediate import (
    IntermediateResult,
    align_runs,
    compute_intermediate,
)
from anndata_proteomics.proteobench.metrics import (
    ProteoBenchScores,
    ScoreConfig,
    build_scores,
)
from anndata_proteomics.rules.schema import QuantificationLevel


@dataclass(frozen=True)
class ProteoBenchResult:
    """Complete calculated result for one quantification level."""

    intermediate: IntermediateResult
    scores: ProteoBenchScores


def score_level(
    observations: pd.DataFrame,
    matrix: QuantMatrix,
    feature_ids: pd.Index,
    reported_proteins: pd.Series,
    module_settings: ModuleSettings,
    level: QuantificationLevel,
) -> ProteoBenchResult:
    """Calculate ProteoBench intermediates and scores for one level."""
    design = align_runs(observations, module_settings)
    intermediate = compute_intermediate(
        matrix,
        feature_ids,
        reported_proteins,
        module_settings,
        design,
        level,
    )
    scores = build_scores(
        intermediate.legacy,
        intermediate.intermediate_hash,
        ScoreConfig(
            default_cutoff=module_settings.general.default_cutoff_min_feature,
            max_nr_observed=module_settings.general.max_nr_observed,
        ),
    )
    return ProteoBenchResult(intermediate=intermediate, scores=scores)
