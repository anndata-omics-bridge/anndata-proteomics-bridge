"""Backend-independent orchestration for ProteoBench scoring."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from anndata_proteomics.proteobench.config import ModuleSettings
from anndata_proteomics.proteobench.contracts import QuantMatrix
from anndata_proteomics.proteobench.pipeline import ProteoBenchResult, score_level
from anndata_proteomics.rules.schema import QuantificationLevel


@dataclass(frozen=True, slots=True)
class ProteoBenchLevelInput:
    """Exact values required to score one quantification level."""

    observations: pd.DataFrame
    matrix: QuantMatrix
    feature_ids: pd.Index
    reported_proteins: pd.Series
    level: QuantificationLevel


def score_levels(
    levels: tuple[ProteoBenchLevelInput, ...],
    module_settings: ModuleSettings,
) -> tuple[ProteoBenchResult, ...]:
    """Score every extracted level in input order."""
    return tuple(
        score_level(
            level.observations,
            level.matrix,
            level.feature_ids,
            level.reported_proteins,
            module_settings,
            level.level,
        )
        for level in levels
    )
