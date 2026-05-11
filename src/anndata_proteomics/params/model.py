"""Parameter model for proteomics search-engine settings."""

from __future__ import annotations

from typing import Any

import pandas as pd
from pydantic import BaseModel, ConfigDict


class Parameters(BaseModel):
    """Proteomics search-parameter record.

    Field set mirrors ProteoBench's ``ProteoBenchParameters`` JSON
    definitions (``proteobench/io/params/json/Quant/*.json``). All fields
    are optional; values stay loosely typed so that vendor-specific
    encodings (e.g. ``fixed_mods`` as a stringified dict) round-trip
    unchanged for ProteoBench-equivalence tests.
    """

    model_config = ConfigDict(extra="allow")

    software_name: Any | None = None
    software_version: Any | None = None
    search_engine: Any | None = None
    search_engine_version: Any | None = None
    ident_fdr_psm: Any | None = None
    ident_fdr_peptide: Any | None = None
    ident_fdr_protein: Any | None = None
    enable_match_between_runs: Any | None = None
    precursor_mass_tolerance: Any | None = None
    fragment_mass_tolerance: Any | None = None
    enzyme: Any | None = None
    semi_enzymatic: Any | None = None
    allowed_miscleavages: Any | None = None
    min_peptide_length: Any | None = None
    max_peptide_length: Any | None = None
    fixed_mods: Any | None = None
    variable_mods: Any | None = None
    max_mods: Any | None = None
    min_precursor_charge: Any | None = None
    max_precursor_charge: Any | None = None
    min_precursor_mz: Any | None = None
    max_precursor_mz: Any | None = None
    min_fragment_mz: Any | None = None
    max_fragment_mz: Any | None = None
    quantification_method: Any | None = None
    protein_inference: Any | None = None
    abundance_normalization_ions: Any | None = None
    predictors_library: Any | None = None
    scan_window: Any | None = None

    def to_series(self) -> pd.Series:
        """Serialize to a pandas Series matching ProteoBench's CSV layout."""
        return pd.Series(self.model_dump())

    @classmethod
    def from_series(cls, series: pd.Series) -> "Parameters":
        """Build a ``Parameters`` instance from a ProteoBench-style Series.

        Empty strings become ``None``; literal ``"None"`` strings are kept
        as-is because ProteoBench writes them that way.
        """
        data: dict[str, Any] = {}
        for key, value in series.items():
            if isinstance(value, float) and pd.isna(value):
                data[str(key)] = None
            elif value == "":
                data[str(key)] = None
            else:
                data[str(key)] = value
        return cls(**data)
