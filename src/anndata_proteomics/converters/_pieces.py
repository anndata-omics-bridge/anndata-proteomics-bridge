"""Shared dataclass for converter outputs."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd


@dataclass
class ConversionPieces:
    """The pieces a converter produces, ready for assemble.to_anndata()."""

    X: np.ndarray
    obs: pd.DataFrame
    var: pd.DataFrame
    layers: dict[str, np.ndarray] = field(default_factory=dict)
    uns: dict[str, Any] = field(default_factory=dict)
