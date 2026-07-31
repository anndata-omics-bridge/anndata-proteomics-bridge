"""Load a converted APB result (``.h5ad`` AnnData or ``.h5mu`` MuData) back into memory."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import mudata
from anndata import AnnData
from mudata import MuData


def load_converted_result(result_path: Path | str) -> AnnData | MuData:
    """Load a converted ``result.h5ad`` or ``result.h5mu`` file."""
    path = Path(result_path).expanduser()
    if path.suffix == ".h5ad":
        return ad.read_h5ad(path)
    if path.suffix == ".h5mu":
        # Adopt the mudata 0.4 default now (no auto-pull of per-modality obs/var into the
        # global frames); modalities keep their own obs/var. Silences the 0.3 FutureWarning.
        with mudata.set_options(pull_on_update=False):
            return mudata.read_h5mu(path)
    raise ValueError(f"unsupported converted result type: {path}")
