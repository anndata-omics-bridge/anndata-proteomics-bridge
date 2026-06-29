"""Load a converted APB result (``.h5ad`` AnnData or ``.h5mu`` MuData) back into memory."""

from __future__ import annotations

from pathlib import Path
from typing import Any


def load_converted_result(result_path: Path | str) -> Any:
    """Load a converted ``result.h5ad`` or ``result.h5mu`` file."""
    path = Path(result_path).expanduser()
    if path.suffix == ".h5ad":
        import anndata as ad

        return ad.read_h5ad(path)
    if path.suffix == ".h5mu":
        import mudata

        # Adopt the mudata 0.4 default now (no auto-pull of per-modality obs/var into the
        # global frames); modalities keep their own obs/var. Silences the 0.3 FutureWarning.
        with mudata.set_options(pull_on_update=False):
            return mudata.read_h5mu(path)
    raise ValueError(f"unsupported converted result type: {path}")
