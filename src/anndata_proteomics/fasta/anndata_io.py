"""Round-trip resolved FASTA identifier configuration through ``uns``."""

from __future__ import annotations

from anndata import AnnData
from mudata import MuData

from anndata_proteomics.fasta.config import ResolvedFastaConfig

_NAMESPACE = "anndata_proteomics"
_CONFIG_KEY = "fasta_config"


def read_fasta_config(obj: AnnData | MuData) -> ResolvedFastaConfig | None:
    """Return the stored resolved FASTA configuration, when present."""
    namespace = obj.uns.get(_NAMESPACE)
    if not namespace:
        return None
    payload = namespace.get(_CONFIG_KEY)
    if payload is None:
        return None
    if isinstance(payload, bytes):
        payload = payload.decode("utf-8")
    return ResolvedFastaConfig.model_validate_json(str(payload))


def write_fasta_config(obj: AnnData | MuData, config: ResolvedFastaConfig) -> None:
    """Store the canonical live FASTA configuration on an AnnData or MuData."""
    namespace = dict(obj.uns.get(_NAMESPACE, {}))
    namespace[_CONFIG_KEY] = config.model_dump_json()
    obj.uns[_NAMESPACE] = namespace
