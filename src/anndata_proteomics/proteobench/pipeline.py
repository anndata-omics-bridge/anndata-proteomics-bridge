"""APB orchestration for ProteoBench quantitative scoring."""

from __future__ import annotations

from anndata import AnnData
from mudata import MuData

from anndata_proteomics.proteobench.config import ModuleSettings
from anndata_proteomics.proteobench.intermediate import align_runs, compute_intermediate
from anndata_proteomics.proteobench.metrics import (
    PROTEOBENCH_COMPATIBILITY_VERSION,
    PROTEOBENCH_SOURCE_REVISION,
    build_scores,
)
from anndata_proteomics.proteobench.resolve import resolve_roles, resolve_targets

_STORAGE_SCHEMA_VERSION = "0.1"


def score_quantification(
    obj: AnnData | MuData,
    module_settings: ModuleSettings,
) -> AnnData | MuData:
    """Compute and store ProteoBench HYE intermediates and scores.

    Scores are per quantification level, so a MuData is scored modality by modality and each
    modality keeps its own scores; nothing is written to the container's own ``uns``.

    Args:
        obj: Annotated APB AnnData or MuData object.
        module_settings: ProteoBench module experiment/scoring settings.

    Returns:
        The same container, enriched in every quantitative modality it holds.

    Raises:
        ValueError: Input/configuration contracts fail or scoring output exists.
    """
    for target in resolve_targets(obj):
        _score_target(target, module_settings)
    return obj


def _score_target(
    target: AnnData,
    module_settings: ModuleSettings,
) -> None:
    """Compute and store one level's intermediates and scores in place."""
    apb_namespace = dict(target.uns.get("anndata_proteomics") or {})
    namespace = dict(apb_namespace.get("proteobench") or {})
    if "proteobench" in target.varm:
        raise ValueError("varm['proteobench'] already exists; refusing to overwrite scores")
    if "scores" in namespace:
        raise ValueError(
            "uns['anndata_proteomics']['proteobench']['scores'] already exists; "
            "refusing to overwrite scores"
        )

    rule, roles = resolve_roles(target)
    design = align_runs(target, module_settings)
    intermediate = compute_intermediate(
        target,
        module_settings,
        roles,
        design,
        level=rule.quantification_level,
    )
    scores = build_scores(
        intermediate.legacy,
        intermediate.intermediate_hash,
        default_cutoff=module_settings.general.default_cutoff_min_feature,
        max_nr_observed=module_settings.general.max_nr_observed,
    )

    target.varm["proteobench"] = intermediate.varm
    namespace.update(
        {
            "schema_version": _STORAGE_SCHEMA_VERSION,
            "compatibility_version": PROTEOBENCH_COMPATIBILITY_VERSION,
            "source_revision": PROTEOBENCH_SOURCE_REVISION,
            "column_roles": roles.as_dict(),
            "protein_mapping": intermediate.protein_mapping,
            "scores": scores,
        }
    )
    apb_namespace["proteobench"] = namespace
    target.uns["anndata_proteomics"] = apb_namespace
