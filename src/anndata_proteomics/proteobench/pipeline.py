"""APB orchestration for ProteoBench quantitative scoring."""

from __future__ import annotations

from typing import Any

from anndata_proteomics.proteobench.config import ModuleSettings, ToolSettings
from anndata_proteomics.proteobench.intermediate import align_runs, compute_intermediate
from anndata_proteomics.proteobench.metrics import (
    PROTEOBENCH_COMPATIBILITY_VERSION,
    PROTEOBENCH_SOURCE_REVISION,
    build_scores,
)
from anndata_proteomics.proteobench.resolve import resolve_roles, resolve_target

_STORAGE_SCHEMA_VERSION = "0.1"


def score_quantification(
    obj: Any,
    module_settings: ModuleSettings,
    tool_settings: ToolSettings,
) -> Any:
    """Compute and store ProteoBench HYE intermediates and scores.

    Args:
        obj: Converted APB AnnData or MuData object.
        module_settings: ProteoBench module experiment/scoring settings.
        tool_settings: ProteoBench per-tool raw-column settings.

    Returns:
        The same container, enriched in its selected quantitative modality.

    Raises:
        ValueError: Input/configuration contracts fail or scoring output exists.
    """
    target = resolve_target(obj, module_settings.general.level)
    namespace = target.uns.get("proteobench") or {}
    if "proteobench" in target.varm:
        raise ValueError("varm['proteobench'] already exists; refusing to overwrite scores")
    if "scores" in namespace:
        raise ValueError(
            "uns['proteobench']['scores'] already exists; refusing to overwrite scores"
        )

    rule, roles = resolve_roles(target, module_settings, tool_settings)
    design = align_runs(target, rule, roles, module_settings, tool_settings)
    intermediate = compute_intermediate(
        target,
        module_settings,
        tool_settings,
        roles,
        design,
    )
    scores = build_scores(
        intermediate.legacy,
        intermediate.intermediate_hash,
        default_cutoff=module_settings.general.default_cutoff_min_feature,
        max_nr_observed=module_settings.general.max_nr_observed,
    )

    target.varm["proteobench"] = intermediate.varm
    namespace = dict(namespace)
    namespace.update(
        {
            "column_roles": roles.as_dict(),
            "protein_mapping": intermediate.protein_mapping,
            "scores": scores,
        }
    )
    target.uns["proteobench"] = namespace

    apb_namespace = dict(target.uns.get("anndata_proteomics") or {})
    apb_namespace["proteobench"] = {
        "schema_version": _STORAGE_SCHEMA_VERSION,
        "compatibility_version": PROTEOBENCH_COMPATIBILITY_VERSION,
        "source_revision": PROTEOBENCH_SOURCE_REVISION,
    }
    target.uns["anndata_proteomics"] = apb_namespace
    return obj
