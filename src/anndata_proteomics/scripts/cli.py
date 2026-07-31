"""apb CLI dispatcher.

Subcommands:
- validate [path ...]        validate one or more JSON rules; defaults to all packaged
- list                       list packaged rules
- export-schema              regenerate parse_rule.schema.json
- convert <data> [level]     convert a vendor file to MuData (.h5mu) or one level to AnnData (.h5ad)
- summary <path>             print a lightweight container metadata view
- annotate <data> <json>     join sample annotations onto obs
- fasta <data> <fasta...>    annotate proteins and validate peptide identifications
- proteobench <data> <module.toml>  compute quantitative benchmark scores
"""

from __future__ import annotations

import json as jsonlib
import sys
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Annotated

from anndata import AnnData
from cyclopts import App, Parameter
from loguru import logger
from mudata import MuData
from pandas import DataFrame

from anndata_proteomics._logging import configure_default_sink
from anndata_proteomics.adapters.anndata import annotation as annotation_adapter
from anndata_proteomics.adapters.anndata import conversion as conversion_adapter
from anndata_proteomics.adapters.anndata import description as description_adapter
from anndata_proteomics.adapters.anndata import fasta as fasta_adapter
from anndata_proteomics.adapters.anndata import proteobench as proteobench_adapter
from anndata_proteomics.adapters.anndata import result as result_adapter
from anndata_proteomics.adapters.anndata import rules as rules_adapter
from anndata_proteomics.adapters.anndata.matrix import layer_names
from anndata_proteomics.annotation import loader as annotation_loader
from anndata_proteomics.annotation.validate_fasta import PeptideFastaMatchingConfig
from anndata_proteomics.annotation.var_fasta import annotate_proteins_from_fasta
from anndata_proteomics.converters import pipeline as conversion_pipeline
from anndata_proteomics.fasta.config import (
    DEFAULT_CONTAMINANT_POLICY,
    DEFAULT_DECOY_POLICY,
    ExplicitPatterns,
    FastaConfig,
)
from anndata_proteomics.proteobench import config as proteobench_config
from anndata_proteomics.readers.dispatch import (
    read_table_columns,
    read_table_preserving_strings,
)
from anndata_proteomics.rules import _export_schema
from anndata_proteomics.rules.loader import load_rule, load_rule_document
from anndata_proteomics.rules.registry import iter_packaged_rules
from anndata_proteomics.rules.schema import PEPTIDE_LEVELS, QuantificationLevel
from anndata_proteomics.rules.validate import (
    log_and_exit_code,
    validate_all_packaged,
    validate_file,
)
from anndata_proteomics.workflows import conversion as conversion_workflow
from anndata_proteomics.workflows import fasta as fasta_workflow
from anndata_proteomics.workflows import proteobench as proteobench_workflow
from anndata_proteomics.workflows import sample_annotation as annotation_workflow

app = App(name="apb", help="anndata_proteomics (APB) CLI", help_on_error=True)


@dataclass(frozen=True, slots=True)
class FastaCliOptions:
    """Flat Cyclopts option group for ``apb fasta``."""

    output: Path | None = None
    match_on: str | None = None
    is_uniprot: bool = True
    decoy_pattern: str | None = None
    contaminant_pattern: str | None = None
    cleavage: str | None = None
    min_length: int | None = None
    max_length: int | None = None
    validate: bool = True
    sequence_field: str = "ProForma_peptide"
    leading_protein_field: str | None = None
    backend: str = "auto"
    il_equivalent: bool = False


DEFAULT_FASTA_CLI_OPTIONS = FastaCliOptions()


@dataclass(frozen=True, slots=True)
class ConvertCliOptions:
    """Flat Cyclopts option group for ``apb convert``."""

    params: Path | None = None
    rule_config: Path | None = None
    software: str | None = None
    params_software: str | None = None
    output: Path | None = None
    strict: bool = False


DEFAULT_CONVERT_CLI_OPTIONS = ConvertCliOptions()


@dataclass(frozen=True, slots=True)
class ConversionContext:
    """Resolved file inputs shared by one CLI conversion route."""

    data_path: Path
    headers: tuple[str, ...]
    strict: bool
    output: Path


@dataclass(frozen=True, slots=True)
class PackagedConversionContext:
    """Resolved packaged-rule inputs shared by one conversion route."""

    slug: str
    parameters: conversion_pipeline.ParameterResolution


@dataclass(frozen=True, slots=True)
class PackagedConversionUnavailable:
    """Packaged-rule inputs could not be resolved; the diagnostic was logged."""


type PackagedConversionResolution = PackagedConversionContext | PackagedConversionUnavailable


@app.command
def validate(*paths: Path) -> int:
    """Validate one or more JSON rule files.

    With no paths, walks all packaged rules (same as `validate-rules`).
    """
    if not paths:
        results = validate_all_packaged()
    else:
        results = [validate_file(p) for p in paths]
    return log_and_exit_code(results)


@app.command(name="list")
def list_rules() -> int:
    """List packaged parsing rules: software, level, file_version, version pattern, path."""
    for locator in iter_packaged_rules():
        rule = load_rule(locator)
        logger.info(
            f"{rule.software_name:14}  {rule.quantification_level:12}  "
            f"v{rule.file_version:<3}  {rule.software_version:14}  "
            f"{locator.path}#{locator.level}"
        )
    return 0


@app.command(name="export-schema")
def export_schema_cmd() -> int:
    """Regenerate parse_rule.schema.json from the pydantic models."""
    _export_schema.main()
    return 0


@app.command
def convert(
    data: Path,
    level: QuantificationLevel | None = None,
    options: Annotated[ConvertCliOptions, Parameter(name="*")] = DEFAULT_CONVERT_CLI_OPTIONS,
) -> int:
    """Convert a vendor file to a multi-level MuData (.h5mu) or one level to an AnnData (.h5ad).

    With no LEVEL, every quantification level the file/version provides is wrapped into a MuData
    (.h5mu) on a shared run axis, including a one-modality MuData for single-level vendors.
    Pass a LEVEL (ion / fragment / peptidoform / peptide / protein) to emit an AnnData (.h5ad).

    --params is the vendor parameter file and is required unless --rule-config is given. The
    result-table vendor is auto-detected from the column headers; override it with --software (the
    rule folder slug, e.g. "diann"). For compound workflows such as FragPipe with DIA-NN output,
    --params-software selects the parameter parser independently (e.g. "fragpipe"). --rule-config
    selects an explicit software-version document; LEVEL chooses one section, while omitting LEVEL
    converts every matching section in that document.
    --output is an extensionless basename. APB appends .h5mu for MuData or .h5ad for AnnData.
    Without --output, the result is written next to the input using the input stem.
    --strict promotes layer-contract warnings to errors. An empty X layer is always an
    error; --strict extends that to every other declared layer.
    """
    if options.output is not None and options.output.suffix:
        logger.error(
            f"--output must be an extensionless basename, got {options.output}; "
            "APB chooses .h5ad or .h5mu"
        )
        return 2
    headers = read_table_columns(data)
    suffix = ".h5ad" if level is not None else ".h5mu"
    output = (
        data.with_suffix(suffix) if options.output is None else Path(f"{options.output}{suffix}")
    )
    context = ConversionContext(data, tuple(headers), options.strict, output)
    if options.rule_config is not None:
        if level is not None:
            return _convert_level_from_rule_config(
                context,
                options.rule_config,
                level,
                options,
            )
        return _convert_levels_from_rule_config(context, options.rule_config, options)
    if level is not None:
        return _convert_level_from_packaged_rules(context, level, options)
    return _convert_levels_from_packaged_rules(context, options)


def _convert_level_from_rule_config(
    context: ConversionContext,
    rule_config: Path,
    level: QuantificationLevel,
    options: ConvertCliOptions,
) -> int:
    """Select and execute one level from an explicit rule document."""
    document = load_rule_document(rule_config)
    if level not in document.levels:
        logger.error(f"{rule_config} has no level {level!r}; available: {list(document.levels)}")
        return 1
    if options.params is not None:
        resolution = conversion_pipeline.resolve_parameters(
            options.params,
            options.params_software or conversion_pipeline.software_slug(document.software_name),
        )
        selection = conversion_pipeline.RuleSelection(
            document.parameterized_effective_rule(level, resolution.parameters),
            "rule_config",
        )
        return _execute_parameterized_level(context, level, selection, resolution)
    selection = conversion_pipeline.RuleSelection(document.effective_rule(level), "rule_config")
    return _execute_level(context, level, selection)


def _convert_levels_from_rule_config(
    context: ConversionContext,
    rule_config: Path,
    options: ConvertCliOptions,
) -> int:
    """Select and execute all matching levels from an explicit rule document."""
    document = load_rule_document(rule_config)
    if options.params is None:
        return _execute_matching_levels(context, document.effective_rules(), rule_config)
    resolution = conversion_pipeline.resolve_parameters(
        options.params,
        options.params_software or conversion_pipeline.software_slug(document.software_name),
    )
    rules = document.parameterized_effective_rules(resolution.parameters)
    return _execute_matching_parameterized_levels(
        context,
        rules,
        rule_config,
        resolution,
    )


def _convert_level_from_packaged_rules(
    context: ConversionContext,
    level: QuantificationLevel,
    options: ConvertCliOptions,
) -> int:
    """Select and execute one packaged quantification level."""
    resolved = _resolve_packaged_conversion(context, options)
    if isinstance(resolved, PackagedConversionUnavailable):
        return 1
    selection = conversion_workflow.select_rule_from_parameters(
        context.headers,
        resolved.slug,
        level,
        resolved.parameters,
    )
    return _execute_parameterized_level(
        context,
        level,
        selection,
        resolved.parameters,
    )


def _convert_levels_from_packaged_rules(
    context: ConversionContext,
    options: ConvertCliOptions,
) -> int:
    """Select and execute every available packaged quantification level."""
    resolved = _resolve_packaged_conversion(context, options)
    if isinstance(resolved, PackagedConversionUnavailable):
        return 1
    selections = conversion_workflow.select_rules_from_parameters(
        context.headers,
        resolved.slug,
        resolved.parameters,
    )
    if not selections:
        return _log_no_packaged_levels(resolved.slug, resolved.parameters)
    return _execute_parameterized_levels(context, selections, resolved.parameters)


def _resolve_packaged_conversion(
    context: ConversionContext,
    options: ConvertCliOptions,
) -> PackagedConversionResolution:
    """Resolve vendor identity and parameters for a packaged-rule route."""
    if options.software is not None:
        slug = options.software
    else:
        recognition = conversion_pipeline.recognize_software(context.headers)
        if not isinstance(recognition, conversion_pipeline.RecognizedSoftware):
            logger.error(
                f"could not auto-detect the vendor for {context.data_path}; "
                "pass --software SLUG or --rule-config PATH"
            )
            return PackagedConversionUnavailable()
        slug = recognition.slug
    if options.params is None:
        logger.error("pass --params (it gives the software version) or --rule-config PATH")
        return PackagedConversionUnavailable()
    resolution = conversion_pipeline.resolve_parameters(
        options.params,
        options.params_software or slug,
    )
    _log_parameter_resolution(slug, resolution)
    return PackagedConversionContext(slug, resolution)


def _log_no_packaged_levels(
    slug: str,
    resolution: conversion_pipeline.ParameterResolution,
) -> int:
    """Report that parameter evidence selected no packaged level."""
    version = conversion_pipeline.resolve_rule_version(resolution, slug)
    version_label = (
        version.value if isinstance(version, conversion_pipeline.PresentRuleVersion) else "missing"
    )
    logger.error(
        f"no quantification level resolves for {slug} at software version "
        f"{version_label!r}; check --params / --software"
    )
    return 1


def _execute_level(
    context: ConversionContext,
    level: QuantificationLevel,
    selection: conversion_pipeline.RuleSelection,
) -> int:
    """Read, convert, and write one already-selected level."""
    data = _read_table_for_selections(context.data_path, (selection,))
    conversion = conversion_workflow.convert_selected_level(
        data,
        level,
        selection,
        strict=context.strict,
    )
    return _write_anndata(conversion_adapter.to_anndata(conversion), context.output)


def _execute_parameterized_level(
    context: ConversionContext,
    level: QuantificationLevel,
    selection: conversion_pipeline.RuleSelection,
    resolution: conversion_pipeline.ParameterResolution,
) -> int:
    """Read, convert, and write one selected level with parameter provenance."""
    data = _read_table_for_selections(context.data_path, (selection,))
    conversion = conversion_workflow.convert_selected_level(
        data,
        level,
        selection,
        strict=context.strict,
    )
    adata = conversion_adapter.to_anndata(conversion)
    _write_parameter_resolution(adata, resolution)
    return _write_anndata(adata, context.output)


def _execute_levels(
    context: ConversionContext,
    selections: dict[QuantificationLevel, conversion_pipeline.RuleSelection],
) -> int:
    """Read, convert, and write several already-selected levels."""
    data = _read_table_for_selections(context.data_path, selections.values())
    conversions = conversion_workflow.convert_selected_levels(
        data,
        selections,
        strict=context.strict,
    )
    return _write_mudata(conversion_adapter.to_mudata(conversions), context.output)


def _execute_parameterized_levels(
    context: ConversionContext,
    selections: dict[QuantificationLevel, conversion_pipeline.RuleSelection],
    resolution: conversion_pipeline.ParameterResolution,
) -> int:
    """Read, convert, and write selected levels with parameter provenance."""
    data = _read_table_for_selections(context.data_path, selections.values())
    conversions = conversion_workflow.convert_selected_levels(
        data,
        selections,
        strict=context.strict,
    )
    md = conversion_adapter.to_mudata(conversions)
    _write_parameter_resolution(md, resolution)
    return _write_mudata(md, context.output)


def _execute_matching_levels(
    context: ConversionContext,
    rules: dict[QuantificationLevel, conversion_pipeline.ParseRule],
    rule_config: Path,
) -> int:
    """Select configured rules matching the headers, then execute them."""
    matching = conversion_pipeline.matching_rules(rules, context.headers)
    if not matching:
        logger.error(f"no level in {rule_config} matches the input columns")
        return 1
    return _execute_levels(context, _configured_selections(matching))


def _execute_matching_parameterized_levels(
    context: ConversionContext,
    rules: dict[QuantificationLevel, conversion_pipeline.ParseRule],
    rule_config: Path,
    resolution: conversion_pipeline.ParameterResolution,
) -> int:
    """Select configured parameterized rules matching headers, then execute them."""
    matching = conversion_pipeline.matching_rules(rules, context.headers)
    if not matching:
        logger.error(f"no level in {rule_config} matches the input columns")
        return 1
    return _execute_parameterized_levels(
        context,
        _configured_selections(matching),
        resolution,
    )


def _configured_selections(
    rules: dict[QuantificationLevel, conversion_pipeline.ParseRule],
) -> dict[QuantificationLevel, conversion_pipeline.RuleSelection]:
    """Attach explicit-config provenance to concrete effective rules."""
    return {
        level: conversion_pipeline.RuleSelection(rule, "rule_config")
        for level, rule in rules.items()
    }


def _read_table_for_selections(
    data_path: Path,
    selections: Iterable[conversion_pipeline.RuleSelection],
) -> DataFrame:
    """Read one vendor table using the logical contracts of selected rules."""
    rules = tuple(selection.rule for selection in selections)
    string_sources = conversion_pipeline.string_sources_for_rules(rules)
    return read_table_preserving_strings(data_path, string_sources)


def _log_parameter_resolution(
    slug: str,
    resolution: conversion_pipeline.ParameterResolution,
) -> None:
    """Log the concrete software-version result used for packaged rule selection."""
    version = conversion_pipeline.resolve_rule_version(resolution, slug)
    version_label = (
        version.value if isinstance(version, conversion_pipeline.PresentRuleVersion) else "missing"
    )
    logger.info(
        "vendor={} software_version={} version_status={}",
        slug,
        version_label,
        conversion_pipeline.version_status(version),
    )


def _write_parameter_resolution(
    target: AnnData | MuData,
    resolution: conversion_pipeline.ParameterResolution,
) -> None:
    """Persist one parsed parameter result across one converted artifact."""
    conversion_adapter.write_parameter_resolution(target, resolution)
    if isinstance(target, MuData):
        for modality in target.mod.values():
            conversion_adapter.write_parameter_resolution(modality, resolution)
    if isinstance(resolution.version, conversion_pipeline.MissingRuleVersion):
        logger.warning(
            "no software version in search parameters {}; selected rule by columns",
            resolution.source_path,
        )


def _write_mudata(md: MuData, output: Path) -> int:
    _write_atomically(output, md.write_h5mu)
    _remove_stale_sibling(output, ".h5ad")
    logger.info(f"wrote {output}  obs={md.n_obs}  modalities={list(md.mod)}")
    return 0


def _write_container(obj: AnnData | MuData, out: Path) -> None:
    """Write a loaded container back to ``out`` in its own on-disk format."""
    if isinstance(obj, MuData):
        obj.write_h5mu(out)
    else:
        obj.write_h5ad(out)


def _write_anndata(adata: AnnData, output: Path) -> int:
    """Write a single-level AnnData to .h5ad and log a one-line summary."""
    _write_atomically(output, adata.write_h5ad)
    _remove_stale_sibling(output, ".h5mu")
    logger.info(f"wrote {output}  shape={adata.shape}  layers={list(layer_names(adata))}")
    return 0


def _write_atomically(output: Path, writer: Callable[[Path], None]) -> None:
    """Write beside the destination and replace it only after a complete write."""
    with TemporaryDirectory(dir=output.parent, prefix=f".{output.name}.") as folder:
        temporary = Path(folder) / output.name
        writer(temporary)
        temporary.replace(output)


def _remove_stale_sibling(output: Path, stale_suffix: str) -> None:
    """Remove the prior alternate container type after a successful replacement."""
    stale = output.with_suffix(stale_suffix)
    if stale.exists():
        stale.unlink()
        logger.info(f"removed stale {stale}")


@app.command(name="summary")
def summary_cmd(
    path: Path,
    *,
    modality: str | None = None,
    json: bool = False,
) -> int:
    """Print APB's lightweight shape and stored-metadata view."""
    result = (
        description_adapter.describe_path(path)
        if modality is None
        else description_adapter.describe_modality_path(path, modality)
    )
    print(jsonlib.dumps(result, indent=None if json else 2, sort_keys=True))
    return 0


@app.command
def annotate(
    data: Path,
    annotations: Path,
    output: Path | None = None,
) -> int:
    """Join a sample-annotation table onto obs and write the result.

    Reads ProteoBench ``module_settings.toml`` directly, as well as APB annotation TOML,
    CSV, and TSV tables. Joins sample records onto obs by run/file name and writes the
    enriched object. --output defaults to ``<stem>.annotated<suffix>`` next to the input.
    """
    obj = result_adapter.load_converted_result(data)
    loaded = annotation_loader.load_annotation(annotations)
    observation_frames = annotation_adapter.read_observation_frames(obj)
    result = annotation_workflow.run_sample_annotation(
        observation_frames,
        loaded.table,
        loaded.origin,
    )
    annotation_adapter.write_sample_annotation(obj, result)

    out = output or data.with_name(f"{data.stem}.annotated{data.suffix}")
    _write_container(obj, out)
    logger.info(f"wrote {out}")
    return 0


@app.command
def fasta(
    data: Path,
    *fasta_files: Path,
    options: Annotated[FastaCliOptions, Parameter(name="*")] = DEFAULT_FASTA_CLI_OPTIONS,
) -> int:
    """Annotate proteins and validate peptide identifications against FASTA.

    Protein-derived layers receive ``varm['fasta']``. By default, every
    peptide-derived layer is also checked with one Aho--Corasick scan and receives
    ``varm['fasta_validation']``; use ``--no-validate`` to disable that check.
    MuData peptide-to-protein relationships are added to
    ``varp['feature_mapping']`` in MuLink format. Decoy and contaminant patterns
    are inferred from the supplied FASTA unless explicitly provided. They classify
    records but never filter quantified features. Unmatched peptides are retained.
    Protein joins use the stored ``column_roles.fasta_accessions`` declaration;
    ``--match-on`` and ``--leading-protein-field`` are explicit overrides.
    """
    if not fasta_files:
        logger.error("no FASTA file given; usage: apb fasta DATA FASTA [FASTA ...]")
        return 1

    sources = tuple(fasta_files)
    obj = result_adapter.load_converted_result(data)
    if isinstance(obj, MuData):
        has_protein = "protein" in obj.mod
        peptide_modalities = [
            name
            for name, target in obj.mod.items()
            if isinstance(target, AnnData)
            and rules_adapter.require_quantification_level(target) in PEPTIDE_LEVELS
        ]
    else:
        level = rules_adapter.require_quantification_level(obj)
        has_protein = level == "protein"
        peptide_modalities = [level] if level in PEPTIDE_LEVELS else []
    decoy_policy = (
        DEFAULT_DECOY_POLICY
        if options.decoy_pattern is None
        else ExplicitPatterns(patterns=(options.decoy_pattern,) if options.decoy_pattern else ())
    )
    contaminant_policy = (
        DEFAULT_CONTAMINANT_POLICY
        if options.contaminant_pattern is None
        else ExplicitPatterns(
            patterns=(options.contaminant_pattern,) if options.contaminant_pattern else ()
        )
    )
    identifier_config = FastaConfig(
        decoy=decoy_policy,
        contaminant=contaminant_policy,
    )
    match_on = (
        fasta_adapter.STORED_FASTA_ACCESSIONS
        if options.match_on is None
        else (
            fasta_adapter.FASTA_ACCESSIONS_INDEX
            if options.match_on == "index"
            else fasta_adapter.FastaAccessionsColumn(options.match_on)
        )
    )
    cleavage = (
        fasta_workflow.STORED_CLEAVAGE_OR_TRYPSIN
        if options.cleavage is None
        else fasta_workflow.NamedCleavage(options.cleavage)
    )
    minimum_length = (
        fasta_workflow.STORED_MINIMUM_LENGTH_OR_DEFAULT
        if options.min_length is None
        else fasta_workflow.MinimumPeptideLength(options.min_length)
    )
    maximum_length = (
        fasta_workflow.STORED_MAXIMUM_LENGTH_OR_DEFAULT
        if options.max_length is None
        else fasta_workflow.MaximumPeptideLength(options.max_length)
    )

    if has_protein:
        protein_config = fasta_adapter.AnnDataProteinFastaConfig(
            match_on=match_on,
            is_uniprot=options.is_uniprot,
            identifiers=identifier_config,
            cleavage=cleavage,
            minimum_length=minimum_length,
            maximum_length=maximum_length,
        )
        if isinstance(obj, MuData):
            _annotate_protein_mudata(obj, sources, protein_config)
        else:
            _annotate_protein_anndata(obj, sources, protein_config)

    if options.validate and peptide_modalities:
        reported_proteins = (
            fasta_adapter.STORED_REPORTED_PROTEINS
            if options.leading_protein_field is None
            else fasta_adapter.ReportedProteinsColumn(options.leading_protein_field)
        )
        validation_config = fasta_adapter.AnnDataPeptideFastaConfig(
            sequence_field=options.sequence_field,
            matching=PeptideFastaMatchingConfig(
                backend=options.backend,
                identifiers=identifier_config,
                il_equivalent=options.il_equivalent,
                is_uniprot=options.is_uniprot,
            ),
            reported_proteins=reported_proteins,
        )
        if isinstance(obj, MuData):
            _validate_peptide_mudata(
                obj,
                sources,
                fasta_adapter.MuDataPeptideFastaConfig(
                    validation=validation_config,
                    protein_match_on=match_on,
                ),
            )
        else:
            _validate_peptide_anndata(
                obj,
                sources,
                validation_config,
            )

    if not has_protein and (not options.validate or not peptide_modalities):
        logger.error(
            "input has no protein layer to annotate and no enabled peptide-derived "
            "layer to validate"
        )
        return 1

    out = options.output or data.with_name(f"{data.stem}.fasta{data.suffix}")
    _write_container(obj, out)
    logger.info(f"wrote {out}")
    return 0


def _annotate_protein_anndata(
    target: AnnData,
    sources: tuple[Path, ...],
    config: fasta_adapter.AnnDataProteinFastaConfig,
) -> None:
    """Compose protein extraction, calculation, and AnnData persistence."""
    extracted = fasta_workflow.resolve_protein_annotation_input(
        fasta_adapter.read_protein_annotation_input(target, config)
    )
    result = annotate_proteins_from_fasta(
        extracted.protein_groups,
        sources,
        extracted.config,
    )
    fasta_adapter.write_protein_annotation(target, result, extracted.provenance)
    fasta_adapter.write_anndata_fasta_config(target, result.fasta_config)


def _annotate_protein_mudata(
    target: MuData,
    sources: tuple[Path, ...],
    config: fasta_adapter.AnnDataProteinFastaConfig,
) -> None:
    """Compose protein extraction, calculation, and MuData persistence."""
    protein = fasta_adapter.require_protein_mudata_target(target)
    extracted = fasta_workflow.resolve_protein_annotation_input(
        fasta_adapter.read_protein_annotation_input(protein, config)
    )
    result = annotate_proteins_from_fasta(
        extracted.protein_groups,
        sources,
        extracted.config,
    )
    fasta_adapter.write_protein_annotation(protein, result, extracted.provenance)
    fasta_adapter.write_mudata_fasta_config(target, result.fasta_config)


def _validate_peptide_anndata(
    target: AnnData,
    sources: tuple[Path, ...],
    config: fasta_adapter.AnnDataPeptideFastaConfig,
) -> None:
    """Compose peptide extraction, shared-scan calculation, and AnnData persistence."""
    extracted = fasta_adapter.read_peptide_anndata_input(target, config)
    workflow = fasta_workflow.validate_peptide_levels((extracted,), sources, config.matching)
    result = workflow.levels[extracted.name]
    fasta_adapter.write_peptide_validation(target, result, config.sequence_field)
    fasta_adapter.write_anndata_fasta_config(target, workflow.fasta_config)
    _log_peptide_validation(extracted.name, result)


def _validate_peptide_mudata(
    target: MuData,
    sources: tuple[Path, ...],
    config: fasta_adapter.MuDataPeptideFastaConfig,
) -> None:
    """Compose shared-scan peptide validation and optional MuLink persistence."""
    extracted = fasta_adapter.read_peptide_mudata_inputs(target, config)
    workflow = fasta_workflow.validate_peptide_levels(
        extracted,
        sources,
        config.validation.matching,
    )
    if fasta_adapter.has_protein_mudata_target(target):
        mapping_input = fasta_adapter.read_feature_mapping_input(
            target,
            tuple(workflow.levels),
            config.protein_match_on,
        )
        mapping = fasta_workflow.build_mulink_feature_mapping(mapping_input.axes, workflow)
        state = fasta_adapter.read_feature_mapping_state(
            target,
            mapping.feature_mapping.mapping.shape,
        )
        update = fasta_workflow.calculate_feature_mapping_update(
            state,
            mapping,
            mapping_input.protein_match_on,
        )
        fasta_adapter.write_feature_mapping(target, update)
        for name, result in workflow.levels.items():
            modality = fasta_adapter.require_peptide_mudata_target(target, name)
            fasta_adapter.write_peptide_validation_with_feature_mapping(
                modality,
                result,
                config.validation.sequence_field,
                update.stats,
            )
            _log_peptide_validation(name, result)
    else:
        for name, result in workflow.levels.items():
            modality = fasta_adapter.require_peptide_mudata_target(target, name)
            fasta_adapter.write_peptide_validation(
                modality,
                result,
                config.validation.sequence_field,
            )
            _log_peptide_validation(name, result)
    fasta_adapter.write_mudata_fasta_config(target, workflow.fasta_config)


def _log_peptide_validation(
    name: str,
    result: fasta_workflow.PeptideValidation,
) -> None:
    """Log the per-level match count returned by the FASTA workflow."""
    logger.info(
        "{}: {}/{} peptide-derived features occur in FASTA",
        name,
        result.matching.n_matched_features,
        result.matching.n_features,
    )


@app.command
def proteobench(
    data: Path,
    module_settings: Path,
    *,
    output: Path | None = None,
) -> int:
    """Compute ProteoBench scores from an APB-annotated quantitative object.

    Every quantification level the object holds is scored into its own ``uns``/``varm``: an
    AnnData scores itself, a MuData scores each modality. Run ``apb annotate`` before this
    command; scoring requires ``sample_name`` and ``condition`` in each scored observation
    table. The default output is ``<stem>.proteobench<suffix>`` beside the input.
    """
    obj = result_adapter.load_converted_result(data)
    targets = tuple(proteobench_adapter.resolve_targets(obj))
    extracted = tuple(proteobench_adapter.read_level(target) for target in targets)
    results = proteobench_workflow.score_levels(
        tuple(level.calculation for level in extracted),
        proteobench_config.load_module_settings(module_settings),
    )
    for target, level, result in zip(targets, extracted, results, strict=True):
        proteobench_adapter.store_result(target, result, level.roles)

    out = output or data.with_name(f"{data.stem}.proteobench{data.suffix}")
    if out.suffix != data.suffix:
        raise ValueError(f"ProteoBench output must keep the input container suffix {data.suffix!r}")
    if out.resolve() == data.resolve():
        raise ValueError("ProteoBench output must differ from the input path")
    _write_atomically(out, lambda path: _write_container(obj, path))
    logger.info(f"wrote {out}")
    return 0


def main() -> int:
    """Console-script entry point."""
    configure_default_sink()
    rc = app()
    return int(rc) if rc is not None else 0


if __name__ == "__main__":
    sys.exit(main())
