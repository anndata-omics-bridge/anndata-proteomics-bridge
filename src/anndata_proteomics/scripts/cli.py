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
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Annotated, Any

from anndata import AnnData
from cyclopts import App, Parameter
from loguru import logger
from mudata import MuData
from pandas import DataFrame

from anndata_proteomics._logging import configure_default_sink
from anndata_proteomics._matrix_types import named_layers
from anndata_proteomics.converters.assemble import convert as _run_convert
from anndata_proteomics.readers.dispatch import read_table
from anndata_proteomics.rules import _export_schema
from anndata_proteomics.rules.loader import load_rule, load_rule_document
from anndata_proteomics.rules.registry import iter_packaged_rules
from anndata_proteomics.rules.schema import QuantificationLevel
from anndata_proteomics.rules.validate import (
    log_and_exit_code,
    validate_all_packaged,
    validate_file,
)

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
class _ConvertRequest:
    data: Path
    level: QuantificationLevel | None
    params: Path | None
    rule_config: Path | None
    software: str | None
    params_software: str | None
    output: Path | None


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
    *,
    params: Path | None = None,
    rule_config: Path | None = None,
    software: str | None = None,
    params_software: str | None = None,
    output: Path | None = None,
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
    """
    request = _ConvertRequest(
        data=data,
        level=level,
        params=params,
        rule_config=rule_config,
        software=software,
        params_software=params_software,
        output=output,
    )
    if output is not None and output.suffix:
        logger.error(
            f"--output must be an extensionless basename, got {output}; APB chooses .h5ad or .h5mu"
        )
        return 2
    df = read_table(data)
    if rule_config is not None:
        return _convert_with_rule_document(df, request)
    return _convert_with_packaged_rules(df, request)


def _convert_with_rule_document(df: DataFrame, request: _ConvertRequest) -> int:
    from anndata_proteomics.converters.pipeline import (
        attach_parameter_resolution,
        build_mudata_from_rules,
        matching_rules,
        resolve_parameters,
        set_rule_selection_method,
        software_slug,
    )

    assert request.rule_config is not None
    document = load_rule_document(request.rule_config)
    rule_slug = software_slug(document.software_name)
    parameter_resolution = (
        resolve_parameters(request.params, request.params_software or rule_slug)
        if request.params is not None
        else None
    )
    search_parameters = (
        parameter_resolution.parameters if parameter_resolution is not None else None
    )
    if request.level is not None:
        if request.level not in document.levels:
            logger.error(
                f"{request.rule_config} has no level {request.level!r}; "
                f"available: {list(document.levels)}"
            )
            return 1
        adata = _run_convert(
            df,
            document.effective_rule(request.level, search_parameters),
            params_path=None if parameter_resolution is not None else request.params,
        )
        if parameter_resolution is not None:
            attach_parameter_resolution(
                adata,
                parameter_resolution,
                selection_method="rule_config",
                warn_missing=True,
            )
        else:
            set_rule_selection_method(adata, "rule_config")
        return _write_anndata(adata, request.output, request.data)
    rules = matching_rules(document.effective_rules(search_parameters), df.columns)
    if not rules:
        logger.error(f"no level in {request.rule_config} matches the columns in {request.data}")
        return 1
    md = build_mudata_from_rules(
        df,
        rules,
        params_path=request.params,
        parameter_resolution=parameter_resolution,
        rule_selection_method="rule_config",
        software=rule_slug,
    )
    return _write_mudata(md, request.output, request.data)


def _convert_with_packaged_rules(df: DataFrame, request: _ConvertRequest) -> int:
    from anndata_proteomics.converters.pipeline import (
        build_mudata,
        convert_level,
        convertible_levels,
        recognize_software,
        resolve_parameters,
        resolve_rule_version,
    )

    slug = request.software or recognize_software(df.columns)
    if slug is None:
        logger.error(
            f"could not auto-detect the vendor for {request.data}; "
            "pass --software SLUG or --rule-config PATH"
        )
        return 1
    if request.params is None:
        logger.error("pass --params (it gives the software version) or --rule-config PATH")
        return 1
    parameter_resolution = resolve_parameters(
        request.params,
        request.params_software or slug,
    )
    version, version_status = resolve_rule_version(parameter_resolution, slug)
    logger.info(f"vendor={slug} software_version={version!r} version_status={version_status}")
    if request.level is not None:
        adata = convert_level(
            df,
            slug,
            request.level,
            version,
            params_path=request.params,
            parameter_resolution=parameter_resolution,
        )
        return _write_anndata(adata, request.output, request.data)
    levels = convertible_levels(
        slug,
        version,
        df.columns,
        version_status=version_status,
        search_parameters=parameter_resolution.parameters,
    )
    if not levels:
        logger.error(
            f"no quantification level resolves for {slug} at software version {version!r}; "
            "check --params / --software"
        )
        return 1
    md = build_mudata(
        df,
        slug,
        version,
        params_path=request.params,
        parameter_resolution=parameter_resolution,
    )
    return _write_mudata(md, request.output, request.data)


def _write_mudata(md: MuData, output: Path | None, data: Path) -> int:
    out = _output_path(output, data, ".h5mu")
    _write_atomically(out, md.write_h5mu)
    _remove_stale_sibling(out, ".h5ad")
    logger.info(f"wrote {out}  obs={md.n_obs}  modalities={list(md.mod)}")
    return 0


def _write_anndata(adata: AnnData, output: Path | None, data: Path) -> int:
    """Write a single-level AnnData to .h5ad and log a one-line summary."""
    out = _output_path(output, data, ".h5ad")
    _write_atomically(out, adata.write_h5ad)
    _remove_stale_sibling(out, ".h5mu")
    logger.info(f"wrote {out}  shape={adata.shape}  layers={list(named_layers(adata))}")
    return 0


def _output_path(output: Path | None, data: Path, suffix: str) -> Path:
    """Resolve a result path, appending APB's chosen suffix to explicit basenames."""
    return data.with_suffix(suffix) if output is None else Path(f"{output}{suffix}")


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
    from anndata_proteomics.readers.summary import describe_path

    result = describe_path(path, modality=modality)
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
    from anndata_proteomics.annotation.apply import annotate_obs
    from anndata_proteomics.annotation.loader import load_annotation
    from anndata_proteomics.readers.result import load_converted_result

    obj = load_converted_result(data)
    annotation = load_annotation(annotations)
    annotate_obs(obj, annotation)

    out = output or data.with_name(f"{data.stem}.annotated{data.suffix}")
    if hasattr(obj, "mod"):
        obj.write_h5mu(out)
    else:
        obj.write_h5ad(out)
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
    from mudata import MuData

    from anndata_proteomics.annotation.validate_fasta import (
        FastaValidationConfig,
        validate_peptide_modalities_against_fasta,
        validate_peptides_against_fasta,
    )
    from anndata_proteomics.annotation.var_fasta import (
        ProteinFastaAnnotationConfig,
        annotate_var_from_fasta,
    )
    from anndata_proteomics.fasta.config import FastaConfig
    from anndata_proteomics.readers.result import load_converted_result

    if not fasta_files:
        logger.error("no FASTA file given; usage: apb fasta DATA FASTA [FASTA ...]")
        return 1

    sources = list(fasta_files)
    obj = load_converted_result(data)
    has_protein = (
        "protein" in obj.mod if isinstance(obj, MuData) else _quantification_level(obj) == "protein"
    )
    peptide_modalities = (
        [
            name
            for name, target in obj.mod.items()
            if _quantification_level(target) in {"ion", "fragment", "peptidoform", "peptide"}
        ]
        if isinstance(obj, MuData)
        else (
            [_quantification_level(obj)]
            if _quantification_level(obj) in {"ion", "fragment", "peptidoform", "peptide"}
            else []
        )
    )
    identifier_config = FastaConfig.from_single_patterns(
        options.decoy_pattern,
        options.contaminant_pattern,
    )

    if has_protein:
        annotate_var_from_fasta(
            obj,
            sources,
            ProteinFastaAnnotationConfig(
                match_on=options.match_on,
                is_uniprot=options.is_uniprot,
                identifiers=identifier_config,
                cleavage=options.cleavage,
                min_length=options.min_length,
                max_length=options.max_length,
            ),
        )

    if options.validate and peptide_modalities:
        validation_config = FastaValidationConfig(
            sequence_field=options.sequence_field,
            backend=options.backend,
            identifiers=identifier_config,
            leading_protein_field=options.leading_protein_field,
            protein_match_on=options.match_on,
            il_equivalent=options.il_equivalent,
            is_uniprot=options.is_uniprot,
        )
        if isinstance(obj, MuData):
            results = validate_peptide_modalities_against_fasta(
                obj,
                sources,
                validation_config,
            )
        else:
            result = validate_peptides_against_fasta(
                obj,
                sources,
                validation_config,
            )
            results = {_quantification_level(obj) or "features": result}
        for name, result in results.items():
            logger.info(
                "{}: {}/{} peptide-derived features occur in FASTA",
                name,
                result.n_matched_features,
                result.n_features,
            )

    if not has_protein and (not options.validate or not peptide_modalities):
        logger.error(
            "input has no protein layer to annotate and no enabled peptide-derived "
            "layer to validate"
        )
        return 1

    out = options.output or data.with_name(f"{data.stem}.fasta{data.suffix}")
    if isinstance(obj, MuData):
        obj.write_h5mu(out)
    else:
        obj.write_h5ad(out)
    logger.info(f"wrote {out}")
    return 0


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
    from anndata_proteomics.proteobench.config import load_module_settings
    from anndata_proteomics.proteobench.pipeline import score_quantification
    from anndata_proteomics.readers.result import load_converted_result

    obj = load_converted_result(data)
    score_quantification(
        obj,
        load_module_settings(module_settings),
    )

    out = output or data.with_name(f"{data.stem}.proteobench{data.suffix}")
    if out.suffix != data.suffix:
        raise ValueError(f"ProteoBench output must keep the input container suffix {data.suffix!r}")
    if out.resolve() == data.resolve():
        raise ValueError("ProteoBench output must differ from the input path")
    writer = obj.write_h5mu if hasattr(obj, "mod") else obj.write_h5ad
    _write_atomically(out, writer)
    logger.info(f"wrote {out}")
    return 0


def _quantification_level(obj: Any) -> str | None:
    """Return APB's stored quantification level for one AnnData-like object."""
    return (obj.uns.get("anndata_proteomics") or {}).get("quantification_level")


def main() -> int:
    """Console-script entry point."""
    configure_default_sink()
    rc = app()
    return int(rc) if rc is not None else 0


if __name__ == "__main__":
    sys.exit(main())
