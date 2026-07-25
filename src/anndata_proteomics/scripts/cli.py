"""apb CLI dispatcher.

Subcommands:
- validate [path ...]        validate one or more JSON rules; defaults to all packaged
- list                       list packaged rules
- export-schema              regenerate parse_rule.schema.json
- convert <data> [level]     convert a vendor file to MuData (.h5mu) or one level to AnnData (.h5ad)
- summary <path>             print a stored descriptive summary
- annotate <data> <json>     join sample annotations onto obs
- fasta <data> <fasta...>    annotate proteins and validate peptide identifications
- proteobench <data> <module.toml>  compute quantitative benchmark scores
"""

from __future__ import annotations

import json as jsonlib
import sys
from collections.abc import Callable
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

from anndata import AnnData
from cyclopts import App
from loguru import logger

from anndata_proteomics._logging import configure_default_sink
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
def convert(  # noqa: C901, PLR0911 - CLI maps validation failures to exit statuses
    data: Path,
    level: QuantificationLevel | None = None,
    *,
    params: Path | None = None,
    rule_config: Path | None = None,
    software: str | None = None,
    output: Path | None = None,
) -> int:
    """Convert a vendor file to a multi-level MuData (.h5mu) or one level to an AnnData (.h5ad).

    With no LEVEL, every quantification level the file/version provides is wrapped into a MuData
    (.h5mu) on a shared run axis, including a one-modality MuData for single-level vendors.
    Pass a LEVEL (ion / fragment / peptidoform / peptide / protein) to emit an AnnData (.h5ad).

    --params is the vendor parameter file; it supplies the software version that selects the rule
    variant (e.g. DIA-NN v1 vs v2) and is required unless --rule-config is given. The vendor is
    auto-detected from the column headers; override with --software (the rule folder slug, e.g.
    "diann"). --rule-config selects an explicit software-version document; LEVEL chooses one
    section, while omitting LEVEL converts every matching section in that document.
    --output is an extensionless basename. APB appends .h5mu for MuData or .h5ad for AnnData.
    Without --output, the result is written next to the input using the input stem.
    """
    from anndata_proteomics.converters.pipeline import (
        attach_parameter_resolution,
        build_mudata,
        build_mudata_from_rules,
        convert_level,
        convertible_levels,
        matching_rules,
        recognize_software,
        resolve_parameters,
        set_rule_selection_method,
        software_slug,
    )

    if output is not None and output.suffix:
        logger.error(
            f"--output must be an extensionless basename, got {output}; APB chooses .h5ad or .h5mu"
        )
        return 2

    df = read_table(data)

    # --rule-config: explicit software-version document, bypassing packaged selection.
    if rule_config is not None:
        document = load_rule_document(rule_config)
        rule_slug = software_slug(document.software_name)
        parameter_resolution = resolve_parameters(params, rule_slug) if params is not None else None
        search_parameters = (
            parameter_resolution.parameters if parameter_resolution is not None else None
        )
        if level is not None:
            if level not in document.levels:
                logger.error(
                    f"{rule_config} has no level {level!r}; available: {list(document.levels)}"
                )
                return 1
            adata = _run_convert(
                df,
                document.effective_rule(level, search_parameters),
                params_path=None if parameter_resolution is not None else params,
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
            return _write_anndata(adata, output, data)
        rules = matching_rules(
            document.effective_rules(search_parameters),
            df.columns,
        )
        if rules:
            md = build_mudata_from_rules(
                df,
                rules,
                params_path=params,
                parameter_resolution=parameter_resolution,
                rule_selection_method="rule_config",
                software=rule_slug,
            )
            out = _output_path(output, data, ".h5mu")
            _write_atomically(out, md.write_h5mu)
            _remove_stale_sibling(out, ".h5ad")
            logger.info(f"wrote {out}  obs={md.n_obs}  modalities={list(md.mod)}")
            return 0
        logger.error(f"no level in {rule_config} matches the columns in {data}")
        return 1

    slug = software or recognize_software(df.columns)
    if slug is None:
        logger.error(
            f"could not auto-detect the vendor for {data}; "
            "pass --software SLUG or --rule-config PATH"
        )
        return 1
    if params is None:
        logger.error("pass --params (it gives the software version) or --rule-config PATH")
        return 1
    parameter_resolution = resolve_parameters(params, slug)
    version = parameter_resolution.version
    logger.info(
        f"vendor={slug} software_version={version!r} "
        f"version_status={parameter_resolution.version_status}"
    )

    if level is not None:
        adata = convert_level(
            df,
            slug,
            level,
            version,
            params_path=params,
            parameter_resolution=parameter_resolution,
        )
        return _write_anndata(adata, output, data)

    levels = convertible_levels(
        slug,
        version,
        df.columns,
        version_status=parameter_resolution.version_status,
        search_parameters=parameter_resolution.parameters,
    )
    if levels:
        md = build_mudata(
            df,
            slug,
            version,
            params_path=params,
            parameter_resolution=parameter_resolution,
        )
        out = _output_path(output, data, ".h5mu")
        _write_atomically(out, md.write_h5mu)
        _remove_stale_sibling(out, ".h5ad")
        logger.info(f"wrote {out}  obs={md.n_obs}  modalities={list(md.mod)}")
        return 0
    logger.error(
        f"no quantification level resolves for {slug} at software version {version!r}; "
        "check --params / --software"
    )
    return 1


def _write_anndata(adata: AnnData, output: Path | None, data: Path) -> int:
    """Write a single-level AnnData to .h5ad and log a one-line summary."""
    out = _output_path(output, data, ".h5ad")
    _write_atomically(out, adata.write_h5ad)
    _remove_stale_sibling(out, ".h5mu")
    logger.info(f"wrote {out}  shape={adata.shape}  layers={list(adata.layers)}")
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
    """Print APB's stored descriptive summary for a container or MuData modality."""
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
def fasta(  # noqa: PLR0913 - stable CLI option surface
    data: Path,
    *fasta_files: Path,
    output: Path | None = None,
    match_on: str | None = None,
    is_uniprot: bool = True,
    decoy_pattern: str | None = None,
    contaminant_pattern: str | None = None,
    cleavage: str | None = None,
    min_length: int | None = None,
    max_length: int | None = None,
    validate: bool = True,
    sequence_field: str = "ProForma_peptide",
    leading_protein_field: str | None = None,
    backend: str = "auto",
    il_equivalent: bool = False,
) -> int:
    """Annotate proteins and validate peptide identifications against FASTA.

    Protein-derived layers receive ``varm['fasta']``. By default, every
    peptide-derived layer is also checked with one Aho--Corasick scan and receives
    ``varm['fasta_validation']``; use ``--no-validate`` to disable that check.
    MuData peptide-to-protein relationships are added to
    ``varp['feature_mapping']`` in MuLink format. Decoy and contaminant patterns
    are inferred from the supplied FASTA unless explicitly provided. They classify
    records but never filter quantified features. Unmatched peptides are retained.
    """
    from anndata_proteomics.annotation.validate_fasta import (
        validate_peptide_modalities_against_fasta,
    )
    from anndata_proteomics.annotation.var_fasta import annotate_var_from_fasta
    from anndata_proteomics.fasta.anndata_io import read_fasta_config
    from anndata_proteomics.readers.result import load_converted_result
    from anndata_proteomics.readers.summary import store_fasta_summary

    if not fasta_files:
        logger.error("no FASTA file given; usage: apb fasta DATA FASTA [FASTA ...]")
        return 1

    sources = list(fasta_files)
    obj = load_converted_result(data)
    has_protein = (
        "protein" in obj.mod if hasattr(obj, "mod") else _quantification_level(obj) == "protein"
    )
    peptide_modalities = (
        [
            name
            for name, target in obj.mod.items()
            if _quantification_level(target) in {"ion", "fragment", "peptidoform", "peptide"}
        ]
        if hasattr(obj, "mod")
        else (
            [_quantification_level(obj)]
            if _quantification_level(obj) in {"ion", "fragment", "peptidoform", "peptide"}
            else []
        )
    )

    if has_protein:
        annotate_var_from_fasta(
            obj,
            sources,
            match_on=match_on,
            is_uniprot=is_uniprot,
            decoy_pattern=decoy_pattern,
            contaminant_pattern=contaminant_pattern,
            cleavage=cleavage,
            min_length=min_length,
            max_length=max_length,
        )

    if validate and peptide_modalities:
        resolved_config = (
            read_fasta_config(obj)
            if has_protein and decoy_pattern is None and contaminant_pattern is None
            else None
        )
        results = validate_peptide_modalities_against_fasta(
            obj,
            sources,
            sequence_field=sequence_field,
            backend=backend,
            fasta_config=resolved_config,
            decoy_pattern=decoy_pattern,
            contaminant_pattern=contaminant_pattern,
            leading_protein_field=leading_protein_field,
            protein_match_on=match_on,
            il_equivalent=il_equivalent,
            is_uniprot=is_uniprot,
        )
        for name, result in results.items():
            logger.info(
                "{}: {}/{} peptide-derived features occur in FASTA",
                name,
                result.n_matched_features,
                result.n_features,
            )

    if not has_protein and (not validate or not peptide_modalities):
        logger.error(
            "input has no protein layer to annotate and no enabled peptide-derived "
            "layer to validate"
        )
        return 1

    store_fasta_summary(obj)
    out = output or data.with_name(f"{data.stem}.fasta{data.suffix}")
    if hasattr(obj, "mod"):
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
    """Compute ProteoBench scores without requiring annotation or FASTA stages.

    The module TOML supplies the sample design and species ratios. Vendor-specific
    interpretation is already complete in the converted APB object. The default
    output is ``<stem>.proteobench<suffix>`` beside the input.
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
