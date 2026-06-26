"""apb CLI dispatcher.

Subcommands:
- validate [path ...]        validate one or more TOML rules; defaults to all packaged
- list                       list packaged rules
- export-schema              regenerate parse_rule.schema.json
- convert <data> [toml]      convert a vendor file to AnnData (.h5ad)
- annotate <data> <toml>     join sample annotations onto obs
- fasta <data> <fasta...>    annotate the protein layer's var from FASTA file(s)
"""

from __future__ import annotations

import sys
from pathlib import Path

from cyclopts import App
from loguru import logger

from anndata_proteomics._logging import configure_default_sink
from anndata_proteomics.converters.assemble import convert as _run_convert
from anndata_proteomics.converters.recognize import recognize
from anndata_proteomics.readers.dispatch import read_table
from anndata_proteomics.rules import _export_schema
from anndata_proteomics.rules.loader import load_rule
from anndata_proteomics.rules.registry import iter_packaged_rules
from anndata_proteomics.rules.validate import (
    _log_and_exit_code,
    validate_all_packaged,
    validate_file,
)

app = App(name="apb", help="anndata_proteomics (APB) CLI")


@app.command
def validate(*paths: Path) -> int:
    """Validate one or more TOML rule files.

    With no paths, walks all packaged rules (same as `validate-rules`).
    """
    if not paths:
        results = validate_all_packaged()
    else:
        results = [validate_file(p) for p in paths]
    return _log_and_exit_code(results)


@app.command(name="list")
def list_rules() -> int:
    """List packaged parsing rules: software, level, file_version, path."""
    for p in iter_packaged_rules():
        parts = p.stem.split("_")
        # filename: parse_<software_tokens>_<level>_<version>.toml
        software = "_".join(parts[1:-2])
        level = parts[-2]
        version = parts[-1]
        logger.info(f"{software:14}  {level:12}  v{version:<3}  {p}")
    return 0


@app.command(name="export-schema")
def export_schema_cmd() -> int:
    """Regenerate parse_rule.schema.json from the pydantic models."""
    _export_schema.main()
    return 0


@app.command
def convert(
    data: Path,
    rule_toml: Path | None = None,
    output: Path | None = None,
) -> int:
    """Convert a vendor file to AnnData and write a .h5ad.

    If --rule-toml is omitted, the rule is auto-recognized from the data's
    column headers. Use --rule-toml to override or when recognition is ambiguous.
    --output defaults to <data>.h5ad next to the input.
    """
    df = read_table(data)
    if rule_toml is None:
        rule = recognize(list(df.columns))
        if rule is None:
            logger.error(f"could not auto-recognize a rule for {data}; pass --rule-toml PATH")
            return 1
    else:
        rule = load_rule(rule_toml)
    adata = _run_convert(df, rule)
    out = output or data.with_suffix(".h5ad")
    adata.write_h5ad(out)
    logger.info(f"wrote {out}  shape={adata.shape}  layers={list(adata.layers)}")
    return 0


@app.command
def annotate(
    data: Path,
    annotation_toml: Path,
    output: Path | None = None,
) -> int:
    """Join sample annotations from an annotation TOML onto obs and write the result.

    Reads an .h5ad/.h5mu, joins the TOML's ``obs.samples`` table onto obs by run/file
    name (matching obs_names by default; see the TOML's ``match_on``), and writes the
    enriched object. --output defaults to ``<stem>.annotated<suffix>`` next to the input
    (non-destructive); point it back at the input to update in place.
    """
    from anndata_proteomics.annotation.apply import annotate_obs
    from anndata_proteomics.annotation.loader import load_annotation
    from anndata_proteomics.scripts._ui_support import load_converted_result

    obj = load_converted_result(data)
    spec = load_annotation(annotation_toml)
    annotate_obs(obj, spec)

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
    output: Path | None = None,
    match_on: str = "Protein_Group",
    is_uniprot: bool = True,
    decoy_pattern: str = "^REV_|^rev_",
    cleavage: str | None = None,
    min_length: int | None = None,
    max_length: int | None = None,
) -> int:
    """Annotate the protein layer from FASTA file(s) and write the result.

    Reads an .h5ad/.h5mu, builds the prolfquapp-style protein annotation
    (fasta.id, fasta.header, protein_length, nr_peptides, gene_name) from
    one or more FASTA files, and attaches it as a var-aligned DataFrame at
    ``varm['fasta']`` of the protein layer only (a protein-level AnnData, or the
    ``protein`` modality of a MuData). The join matches the leading accession of
    each protein group against the FASTA proteinname. ``nr_peptides`` uses
    the digestion enzyme stored in the object's search parameters;
    --cleavage / --min-length / --max-length override it (needed for objects
    converted without a parameters file). --output defaults to
    ``<stem>.annotated<suffix>`` (non-destructive).
    """
    from anndata_proteomics.annotation.var_fasta import annotate_var_from_fasta
    from anndata_proteomics.scripts._ui_support import load_converted_result

    if not fasta_files:
        logger.error("no FASTA file given; usage: apb fasta DATA FASTA [FASTA ...]")
        return 1

    obj = load_converted_result(data)
    annotate_var_from_fasta(
        obj,
        list(fasta_files),
        match_on=match_on,
        is_uniprot=is_uniprot,
        decoy_pattern=decoy_pattern,
        cleavage=cleavage,
        min_length=min_length,
        max_length=max_length,
    )

    out = output or data.with_name(f"{data.stem}.annotated{data.suffix}")
    if hasattr(obj, "mod"):
        obj.write_h5mu(out)
    else:
        obj.write_h5ad(out)
    logger.info(f"wrote {out}")
    return 0


def main() -> int:
    """Console-script entry point."""
    configure_default_sink()
    rc = app()
    return int(rc) if rc is not None else 0


if __name__ == "__main__":
    sys.exit(main())
