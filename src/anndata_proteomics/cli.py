"""
Command-line interface for anndata_proteomics.

Usage:
    prot2ad convert report.tsv output.h5ad
    prot2ad full report.tsv annotation.csv params.log output.h5ad
"""

import sys
from pathlib import Path
from typing import Annotated, Optional

from cyclopts import App, Parameter

from .builder import ConverterBuilder

app = App(
    name="prot2ad",
    help="Convert proteomics software output to AnnData format.",
)


# Valid levels for validation
VALID_LEVELS = ["ion", "peptidoform", "peptide", "protein"]


def _validate_level(detected_level: str, specified_level: Optional[str]) -> None:
    """Validate that specified level matches detected level."""
    if specified_level is None:
        return

    if specified_level not in VALID_LEVELS:
        print(f"Error: Invalid level '{specified_level}'. Must be one of: {VALID_LEVELS}")
        sys.exit(1)

    if specified_level != detected_level:
        print(
            f"Error: Specified level '{specified_level}' does not match "
            f"detected level '{detected_level}'"
        )
        sys.exit(1)


def _get_converter(quant_file: Path, software: Optional[str]):
    """Get converter, either auto-detected or for specified software."""
    if software:
        try:
            converter = ConverterBuilder.for_software(software)
            # Also check if file matches
            auto_converter = ConverterBuilder.from_file(quant_file)
            if auto_converter.strategy.name.lower() != software.lower():
                print(
                    f"Warning: Specified software '{software}' but file appears to be "
                    f"'{auto_converter.strategy.name}'"
                )
        except ValueError as e:
            print(f"Error: {e}")
            sys.exit(1)
    else:
        try:
            converter = ConverterBuilder.from_file(quant_file)
        except ValueError as e:
            print(f"Error: Could not detect file format. {e}")
            print("Try specifying --software explicitly.")
            sys.exit(1)

    return converter


@app.command
def convert(
    quant_file: Annotated[Path, Parameter(help="Proteomics quantification file")],
    output: Annotated[Path, Parameter(help="Output .h5ad file path")],
    *,
    software: Annotated[
        Optional[str],
        Parameter(name=["--software", "-s"], help="Override software detection"),
    ] = None,
    level: Annotated[
        Optional[str],
        Parameter(name=["--level", "-l"], help="Validate quantification level"),
    ] = None,
) -> None:
    """
    Basic conversion: quantification file to AnnData.

    Minimal conversion without sample annotation or parameter parsing.
    Sample metadata will be auto-generated from file names.
    """
    if not quant_file.exists():
        print(f"Error: File not found: {quant_file}")
        sys.exit(1)

    converter = _get_converter(quant_file, software)
    strategy = converter.strategy

    print(f"Detected: {strategy.name} ({strategy.level} level)")

    # Validate level if specified
    _validate_level(strategy.level, level)

    # For basic convert, we don't have annotation - create minimal obs
    print(f"Loading {quant_file}...")
    df = strategy.load(quant_file)

    # Get var and layers
    var_df = strategy.get_var(df)
    layers_df = strategy.get_layers(df)

    # Pivot to create matrices
    from .core import Converter

    temp_converter = Converter(strategy)

    # Create a minimal annotation from unique obs values
    import pandas as pd

    unique_obs = df[strategy.obs_id].unique()
    annotation_df = pd.DataFrame({"sample": unique_obs})
    annotation_df.to_csv("/tmp/_temp_annotation.csv", index=False)

    try:
        adata = temp_converter.convert(
            quant_file,
            "/tmp/_temp_annotation.csv",
            annotation_id_col="sample",
        )
    finally:
        Path("/tmp/_temp_annotation.csv").unlink(missing_ok=True)

    # Save output
    print(f"Writing {output}...")
    adata.write_h5ad(output)
    print(f"Done! Created {output} with {adata.n_obs} samples x {adata.n_vars} precursors")


@app.command
def full(
    quant_file: Annotated[Path, Parameter(help="Proteomics quantification file")],
    annotation: Annotated[Path, Parameter(help="Sample annotation CSV file")],
    params: Annotated[Path, Parameter(help="Parameter/settings file")],
    output: Annotated[Path, Parameter(help="Output .h5ad file path")],
    *,
    software: Annotated[
        Optional[str],
        Parameter(name=["--software", "-s"], help="Override software detection"),
    ] = None,
    level: Annotated[
        Optional[str],
        Parameter(name=["--level", "-l"], help="Validate quantification level"),
    ] = None,
    annotation_id_col: Annotated[
        str,
        Parameter(name="--annotation-id-col", help="Column for sample IDs in annotation"),
    ] = "sample",
) -> None:
    """
    Full conversion: quantification + annotation + parameters.

    Complete conversion with sample annotation and parameter extraction.
    """
    # Check files exist
    for path, name in [
        (quant_file, "Quantification"),
        (annotation, "Annotation"),
        (params, "Parameters"),
    ]:
        if not path.exists():
            print(f"Error: {name} file not found: {path}")
            sys.exit(1)

    converter = _get_converter(quant_file, software)
    strategy = converter.strategy

    print(f"Detected: {strategy.name} ({strategy.level} level)")

    # Validate level if specified
    _validate_level(strategy.level, level)

    # Parse parameters (placeholder)
    print(f"Parsing parameters from {params}...")
    from .params.diann import extract_params as extract_diann_params
    from .params.maxquant import extract_params as extract_maxquant_params
    from .params.spectronaut import extract_params as extract_spectronaut_params

    param_extractors = {
        "DIA-NN": extract_diann_params,
        "MaxQuant": extract_maxquant_params,
        "Spectronaut": extract_spectronaut_params,
    }

    extractor = param_extractors.get(strategy.name)
    if extractor:
        parsed_params = extractor(params)
        print(f"  Parameters: {parsed_params}")
    else:
        print(f"  Warning: No parameter parser for {strategy.name}")
        parsed_params = None

    # Convert
    print(f"Loading {quant_file}...")
    adata = converter.convert(
        quant_file,
        annotation,
        annotation_id_col=annotation_id_col,
    )

    # Store parsed parameters in uns if available
    if parsed_params:
        adata.uns["params"] = {
            k: v for k, v in parsed_params.__dict__.items() if v is not None
        }

    # Save output
    print(f"Writing {output}...")
    adata.write_h5ad(output)
    print(f"Done! Created {output} with {adata.n_obs} samples x {adata.n_vars} precursors")


def main():
    """Entry point for CLI."""
    app()


if __name__ == "__main__":
    main()
