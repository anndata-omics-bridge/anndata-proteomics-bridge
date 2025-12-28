"""
ConverterBuilder for auto-detecting file formats and creating converters.

Usage:
    # Auto-detect format
    converter = ConverterBuilder.from_file("evidence.txt")
    adata = converter.convert("evidence.txt", "annotation.csv")

    # Explicit software
    converter = ConverterBuilder.for_software("MaxQuant")
    adata = converter.convert("evidence.txt", "annotation.csv")
"""

from pathlib import Path

from .core import Converter
from .strategies.diann import DIANNStrategy
from .strategies.maxquant import MaxQuantStrategy
from .strategies.spectronaut import SpectronautStrategy


def _normalize_software_name(name: str) -> str:
    """Normalize software name for registry lookup."""
    return name.lower().replace("-", "").replace(" ", "")


# Registry of available strategies
STRATEGY_REGISTRY = {
    _normalize_software_name(DIANNStrategy.name): DIANNStrategy,
    _normalize_software_name(MaxQuantStrategy.name): MaxQuantStrategy,
    _normalize_software_name(SpectronautStrategy.name): SpectronautStrategy,
}

# Strategies to try for auto-detection
AUTO_DETECT_STRATEGIES = [cls() for cls in STRATEGY_REGISTRY.values()]


class ConverterBuilder:
    """
    Builder for creating Converter instances.

    Supports auto-detection of file formats and explicit software selection.
    """

    @classmethod
    def from_file(cls, path: str | Path) -> Converter:
        """
        Create a Converter by auto-detecting the file format.

        Parameters
        ----------
        path : str or Path
            Path to the proteomics data file.

        Returns
        -------
        Converter
            Configured converter for the detected format.

        Raises
        ------
        ValueError
            If the file format cannot be detected.

        Examples
        --------
        >>> converter = ConverterBuilder.from_file("evidence.txt")
        >>> adata = converter.convert("evidence.txt", "annotation.csv")
        """
        path = Path(path)

        for strategy in AUTO_DETECT_STRATEGIES:
            if strategy.detect(path):
                # Create fresh instance
                strategy_class = type(strategy)
                strategy = strategy_class()
                print(f"Detected format: {strategy.name}")
                return Converter(strategy)

        raise ValueError(
            f"Could not detect file format for: {path}\n"
            f"Supported formats: {list(STRATEGY_REGISTRY.keys())}"
        )

    @classmethod
    def for_software(cls, software: str) -> Converter:
        """
        Create a Converter for a specific software.

        Parameters
        ----------
        software : str
            Software name (e.g., "MaxQuant", "DIA-NN", "Spectronaut").

        Returns
        -------
        Converter
            Configured converter for the specified software.

        Raises
        ------
        ValueError
            If the software is not supported.

        Examples
        --------
        >>> converter = ConverterBuilder.for_software("MaxQuant")
        >>> adata = converter.convert("evidence.txt", "annotation.csv")
        """
        normalized = _normalize_software_name(software)
        if normalized not in STRATEGY_REGISTRY:
            raise ValueError(
                f"Unknown software: {software}\n"
                f"Supported: {list(STRATEGY_REGISTRY.keys())}"
            )

        strategy_class = STRATEGY_REGISTRY[normalized]
        return Converter(strategy_class())

    @classmethod
    def list_supported(cls) -> list[str]:
        """
        List supported software formats.

        Returns
        -------
        list of str
            Names of supported software formats.
        """
        return list(STRATEGY_REGISTRY.keys())
