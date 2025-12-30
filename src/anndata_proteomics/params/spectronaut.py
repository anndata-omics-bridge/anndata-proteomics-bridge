"""
Spectronaut parameter extraction.

Placeholder implementation - full parsing to be implemented.
See ProteoBench for reference implementation.
"""

from pathlib import Path
from typing import Union

from .parameters import ProteoBenchParameters


def extract_params(path: Union[str, Path]) -> ProteoBenchParameters:
    """
    Extract parameters from Spectronaut settings export file.

    Parameters
    ----------
    path : str or Path
        Path to Spectronaut settings file (.txt)

    Returns
    -------
    ProteoBenchParameters
        Extracted parameters (placeholder: only software_name set)

    Notes
    -----
    TODO: Implement full parsing from Spectronaut settings exports.
    Reference: ProteoBench proteobench/io/params/spectronaut.py
    """
    return ProteoBenchParameters(software_name="Spectronaut")
