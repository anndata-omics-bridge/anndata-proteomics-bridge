"""
Proteomics to AnnData Converter

Convert proteomics software outputs (DIA-NN, Spectronaut, etc.) to AnnData format
following the AnnData Omics Bridge specification.
"""

__version__ = "0.1.0"

from .converter import proteomics_to_anndata, convert_file
from .reader import load_file, load_diann, load_spectronaut
from .annotation import load_annotation, match_samples
from .parse_settings import load_config

__all__ = [
    "proteomics_to_anndata",
    "convert_file",
    "load_file",
    "load_diann",
    "load_spectronaut",
    "load_annotation",
    "match_samples",
    "load_config",
]
