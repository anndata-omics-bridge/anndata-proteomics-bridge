"""
Base parameter dataclass for proteomics software.

Copied from ProteoBench (https://github.com/Proteobench/proteobench).
"""

from dataclasses import dataclass
from typing import Optional


@dataclass
class ProteoBenchParameters:
    """
    Standardized parameter container for proteomics software settings.

    Attributes match ProteoBench convention for interoperability.
    """

    # Software identification
    software_name: Optional[str] = None
    software_version: Optional[str] = None
    search_engine: Optional[str] = None
    search_engine_version: Optional[str] = None

    # FDR thresholds
    ident_fdr_psm: Optional[float] = None
    ident_fdr_peptide: Optional[float] = None
    ident_fdr_protein: Optional[float] = None

    # Mass tolerances
    precursor_mass_tolerance: Optional[float] = None
    fragment_mass_tolerance: Optional[float] = None

    # Enzyme settings
    enzyme: Optional[str] = None
    allowed_miscleavages: Optional[int] = None
    semi_enzymatic: Optional[bool] = None

    # Peptide length
    min_peptide_length: Optional[int] = None
    max_peptide_length: Optional[int] = None

    # Modifications
    fixed_mods: Optional[str] = None
    variable_mods: Optional[str] = None
    max_mods: Optional[int] = None

    # Charge range
    min_precursor_charge: Optional[int] = None
    max_precursor_charge: Optional[int] = None

    # Quantification
    enable_match_between_runs: Optional[bool] = None
    quantification_method: Optional[str] = None
    protein_inference: Optional[str] = None

    def __repr__(self) -> str:
        """Show only non-None attributes."""
        attrs = {k: v for k, v in self.__dict__.items() if v is not None}
        return f"ProteoBenchParameters({attrs})"
