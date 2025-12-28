"""
ProForma modification conversion utilities.

Adapted from ProteoBench (proteobench/io/parsing/parse_ion.py)
https://github.com/Proteobench/ProteoBench
"""

import math
import re
import tomllib
from pathlib import Path
from typing import Optional


class ProFormaConverter:
    """
    Convert software-specific modification notation to ProForma format.

    Parameters
    ----------
    pattern : str
        Regex pattern to match modifications.
    modification_dict : dict
        Mapping from software notation to ProForma names.
    before_aa : bool
        Whether modifications appear before the amino acid in input.
    isalpha : bool
        Whether to filter for alphabetic characters.
    isupper : bool
        Whether to filter for uppercase characters.

    Example
    -------
    >>> converter = ProFormaConverter(
    ...     pattern=r"\\(([^()]*)\\)",
    ...     modification_dict={"(unimod:35)": "Oxidation"}
    ... )
    >>> converter.convert("PEPTM(unimod:35)IDE")
    'PEPTM[Oxidation]IDE'
    """

    def __init__(
        self,
        pattern: str = r"\[([^]]+)\]",
        modification_dict: Optional[dict] = None,
        before_aa: bool = False,
        isalpha: bool = True,
        isupper: bool = True,
    ):
        self.pattern = pattern
        self.modification_dict = modification_dict or {
            "+57.0215": "Carbamidomethyl",
            "+15.9949": "Oxidation",
            "-17.026548": "Gln->pyro-Glu",
            "-18.010565": "Glu->pyro-Glu",
            "+42": "Acetyl",
        }
        self.before_aa = before_aa
        self.isalpha = isalpha
        self.isupper = isupper

    @classmethod
    def from_config(cls, config: dict) -> "ProFormaConverter":
        """Create converter from a modifications_parser config dict."""
        return cls(
            pattern=config.get("pattern", r"\[([^]]+)\]"),
            modification_dict=config.get("modification_dict", {}),
            before_aa=config.get("before_aa", False),
            isalpha=config.get("isalpha", True),
            isupper=config.get("isupper", True),
        )

    @classmethod
    def from_toml(cls, path: str | Path) -> "ProFormaConverter":
        """
        Load ProForma converter from a TOML configuration file.

        Parameters
        ----------
        path : str or Path
            Path to the TOML configuration file.

        Returns
        -------
        ProFormaConverter
            Configured converter instance.

        Example
        -------
        >>> converter = ProFormaConverter.from_toml("configs/diann.toml")
        >>> converter.convert("PEPTM(unimod:35)IDE")
        'PEPTM[Oxidation]IDE'
        """
        path = Path(path)
        with open(path, "rb") as f:
            config = tomllib.load(f)

        # Extract patterns section
        patterns = config.get("patterns", {})

        # Build the config dict for from_config
        converter_config = {
            "pattern": patterns.get("pattern", r"\[([^]]+)\]"),
            "modification_dict": config.get("modifications", {}),
            "before_aa": patterns.get("before_aa", False),
            "isalpha": patterns.get("isalpha", True),
            "isupper": patterns.get("isupper", True),
        }

        return cls.from_config(converter_config)

    def _count_chars(self, s: str) -> int:
        """Count characters matching filter criteria."""
        if self.isalpha and self.isupper:
            return sum(1 for c in s if c.isalpha() and c.isupper())
        if self.isalpha:
            return sum(1 for c in s if c.isalpha())
        if self.isupper:
            return sum(1 for c in s if c.isupper())
        return len(s)

    def get_stripped_sequence(self, s: str) -> str:
        """Extract bare sequence without modifications."""
        if self.isalpha and self.isupper:
            return "".join(c for c in s if c.isalpha() and c.isupper())
        if self.isalpha:
            return "".join(c for c in s if c.isalpha())
        if self.isupper:
            return "".join(c for c in s if c.isupper())
        return s

    def _match_modifications(self, s: str) -> tuple[list[str], list[int]]:
        """Find modifications and their positions in the sequence."""
        matches = [(m.group(), m.start(), m.end()) for m in re.finditer(self.pattern, s)]
        positions = [self._count_chars(s[:m[1]]) for m in matches]
        mods = [m[0] for m in matches]
        return mods, positions

    def _map_modification(self, mod: str) -> str:
        """Map software-specific modification to ProForma name."""
        mod_lower = mod.lower()
        return self.modification_dict.get(mod_lower, self.modification_dict.get(mod, mod.strip("[](){}")))

    def convert(self, sequence: str) -> str:
        """
        Convert sequence to ProForma notation.

        Parameters
        ----------
        sequence : str
            Input sequence with software-specific modification notation.

        Returns
        -------
        str
            Sequence in ProForma notation.
        """
        if not sequence or (isinstance(sequence, float) and math.isnan(sequence)):
            return ""

        # Strip common wrappers (e.g., MaxQuant underscores)
        sequence = str(sequence).strip("_")

        # Lowercase for matching
        seq_lower = re.sub(self.pattern, lambda m: m.group(0).lower(), sequence)

        # Find modifications
        mods, positions = self._match_modifications(seq_lower)
        mapped_mods = [self._map_modification(m) for m in mods]
        pos_mod_dict = dict(zip(positions, mapped_mods))

        # Build ProForma sequence
        stripped = self.get_stripped_sequence(seq_lower)
        result = []

        for idx, aa in enumerate(stripped):
            if self.before_aa:
                result.append(aa)

            if idx in pos_mod_dict:
                mod = pos_mod_dict[idx]
                if idx == 0:
                    result.append(f"[{mod}]-")
                elif idx == len(stripped):
                    result.append(f"-[{mod}]")
                else:
                    result.append(f"[{mod}]")

            if not self.before_aa:
                result.append(aa)

        return "".join(result)

    def convert_with_sites(self, sequence: str, modifications: str, sites: str) -> str:
        """
        Convert sequence with separate modification/site columns (AlphaDIA format).

        Parameters
        ----------
        sequence : str
            The bare sequence.
        modifications : str
            Semicolon-separated modification names.
        sites : str
            Semicolon-separated modification positions.

        Returns
        -------
        str
            Sequence in ProForma notation.
        """
        if isinstance(modifications, float) and math.isnan(modifications):
            return sequence
        if not modifications:
            return sequence

        mods_list = modifications.split(";")
        sites_list = list(map(int, str(sites).split(";")))

        # Sort by position (reverse) to insert from end
        mods_and_sites = sorted(zip(mods_list, sites_list), key=lambda x: x[1], reverse=True)

        result = sequence
        for mod, site in mods_and_sites:
            if not mod:
                continue
            mod_name = mod.split("@")[0]
            if site == 0:
                result = f"[{mod_name}]-" + result
            elif site == -1:
                result = result + f"-[{mod_name}]"
            else:
                result = result[:site] + f"[{mod_name}]" + result[site:]

        return result


# Convenience function for simple use cases
def to_proforma(
    sequence: str,
    pattern: str = r"\[([^]]+)\]",
    modification_dict: Optional[dict] = None,
) -> str:
    """
    Convert sequence to ProForma notation (convenience function).

    Parameters
    ----------
    sequence : str
        Input sequence with modifications.
    pattern : str
        Regex pattern to match modifications.
    modification_dict : dict, optional
        Mapping from software notation to ProForma names.

    Returns
    -------
    str
        Sequence in ProForma notation.
    """
    converter = ProFormaConverter(pattern=pattern, modification_dict=modification_dict or {})
    return converter.convert(sequence)
