"""Typed parameter models for proteomics search-engine settings."""

from __future__ import annotations

import math
import re
from typing import Literal

import pandas as pd
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    NonNegativeFloat,
    NonNegativeInt,
    PositiveInt,
    field_validator,
    model_validator,
)

from anndata_proteomics.modifications.model import SearchedModification

ScalarValue = str | int | float | bool | None
ToleranceUnit = Literal["ppm", "Da"]
ToleranceMode = Literal["absolute", "range", "automatic"]

_MISSING_STRINGS = {"", "none", "nan", "n/a", "na", "not specified", "unknown"}
_RANGE_RE = re.compile(
    r"^\[\s*(?P<lower>[+-]?\d+(?:\.\d+)?)\s*(?P<unit1>[A-Za-z]*)\s*,\s*"
    r"(?P<upper>[+-]?\d+(?:\.\d+)?)\s*(?P<unit2>[A-Za-z]*)\s*\]$"
)
_ABSOLUTE_RE = re.compile(r"^(?P<value>[+-]?\d+(?:\.\d+)?)\s*(?P<unit>[A-Za-z]*)$")


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid")


class Probability(_Strict):
    """A probability value constrained to the closed interval [0, 1]."""

    value: float = Field(ge=0, le=1)

    @field_validator("value", mode="before")
    @classmethod
    def _coerce_value(cls, value: object) -> object:
        if isinstance(value, str):
            text = value.strip()
            if text.endswith("%"):
                return float(text[:-1]) / 100
            return float(text)
        return value

    def to_legacy(self) -> float:
        """Return the scalar representation used in ProteoBench CSV fixtures."""
        return self.value


class MassTolerance(_Strict):
    """Mass tolerance as an absolute value, signed range, or automatic mode."""

    lower: float | None = None
    upper: float | None = None
    value: NonNegativeFloat | None = None
    unit: ToleranceUnit | None = None
    mode: ToleranceMode
    label: str | None = None

    @model_validator(mode="after")
    def _validate_shape(self) -> "MassTolerance":
        if self.mode == "absolute":
            if self.value is None:
                raise ValueError("absolute tolerance requires value")
            if self.unit is None:
                raise ValueError("absolute tolerance requires unit")
            if self.lower is not None or self.upper is not None:
                raise ValueError("absolute tolerance cannot also define lower/upper")
        elif self.mode == "range":
            if self.lower is None or self.upper is None:
                raise ValueError("range tolerance requires lower and upper")
            if self.unit is None:
                raise ValueError("range tolerance requires unit")
            if self.lower > self.upper:
                raise ValueError("range tolerance requires lower <= upper")
            if self.value is not None:
                raise ValueError("range tolerance cannot also define value")
        elif self.mode == "automatic":
            if self.unit is not None:
                raise ValueError("automatic tolerance cannot define unit")
            if self.value is not None or self.lower is not None or self.upper is not None:
                raise ValueError("automatic tolerance cannot define numeric bounds")
        return self

    @classmethod
    def parse(cls, value: object) -> "MassTolerance | None":
        """Parse vendor/legacy tolerance values into a typed tolerance."""
        if _is_missing(value):
            return None
        if isinstance(value, MassTolerance):
            return value
        if isinstance(value, int | float):
            raise ValueError("mass tolerance numeric values require an explicit unit")
        if not isinstance(value, str):
            raise TypeError(f"unsupported mass tolerance value: {value!r}")

        text = value.strip()
        if text.lower() in {
            "dynamic",
            "automatic",
            "automatic calibration",
            "auto",
            "auto detected",
        }:
            return cls(mode="automatic", label=text)

        range_match = _RANGE_RE.match(text)
        if range_match:
            unit = _normalize_unit(
                range_match.group("unit1") or range_match.group("unit2")
            )
            return cls(
                mode="range",
                lower=float(range_match.group("lower")),
                upper=float(range_match.group("upper")),
                unit=unit,
            )

        absolute_match = _ABSOLUTE_RE.match(text)
        if absolute_match:
            unit = _normalize_unit(absolute_match.group("unit"))
            return cls(
                mode="absolute",
                value=float(absolute_match.group("value")),
                unit=unit,
            )

        raise ValueError(f"could not parse mass tolerance: {value!r}")

    def to_legacy(self) -> str:
        """Return the scalar representation used in ProteoBench CSV fixtures."""
        if self.mode == "automatic":
            return self.label or "Dynamic"
        if self.mode == "absolute":
            return _format_number(self.value) + _format_unit(self.unit)
        return (
            f"[{_format_number(self.lower)}{_format_unit(self.unit)}, "
            f"{_format_number(self.upper)}{_format_unit(self.unit)}]"
        )


class ChargeRange(_Strict):
    """Optional precursor charge range."""

    minimum: PositiveInt | None = None
    maximum: PositiveInt | None = None

    @model_validator(mode="after")
    def _ordered(self) -> "ChargeRange":
        if self.minimum is not None and self.maximum is not None:
            if self.minimum > self.maximum:
                raise ValueError("minimum charge cannot exceed maximum charge")
        return self


class MzRange(_Strict):
    """Optional m/z range."""

    minimum: NonNegativeFloat | None = None
    maximum: NonNegativeFloat | None = None

    @model_validator(mode="after")
    def _ordered(self) -> "MzRange":
        if self.minimum is not None and self.maximum is not None:
            if self.minimum > self.maximum:
                raise ValueError("minimum m/z cannot exceed maximum m/z")
        return self


class UnparsedParameter(_Strict):
    """Explicitly retained vendor value that is outside the core schema."""

    name: str
    value: ScalarValue
    source: str | None = None


class Parameters(_Strict):
    """Proteomics search-parameter record with typed fields."""

    software_name: str | None = None
    software_version: str | None = None
    search_engine: str | None = None
    search_engine_version: str | None = None
    ident_fdr_psm: Probability | None = None
    ident_fdr_peptide: Probability | None = None
    ident_fdr_protein: Probability | None = None
    enable_match_between_runs: bool | None = None
    precursor_mass_tolerance: MassTolerance | None = None
    fragment_mass_tolerance: MassTolerance | None = None
    enzyme: str | None = None
    semi_enzymatic: bool | None = None
    allowed_miscleavages: NonNegativeInt | None = None
    min_peptide_length: NonNegativeInt | None = None
    max_peptide_length: NonNegativeInt | None = None
    fixed_mods: list[SearchedModification] = Field(default_factory=list)
    variable_mods: list[SearchedModification] = Field(default_factory=list)
    max_mods: NonNegativeInt | None = None
    min_precursor_charge: PositiveInt | None = None
    max_precursor_charge: PositiveInt | None = None
    min_precursor_mz: NonNegativeFloat | None = None
    max_precursor_mz: NonNegativeFloat | None = None
    min_fragment_mz: NonNegativeFloat | None = None
    max_fragment_mz: NonNegativeFloat | None = None
    quantification_method: str | None = None
    protein_inference: str | None = None
    abundance_normalization_ions: str | bool | None = None
    predictors_library: str | None = None
    scan_window: NonNegativeInt | str | None = None
    unparsed_parameters: list[UnparsedParameter] = Field(default_factory=list)

    @field_validator(
        "software_name",
        "software_version",
        "search_engine",
        "search_engine_version",
        "enzyme",
        "quantification_method",
        "protein_inference",
        "abundance_normalization_ions",
        "predictors_library",
        mode="before",
    )
    @classmethod
    def _empty_strings_to_none(cls, value: object) -> object:
        return None if _is_missing(value) else value

    @field_validator(
        "ident_fdr_psm",
        "ident_fdr_peptide",
        "ident_fdr_protein",
        mode="before",
    )
    @classmethod
    def _coerce_probability(cls, value: object) -> object:
        if _is_missing(value):
            return None
        if isinstance(value, Probability):
            return value
        numeric = _coerce_float(value)
        if numeric is None:
            return None
        if numeric > 1:
            numeric /= 100
        return Probability(value=numeric)

    @field_validator("precursor_mass_tolerance", "fragment_mass_tolerance", mode="before")
    @classmethod
    def _coerce_tolerance(cls, value: object) -> object:
        return MassTolerance.parse(value)

    @field_validator(
        "allowed_miscleavages",
        "min_peptide_length",
        "max_peptide_length",
        "max_mods",
        mode="before",
    )
    @classmethod
    def _coerce_non_negative_int(cls, value: object) -> object:
        if _is_missing(value):
            return None
        return int(float(str(value).strip()))

    @field_validator("min_precursor_charge", "max_precursor_charge", mode="before")
    @classmethod
    def _coerce_positive_int(cls, value: object) -> object:
        if _is_missing(value):
            return None
        return int(float(str(value).strip()))

    @field_validator(
        "min_precursor_mz",
        "max_precursor_mz",
        "min_fragment_mz",
        "max_fragment_mz",
        mode="before",
    )
    @classmethod
    def _coerce_non_negative_float(cls, value: object) -> object:
        if _is_missing(value):
            return None
        return float(str(value).strip())

    @field_validator("enable_match_between_runs", "semi_enzymatic", mode="before")
    @classmethod
    def _coerce_bool(cls, value: object) -> object:
        if _is_missing(value):
            return None
        if isinstance(value, bool):
            return value
        if isinstance(value, int | float):
            return bool(value)
        if isinstance(value, str):
            text = value.strip().lower()
            if text in {"true", "1", "yes", "y"}:
                return True
            if text in {"false", "0", "no", "n"}:
                return False
        raise ValueError(f"cannot coerce boolean value: {value!r}")

    @field_validator("scan_window", mode="before")
    @classmethod
    def _coerce_scan_window(cls, value: object) -> object:
        if _is_missing(value):
            return None
        if isinstance(value, int):
            return value
        if isinstance(value, float):
            return int(value) if value.is_integer() else str(value)
        if isinstance(value, str):
            text = value.strip()
            try:
                numeric = float(text)
            except ValueError:
                return text
            return int(numeric) if numeric.is_integer() else text
        return value

    @field_validator("fixed_mods", "variable_mods", mode="before")
    @classmethod
    def _coerce_modifications(cls, value: object) -> object:
        if _is_missing(value):
            return []
        if isinstance(value, SearchedModification):
            return [value]
        if isinstance(value, dict):
            return [
                SearchedModification(
                    name=str(target),
                    target=str(target),
                    mass_delta=_coerce_float(delta),
                    source=f"{target}: {delta}",
                )
                for target, delta in value.items()
            ]
        if isinstance(value, list | tuple | set):
            return [_modification_from_item(item) for item in value if not _is_missing(item)]
        if isinstance(value, str):
            return [_modification_from_item(part) for part in _split_mod_string(value)]
        raise TypeError(f"unsupported modification value: {value!r}")

    @field_validator("unparsed_parameters", mode="before")
    @classmethod
    def _coerce_unparsed(cls, value: object) -> object:
        if _is_missing(value):
            return []
        return value

    @model_validator(mode="after")
    def _validate_ranges(self) -> "Parameters":
        _validate_order(self.min_precursor_charge, self.max_precursor_charge, "charge")
        _validate_order(self.min_peptide_length, self.max_peptide_length, "peptide length")
        _validate_order(self.min_precursor_mz, self.max_precursor_mz, "precursor m/z")
        _validate_order(self.min_fragment_mz, self.max_fragment_mz, "fragment m/z")
        return self

    def to_series(self) -> pd.Series:
        """Serialize to a pandas Series matching ProteoBench's CSV layout."""
        return pd.Series({field: self._legacy_value(field) for field in _SERIES_FIELDS})

    @classmethod
    def from_series(cls, series: pd.Series) -> "Parameters":
        """Build a ``Parameters`` instance from a ProteoBench-style Series."""
        fields = set(cls.model_fields)
        data: dict[str, object] = {}
        unparsed: list[UnparsedParameter] = []
        for key, value in series.items():
            name = str(key)
            normalized = None if _is_missing(value) else value
            if name in fields:
                data[name] = normalized
            else:
                unparsed.append(
                    UnparsedParameter(name=name, value=_to_scalar(normalized), source="series")
                )
        if unparsed:
            data["unparsed_parameters"] = unparsed
        return cls(**data)

    def _legacy_value(self, field: str) -> ScalarValue:
        value = getattr(self, field)
        if isinstance(value, Probability):
            return value.to_legacy()
        if isinstance(value, MassTolerance):
            return value.to_legacy()
        if isinstance(value, list):
            return _serialize_modifications(value)
        return _to_scalar(value)


_SERIES_FIELDS = (
    "software_name",
    "software_version",
    "search_engine",
    "search_engine_version",
    "ident_fdr_psm",
    "ident_fdr_peptide",
    "ident_fdr_protein",
    "enable_match_between_runs",
    "precursor_mass_tolerance",
    "fragment_mass_tolerance",
    "enzyme",
    "semi_enzymatic",
    "allowed_miscleavages",
    "min_peptide_length",
    "max_peptide_length",
    "fixed_mods",
    "variable_mods",
    "max_mods",
    "min_precursor_charge",
    "max_precursor_charge",
    "min_precursor_mz",
    "max_precursor_mz",
    "min_fragment_mz",
    "max_fragment_mz",
    "quantification_method",
    "protein_inference",
    "abundance_normalization_ions",
    "predictors_library",
    "scan_window",
)


def _is_missing(value: object) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    if isinstance(value, str) and value.strip().lower() in _MISSING_STRINGS:
        return True
    return False


def _normalize_unit(unit: str | None) -> ToleranceUnit:
    if not unit:
        raise ValueError("mass tolerance requires unit ppm or Da")
    lookup: dict[str, ToleranceUnit] = {"ppm": "ppm", "da": "Da", "th": "Da"}
    normalized = lookup.get(unit.strip().lower())
    if normalized is None:
        raise ValueError("mass tolerance unit must be ppm or Da")
    return normalized


def _format_unit(unit: ToleranceUnit | None) -> str:
    if unit is None:
        raise ValueError("cannot format mass tolerance without unit")
    return f" {unit}"


def _format_number(value: float | None) -> str:
    if value is None:
        raise ValueError("cannot format missing numeric value")
    return f"{value:g}"


def _coerce_float(value: object) -> float | None:
    if _is_missing(value):
        return None
    if isinstance(value, int | float):
        return float(value)
    if isinstance(value, str):
        text = value.strip().rstrip("%")
        return float(text)
    return None


def _split_mod_string(value: str) -> list[str]:
    text = value.strip()
    if not text:
        return []
    if text.startswith("{") and text.endswith("}"):
        return [text]
    return [part.strip() for part in re.split(r"\s*,\s*", text) if part.strip()]


def _modification_from_item(item: object) -> SearchedModification:
    if isinstance(item, SearchedModification):
        return item
    if isinstance(item, dict):
        return SearchedModification.model_validate(item)
    return SearchedModification(name=str(item), source=str(item))


def _serialize_modifications(value: list[SearchedModification]) -> str | None:
    if not value:
        return None
    return ",".join(mod.source or mod.name for mod in value)


def _to_scalar(value: object) -> ScalarValue:
    if _is_missing(value):
        return None
    if isinstance(value, str | int | float | bool):
        return value
    return str(value)


def _validate_order(
    minimum: int | float | None,
    maximum: int | float | None,
    label: str,
) -> None:
    if minimum is not None and maximum is not None and minimum > maximum:
        raise ValueError(f"minimum {label} cannot exceed maximum {label}")
