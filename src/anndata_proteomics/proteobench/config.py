"""Typed readers for ProteoBench module and per-tool TOML settings."""

from __future__ import annotations

import tomllib
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

QuantificationLevel = Literal["ion", "peptidoform"]


class _SettingsModel(BaseModel):
    model_config = ConfigDict(extra="ignore")


class SampleSettings(_SettingsModel):
    """One run in a ProteoBench experiment design."""

    raw_file: str = Field(min_length=1)
    sample_name: str = Field(min_length=1)
    condition: str = Field(min_length=1)


class ExpectedRatio(_SettingsModel):
    """Expected abundance ratio for one species."""

    a_vs_b: float = Field(alias="A_vs_B", gt=0)
    color: str | None = None


class ModuleGeneral(_SettingsModel):
    """Scoring fields from the module ``[general]`` table."""

    min_count_multispec: int = Field(ge=1)
    level: QuantificationLevel
    default_cutoff_min_feature: int = Field(default=1, ge=1)
    max_nr_observed: int = Field(default=6, ge=1)

    @model_validator(mode="after")
    def _validate_cutoffs(self) -> ModuleGeneral:
        if self.default_cutoff_min_feature > self.max_nr_observed:
            raise ValueError("default_cutoff_min_feature must not exceed max_nr_observed")
        return self


class ModuleSettings(_SettingsModel):
    """ProteoBench experiment design and HYE scoring configuration."""

    species_expected_ratio: dict[str, ExpectedRatio]
    species_mapper: dict[str, str]
    general: ModuleGeneral
    samples: list[SampleSettings] = Field(min_length=1)

    @model_validator(mode="after")
    def _validate_design(self) -> ModuleSettings:
        species = list(self.species_mapper.values())
        if len(species) != len(set(species)):
            raise ValueError("species_mapper values must be unique")
        if set(species) != set(self.species_expected_ratio):
            raise ValueError(
                "species_mapper values must equal species_expected_ratio keys; "
                f"mapper={species}, ratios={list(self.species_expected_ratio)}"
            )

        raw_files = [sample.raw_file for sample in self.samples]
        if len(raw_files) != len(set(raw_files)):
            raise ValueError("[[samples]].raw_file values must be unique")
        sample_names = [sample.sample_name for sample in self.samples]
        if len(sample_names) != len(set(sample_names)):
            raise ValueError("[[samples]].sample_name values must be unique")

        conditions = {sample.condition for sample in self.samples}
        if not {"A", "B"} <= conditions:
            raise ValueError(
                "ProteoBench HYE scoring requires conditions 'A' and 'B'; "
                f"found {sorted(conditions)}"
            )
        return self


class ToolGeneral(_SettingsModel):
    """Filtering and run-name settings from a per-tool TOML."""

    contaminant_flag: str | None = None
    decoy_flag: bool | int | str | None = None
    run_name_cleanup: str = ""


class ModificationParserSettings(_SettingsModel):
    """Feature-label behavior needed for ProteoBench compatibility."""

    parse_column: str = Field(min_length=1)
    before_aa: bool
    isalpha: bool = True
    isupper: bool = True
    pattern: str = r"\[([^]]+)\]"
    modification_dict: dict[str, str] = Field(default_factory=dict)


class ToolSettings(_SettingsModel):
    """ProteoBench raw-vendor column interpretation for one tool."""

    mapper: dict[str, str]
    general: ToolGeneral
    run_mapper: dict[str, str] = Field(default_factory=dict)
    condition_mapper: dict[str, str] = Field(default_factory=dict)
    modifications_parser: ModificationParserSettings | None = None

    @model_validator(mode="after")
    def _require_proteins(self) -> ToolSettings:
        protein_sources = [source for source, role in self.mapper.items() if role == "Proteins"]
        if len(protein_sources) != 1:
            raise ValueError(
                "per-tool [mapper] must map exactly one raw column to 'Proteins'; "
                f"found {protein_sources}"
            )
        return self

    def source_for(self, role: str) -> str | None:
        """Return the raw vendor source mapped to a ProteoBench role."""
        return next((source for source, mapped in self.mapper.items() if mapped == role), None)


def load_module_settings(path: str | Path) -> ModuleSettings:
    """Load the scoring subset of a ProteoBench module TOML."""
    return ModuleSettings.model_validate(_load_toml(path))


def load_tool_settings(path: str | Path) -> ToolSettings:
    """Load the column/filtering subset of a ProteoBench per-tool TOML."""
    return ToolSettings.model_validate(_load_toml(path))


def _load_toml(path: str | Path) -> dict[str, Any]:
    settings_path = Path(path)
    with settings_path.open("rb") as handle:
        return tomllib.load(handle)
