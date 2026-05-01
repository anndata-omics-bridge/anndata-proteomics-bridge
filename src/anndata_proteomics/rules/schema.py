"""Pydantic models for the parsing-rule TOML schema (see ../../../docs/toml_schema.md)."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

InputShape = Literal["long", "wide"]
QuantificationLevel = Literal["ion", "peptidoform", "peptide", "protein"]
EncodingMode = Literal["numeric", "factor"]
DuplicateMode = Literal["error", "aggregate", "keep_first", "keep_all_as_raw_table"]


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid")


class Axis(_Strict):
    obs_keys: list[str] = Field(min_length=1)
    var_keys: list[str] = Field(min_length=1)
    x_layer: str


class Columns(_Strict):
    obs: dict[str, str]
    var: dict[str, str]


class Layer(_Strict):
    name: str
    encoding_mode: EncodingMode = "numeric"
    categories: dict[str, int] | None = None
    source_column: str | None = None
    column_pattern: str | None = None

    @model_validator(mode="after")
    def _factor_requires_categories(self) -> Layer:
        if self.encoding_mode == "factor" and not self.categories:
            raise ValueError(
                f"Layer {self.name!r}: encoding_mode='factor' requires non-empty 'categories'."
            )
        return self


class Duplicates(_Strict):
    mode: DuplicateMode = "error"


class SampleNameCleanup(_Strict):
    pattern: str = ""


class ParseRule(_Strict):
    schema_version: str
    file_version: str
    software_name: str
    software_version: str | None = None
    input_shape: InputShape
    quantification_level: QuantificationLevel
    axis: Axis
    columns: Columns
    layers: list[Layer] = Field(min_length=1)
    duplicates: Duplicates = Field(default_factory=Duplicates)
    sample_name_cleanup: SampleNameCleanup | None = None

    @model_validator(mode="after")
    def _shape_layer_consistency(self) -> ParseRule:
        for layer in self.layers:
            if self.input_shape == "long":
                if layer.source_column is None:
                    raise ValueError(
                        f"Layer {layer.name!r}: long rules require 'source_column'."
                    )
                if layer.column_pattern is not None:
                    raise ValueError(
                        f"Layer {layer.name!r}: 'column_pattern' is only valid for wide rules."
                    )
            else:  # "wide"
                if layer.column_pattern is None:
                    raise ValueError(
                        f"Layer {layer.name!r}: wide rules require 'column_pattern'."
                    )
                if layer.source_column is not None:
                    raise ValueError(
                        f"Layer {layer.name!r}: 'source_column' is only valid for long rules."
                    )
        return self

    @model_validator(mode="after")
    def _x_layer_exists(self) -> ParseRule:
        names = {layer.name for layer in self.layers}
        if self.axis.x_layer not in names:
            raise ValueError(
                f"axis.x_layer={self.axis.x_layer!r} does not match any layer name; "
                f"available: {sorted(names)}."
            )
        return self

    @model_validator(mode="after")
    def _cleanup_only_for_wide(self) -> ParseRule:
        if self.sample_name_cleanup is not None and self.input_shape == "long":
            raise ValueError("sample_name_cleanup is only valid for wide rules.")
        return self
