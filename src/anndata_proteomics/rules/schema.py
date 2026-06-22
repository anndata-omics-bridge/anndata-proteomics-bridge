"""Pydantic models for the parsing-rule TOML schema (see ../../../docs/toml_schema.md)."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

InputShape = Literal["long", "wide"]
QuantificationLevel = Literal["ion", "peptidoform", "peptide", "protein", "fragment"]
EncodingMode = Literal["numeric", "factor"]
DuplicateMode = Literal["error", "aggregate", "keep_first", "keep_all_as_raw_table"]
ModificationParser = Literal["token_regex", "already_proforma", "separate_mod_column"]
TokenPosition = Literal["before_residue", "after_residue", "n_term", "c_term", "embedded", "unknown"]
UnknownPolicy = Literal["preserve", "drop", "error"]
ColumnComputeMode = Literal[
    "proforma_sequence", "stripped_sequence", "proforma_ion", "proforma_fragment"
]

_PROFORMA_COMPUTE_NAME = {
    "stripped_sequence": "ProForma_peptide",
    "proforma_sequence": "ProForma_peptidoform",
    "proforma_ion": "ProForma_ion",
    "proforma_fragment": "ProForma_fragment",
}


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid")


class Duplicates(_Strict):
    mode: DuplicateMode = "error"


class Axis(_Strict):
    obs_keys: list[str] = Field(min_length=1)
    var_keys: list[str] = Field(min_length=1)
    x_layer: str
    duplicates: Duplicates = Field(default_factory=Duplicates)


class ColumnCompute(_Strict):
    name: str
    from_: list[str] = Field(alias="from", min_length=1)
    how: ColumnComputeMode


class ColumnGroup(_Strict):
    select: dict[str, str] = Field(default_factory=dict)
    compute: list[ColumnCompute] = Field(default_factory=list)

    @property
    def names(self) -> list[str]:
        return list(self.select) + [column.name for column in self.compute]


class Columns(_Strict):
    obs: ColumnGroup
    var: ColumnGroup


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


class SampleNameCleanup(_Strict):
    pattern: str = ""


class ModificationMapEntry(_Strict):
    """User-facing TOML entry: a vendor token plus the Unimod accession.

    ``name``, ``target``, ``position`` and ``mass_delta`` are NOT carried
    on the entry itself — they are filled at rule-load time from
    ``modifications/unimod_registry.toml``. This keeps the per-tool TOMLs
    free of duplicated canonical data and guarantees that all tools agree
    on what e.g. ``UNIMOD:35`` means.
    """

    token: str
    accession: str


class Modifications(_Strict):
    source_column: str
    parser: ModificationParser = "token_regex"
    token_pattern: str | None = None
    token_position: TokenPosition = "after_residue"
    case_sensitive: bool = False
    unknown_policy: UnknownPolicy = "preserve"
    sequence_column: str | None = None
    output_column: str = "proforma_sequence"
    map: list[ModificationMapEntry] = Field(default_factory=list)

    @model_validator(mode="after")
    def _parser_consistency(self) -> Modifications:
        if self.parser == "token_regex":
            if not self.token_pattern:
                raise ValueError("parser='token_regex' requires 'token_pattern'.")
            if not self.map:
                raise ValueError("parser='token_regex' requires at least one 'map' entry.")
        elif self.parser == "already_proforma":
            if self.token_pattern is not None:
                raise ValueError("parser='already_proforma' must not set 'token_pattern'.")
            if self.map:
                raise ValueError("parser='already_proforma' must not set 'map' entries.")
        # separate_mod_column: source_column suffices; map optional.
        return self


class Fragments(_Strict):
    """Declares packed parallel-list fragment columns to explode before conversion.

    DIA-NN-style reports pack per-fragment values as ``delimiter``-joined lists inside each
    precursor row (parallel ``Fragment.Quant.*`` lists, aligned by index).
    ``converters._fragments.explode_fragments`` splits these into one row per fragment before
    the normal long-conversion pivot, producing ``label_output`` as a source for a
    ``proforma_fragment`` computed column. Only valid for ``quantification_level="fragment"``.

    ``label_column`` is the packed column carrying fragment identities (e.g. ``Fragment.Info``,
    tokens like ``b4-unknown^1/327.16``); ``label_output`` then = the token before ``/``. When
    ``label_column`` is ``None`` (older DIA-NN with no ``Fragment.Info``), labels are
    **positional** — ``frag_0``, ``frag_1``, … by index within the precursor.
    """

    value_columns: list[str] = Field(min_length=1)
    label_column: str | None = None
    delimiter: str = ";"
    label_output: str = "fragment_label"


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
    sample_name_cleanup: SampleNameCleanup | None = None
    modifications: Modifications | None = None
    fragments: Fragments | None = None

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

    @model_validator(mode="after")
    def _axis_keys_are_declared_columns(self) -> ParseRule:
        obs_columns = set(self.columns.obs.names)
        var_columns = set(self.columns.var.names)
        missing_obs = [key for key in self.axis.obs_keys if key not in obs_columns]
        missing_var = [key for key in self.axis.var_keys if key not in var_columns]
        if missing_obs:
            raise ValueError(f"axis.obs_keys must be declared in columns.obs: {missing_obs}")
        if missing_var:
            raise ValueError(f"axis.var_keys must be declared in columns.var: {missing_var}")
        return self

    @model_validator(mode="after")
    def _fragments_only_for_fragment_level(self) -> ParseRule:
        if self.fragments is not None and self.quantification_level != "fragment":
            raise ValueError(
                "[fragments] is only valid for quantification_level='fragment'."
            )
        return self

    @model_validator(mode="after")
    def _computed_column_consistency(self) -> ParseRule:
        available_var_columns = set(self.columns.var.select)
        if self.fragments is not None:
            # explode_fragments injects this column before materialization, so it is a
            # legal `from` source even though it is not a selected vendor column.
            available_var_columns.add(self.fragments.label_output)
        for column in self.columns.var.compute:
            missing_sources = [
                source for source in column.from_ if source not in available_var_columns
            ]
            if missing_sources:
                raise ValueError(
                    f"computed column {column.name!r} references undeclared "
                    f"var column(s): {missing_sources}"
                )
            expected_name = _PROFORMA_COMPUTE_NAME[column.how]
            if column.name != expected_name:
                raise ValueError(
                    f"computed column with how={column.how!r} must be named "
                    f"{expected_name!r}, got {column.name!r}"
                )
            if column.how in {"proforma_sequence", "stripped_sequence"}:
                if self.modifications is None:
                    raise ValueError(
                        f"how={column.how!r} requires a [modifications] block."
                    )
                if len(column.from_) != 1:
                    raise ValueError(f"how={column.how!r} requires exactly one source column.")
            elif column.how == "proforma_ion":
                # At ion level ProForma_ion is the feature key; at fragment level it is an
                # intermediate used to build ProForma_fragment (not a var key itself).
                if self.quantification_level not in {"ion", "fragment"}:
                    raise ValueError(
                        "how='proforma_ion' is valid only for ion or fragment rules."
                    )
                if len(column.from_) != 2:
                    raise ValueError("how='proforma_ion' requires exactly two source columns.")
                if (
                    self.quantification_level == "ion"
                    and column.name not in self.axis.var_keys
                ):
                    raise ValueError(
                        "computed ProForma ion columns must be used in axis.var_keys."
                    )
            elif column.how == "proforma_fragment":
                if self.quantification_level != "fragment":
                    raise ValueError(
                        "how='proforma_fragment' is valid only for fragment rules."
                    )
                if len(column.from_) != 2:
                    raise ValueError(
                        "how='proforma_fragment' requires exactly two source columns."
                    )
                if column.name not in self.axis.var_keys:
                    raise ValueError(
                        "computed ProForma fragment columns must be used in axis.var_keys."
                    )
            available_var_columns.add(column.name)
        if self.columns.obs.compute:
            raise ValueError("computed columns are currently supported only for columns.var.")
        return self

    @model_validator(mode="after")
    def _derived_columns_are_not_selected(self) -> ParseRule:
        if self.modifications is None:
            return self
        derived = {self.modifications.output_column, "stripped_sequence"}
        selected_sources = list(self.columns.obs.select.values()) + list(
            self.columns.var.select.values()
        )
        selected = {
            source
            for source in selected_sources
            if source in derived
        }
        if selected:
            raise ValueError(
                "APB-derived modification columns must be declared in "
                f"columns.var.compute, not select: {sorted(selected)}"
            )
        return self
