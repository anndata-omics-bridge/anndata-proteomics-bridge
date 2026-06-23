"""Pydantic models for the parsing-rule TOML schema (see ../../../docs/toml_schema.md)."""

from __future__ import annotations

import re
from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

InputShape = Literal["long", "wide"]
QuantificationLevel = Literal["ion", "peptidoform", "peptide", "protein", "fragment"]
EncodingMode = Literal["numeric", "factor"]
DuplicateMode = Literal["error", "aggregate", "keep_first", "keep_all_as_raw_table"]
TokenPosition = Literal[
    "before_residue", "after_residue", "n_term", "c_term", "embedded", "unknown"
]
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

_SAMPLE_GROUP = "sample"


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
    """A quantitative layer fed by one ``source``.

    ``source`` interpretation is owned by the rule-level ``input_shape``:

    - ``input_shape="long"``: ``source`` is an exact vendor column name.
    - ``input_shape="wide"``: ``source`` is a regex over matrix headers and must
      contain a ``(?P<sample>...)`` named group (enforced on ``ParseRule``).
    """

    name: str
    source: str
    encoding_mode: EncodingMode = "numeric"
    categories: dict[str, int] = Field(default_factory=dict)

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


class TokenRegexModifications(_Strict):
    """Extract vendor modification tokens with a regex and map them to Unimod.

    The only parser with a runtime implementation today (see
    ``modifications.pipeline``).
    """

    parser: Literal["token_regex"]
    source_column: str
    token_pattern: str
    token_position: TokenPosition = "after_residue"
    case_sensitive: bool = False
    unknown_policy: UnknownPolicy = "preserve"
    output_column: str = "proforma_sequence"
    map: list[ModificationMapEntry] = Field(min_length=1)


class AlreadyProformaModifications(_Strict):
    """``source_column`` already holds a ProForma string; copy it through."""

    parser: Literal["already_proforma"]
    source_column: str
    output_column: str = "proforma_sequence"


class SeparateModColumnModifications(_Strict):
    """Modification tokens live in a column separate from the stripped sequence."""

    parser: Literal["separate_mod_column"]
    source_column: str
    sequence_column: str
    token_position: TokenPosition = "after_residue"
    case_sensitive: bool = False
    unknown_policy: UnknownPolicy = "preserve"
    output_column: str = "proforma_sequence"
    map: list[ModificationMapEntry] = Field(default_factory=list)


Modifications = Annotated[
    TokenRegexModifications | AlreadyProformaModifications | SeparateModColumnModifications,
    Field(discriminator="parser"),
]


class PositionalFragments(_Strict):
    """Older DIA-NN exports with no per-fragment label column.

    Fragment labels are synthesised positionally (``frag_0``, ``frag_1``, …) by
    index within the precursor.
    """

    label_strategy: Literal["positional"]
    value_columns: list[str] = Field(min_length=1)
    delimiter: str = ";"
    label_output: str = "fragment_label"


class ColumnLabeledFragments(_Strict):
    """DIA-NN exports carrying fragment identities in a packed ``label_column``.

    ``label_column`` (e.g. ``Fragment.Info``, tokens like ``b4-unknown^1/327.16``)
    yields ``label_output`` = the token before ``/``.
    """

    label_strategy: Literal["column"]
    value_columns: list[str] = Field(min_length=1)
    label_column: str
    delimiter: str = ";"
    label_output: str = "fragment_label"


Fragments = Annotated[
    PositionalFragments | ColumnLabeledFragments,
    Field(discriminator="label_strategy"),
]
"""Declares packed parallel-list fragment columns to explode before conversion.

DIA-NN-style reports pack per-fragment values as ``delimiter``-joined lists inside each
precursor row (parallel ``Fragment.Quant.*`` lists, aligned by index).
``converters._fragments.explode_fragments`` splits these into one row per fragment before
the normal long-conversion pivot, producing ``label_output`` as a source for a
``proforma_fragment`` computed column. Only valid for ``quantification_level="fragment"``.
"""


class ParseRule(_Strict):
    schema_version: str
    file_version: str
    software_name: str
    software_version: str
    input_shape: InputShape
    quantification_level: QuantificationLevel
    axis: Axis
    columns: Columns
    layers: list[Layer] = Field(min_length=1)
    sample_name_cleanup: SampleNameCleanup | None = None
    modifications: Modifications | None = None
    fragments: Fragments | None = None

    @model_validator(mode="after")
    def _wide_layer_sources_are_sample_regexes(self) -> ParseRule:
        """Wide rules: every ``layer.source`` must be a regex with a ``sample`` group.

        Long rules need no check — ``source`` is a required exact column name.
        """
        if self.input_shape != "wide":
            return self
        for layer in self.layers:
            try:
                pattern = re.compile(layer.source)
            except re.error as exc:
                raise ValueError(
                    f"Layer {layer.name!r}: wide rule 'source' must be a valid regex: {exc}"
                ) from exc
            if _SAMPLE_GROUP not in pattern.groupindex:
                raise ValueError(
                    f"Layer {layer.name!r}: wide rule 'source' must contain a "
                    f"'(?P<{_SAMPLE_GROUP}>...)' named group; got {layer.source!r}."
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
            raise ValueError("[fragments] is only valid for quantification_level='fragment'.")
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
                    raise ValueError(f"how={column.how!r} requires a [modifications] block.")
                if len(column.from_) != 1:
                    raise ValueError(f"how={column.how!r} requires exactly one source column.")
            elif column.how == "proforma_ion":
                # At ion level ProForma_ion is the feature key; at fragment level it is an
                # intermediate used to build ProForma_fragment (not a var key itself).
                if self.quantification_level not in {"ion", "fragment"}:
                    raise ValueError("how='proforma_ion' is valid only for ion or fragment rules.")
                if len(column.from_) != 2:
                    raise ValueError("how='proforma_ion' requires exactly two source columns.")
                if self.quantification_level == "ion" and column.name not in self.axis.var_keys:
                    raise ValueError("computed ProForma ion columns must be used in axis.var_keys.")
            elif column.how == "proforma_fragment":
                if self.quantification_level != "fragment":
                    raise ValueError("how='proforma_fragment' is valid only for fragment rules.")
                if len(column.from_) != 2:
                    raise ValueError("how='proforma_fragment' requires exactly two source columns.")
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
        selected = {source for source in selected_sources if source in derived}
        if selected:
            raise ValueError(
                "APB-derived modification columns must be declared in "
                f"columns.var.compute, not select: {sorted(selected)}"
            )
        return self
