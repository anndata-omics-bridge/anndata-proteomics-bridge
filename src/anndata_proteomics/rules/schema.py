"""Pydantic models for the effective parsing-rule JSON schema."""

from __future__ import annotations

import re
from typing import Annotated, Any, Literal

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
    "coalesce",
    "join_nonempty",
    "proforma_sequence",
    "stripped_sequence",
    "proforma_ion",
    "proforma_fragment",
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
    separator: str | None = None

    @model_validator(mode="after")
    def _validate_generic_compute(self) -> ColumnCompute:
        if self.how in {"coalesce", "join_nonempty"} and len(self.from_) < 2:
            raise ValueError(f"how={self.how!r} requires at least two source columns.")
        if self.how == "join_nonempty":
            if self.separator is None or not self.separator:
                raise ValueError("how='join_nonempty' requires a non-empty separator.")
        elif self.separator is not None:
            raise ValueError("separator is valid only for how='join_nonempty'.")
        return self


class ColumnGroup(_Strict):
    select: dict[str, str] = Field(default_factory=dict)
    compute: list[ColumnCompute] = Field(default_factory=list)

    @property
    def names(self) -> list[str]:
        return list(dict.fromkeys([*self.select, *(column.name for column in self.compute)]))


class Columns(_Strict):
    obs: ColumnGroup
    var: ColumnGroup


class ColumnRoles(_Strict):
    """Semantic locations needed by downstream canonical-data consumers."""

    protein_accessions: str = Field(min_length=1)


class Layer(_Strict):
    """A quantitative layer fed by one ``source``.

    ``source`` interpretation is owned by the rule-level ``input_shape``:

    - ``input_shape="long"``: ``source`` is an exact vendor column name.
    - ``input_shape="wide"``: ``source`` is a regex over matrix headers and must
      contain a ``(?P<sample>...)`` named group (enforced on ``ParseRule``).

    Layers are optional by default: captured when the ``source`` is present in the
    input and silently skipped when absent, so a rule can declare the full superset
    of a vendor's volatile columns without rejecting exports that omit some. Set
    ``required = true`` to gate a layer (a file missing it is rejected). The
    ``axis.x_layer`` is always required regardless of this flag.
    """

    name: str
    source: str
    encoding_mode: EncodingMode = "numeric"
    categories: dict[str, int] = Field(default_factory=dict)
    missing_values: list[float] = Field(default_factory=list)
    required: bool = False

    @model_validator(mode="after")
    def _validate_encoding(self) -> Layer:
        if self.encoding_mode == "factor" and not self.categories:
            raise ValueError(
                f"Layer {self.name!r}: encoding_mode='factor' requires non-empty 'categories'."
            )
        if self.encoding_mode == "factor" and self.missing_values:
            raise ValueError(
                f"Layer {self.name!r}: 'missing_values' is only valid for numeric layers."
            )
        return self


class SampleNameCleanup(_Strict):
    pattern: str = ""


class ModificationMapEntry(_Strict):
    """User-facing JSON entry: a vendor token plus the Unimod accession.

    ``name``, ``target``, ``position`` and ``mass_delta`` are NOT carried
    on the entry itself — they are filled at rule-load time from
    ``modifications/unimod_registry.toml``. This keeps the per-tool rule files
    free of duplicated canonical data and guarantees that all tools agree
    on what e.g. ``UNIMOD:35`` means.
    """

    token: str
    accession: str


class TokenRegexModifications(_Strict):
    """Extract vendor modification tokens with a regex and map them to Unimod.

    ``token_regex`` is the only modification parser with a runtime
    implementation (see ``modifications.pipeline``), so it is the only
    accepted ``parser`` value.
    """

    parser: Literal["token_regex"]
    source_column: str
    token_pattern: str
    token_position: TokenPosition = "after_residue"
    case_sensitive: bool = False
    unknown_policy: UnknownPolicy = "preserve"
    output_column: str = "proforma_sequence"
    map: list[ModificationMapEntry] = Field(min_length=1)


# Only ``token_regex`` has a runtime implementation; the alias keeps the public
# name stable for ``ParseRule.modifications`` and importers.
Modifications = TokenRegexModifications


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
    column_roles: ColumnRoles | None = None
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

    def layer_required(self, layer: Layer) -> bool:
        """A layer must be present iff it is the ``x_layer`` or explicitly ``required``."""
        return layer.required or layer.name == self.axis.x_layer

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
    def _column_roles_are_declared(self) -> ParseRule:
        if (
            self.column_roles is not None
            and self.column_roles.protein_accessions not in self.columns.var.names
        ):
            raise ValueError(
                "column_roles.protein_accessions must name a declared var column; "
                f"got {self.column_roles.protein_accessions!r}"
            )
        return self

    @model_validator(mode="after")
    def _fragments_only_for_fragment_level(self) -> ParseRule:
        if self.fragments is not None and self.quantification_level != "fragment":
            raise ValueError("[fragments] is only valid for quantification_level='fragment'.")
        return self

    @model_validator(mode="after")
    def _computed_column_consistency(  # noqa: C901, PLR0912 - schema invariant
        self,
    ) -> ParseRule:
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
            expected_name = _PROFORMA_COMPUTE_NAME.get(column.how)
            if expected_name is not None and column.name != expected_name:
                raise ValueError(
                    f"computed column with how={column.how!r} must be named "
                    f"{expected_name!r}, got {column.name!r}"
                )
            if column.how in {"coalesce", "join_nonempty"}:
                pass
            elif column.how in {"proforma_sequence", "stripped_sequence"}:
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
            else:
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


class PartialAxis(_Strict):
    """Axis fragment that becomes complete after merging a document level with its base."""

    obs_keys: list[str] | None = Field(default=None, min_length=1)
    var_keys: list[str] | None = Field(default=None, min_length=1)
    x_layer: str | None = None
    duplicates: Duplicates | None = None


class PartialColumnGroup(_Strict):
    """Column-group fragment with presence preserved for deterministic merging."""

    select: dict[str, str] | None = None
    compute: list[ColumnCompute] | None = None


class PartialColumns(_Strict):
    """Partial obs/var column declarations used inside a source document."""

    obs: PartialColumnGroup | None = None
    var: PartialColumnGroup | None = None


class RuleFragment(_Strict):
    """Strict partial ParseRule body used by a document base or one level."""

    input_shape: InputShape | None = None
    axis: PartialAxis | None = None
    columns: PartialColumns | None = None
    column_roles: ColumnRoles | None = None
    layers: list[Layer] | None = None
    sample_name_cleanup: SampleNameCleanup | None = None
    modifications: Modifications | None = None
    fragments: Fragments | None = None

    def as_merge_dict(self) -> dict[str, Any]:
        """Return only fields explicitly present in this fragment."""
        return self.model_dump(by_alias=True, exclude_unset=True, mode="json")


class ParseRuleDocument(_Strict):
    """One self-contained software-version document with shared and level rules."""

    schema_version: str
    file_version: str
    software_name: str
    software_version: str
    base: RuleFragment
    levels: dict[QuantificationLevel, RuleFragment] = Field(min_length=1)

    def effective_rule(self, level: QuantificationLevel) -> ParseRule:
        """Merge and validate one level as the effective converter contract."""
        if level not in self.levels:
            raise KeyError(level)
        body = _merge_rule_dicts(
            self.base.as_merge_dict(),
            self.levels[level].as_merge_dict(),
        )
        return ParseRule.model_validate(
            {
                "schema_version": self.schema_version,
                "file_version": self.file_version,
                "software_name": self.software_name,
                "software_version": self.software_version,
                "quantification_level": level,
                **body,
            }
        )

    def effective_rules(self) -> dict[QuantificationLevel, ParseRule]:
        """Validate and return every effective level in stable source order."""
        return {level: self.effective_rule(level) for level in self.levels}


def _is_array_of_objects(value: list[Any]) -> bool:
    """Return whether a JSON array contains rule objects rather than scalar values."""
    return bool(value) and all(isinstance(item, dict) for item in value)


def _merge_rule_dicts(base: dict[str, Any], level: dict[str, Any]) -> dict[str, Any]:
    """Deep-merge a level onto its same-document base using APB's rule semantics."""
    merged = dict(base)
    for key, level_value in level.items():
        base_value = merged.get(key)
        if isinstance(base_value, dict) and isinstance(level_value, dict):
            merged[key] = _merge_rule_dicts(base_value, level_value)
        elif isinstance(base_value, list) and isinstance(level_value, list):
            if _is_array_of_objects(base_value) or _is_array_of_objects(level_value):
                merged[key] = base_value + level_value
            else:
                merged[key] = level_value
        else:
            merged[key] = level_value
    return merged
