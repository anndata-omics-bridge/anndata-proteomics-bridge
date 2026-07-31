"""Pydantic models for the effective parsing-rule JSON schema."""

from __future__ import annotations

import re
from collections.abc import Sequence
from itertools import combinations
from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

from anndata_proteomics.modifications import schema as modification_schema
from anndata_proteomics.params.model import Parameters

ModificationMapEntry = modification_schema.ModificationMapEntry
Modifications = modification_schema.Modifications
SiteListModifications = modification_schema.SiteListModifications
TokenPosition = modification_schema.TokenPosition
TokenRegexModifications = modification_schema.TokenRegexModifications
UnknownPolicy = modification_schema.UnknownPolicy

InputShape = Literal["long", "wide"]
QuantificationLevel = Literal["ion", "peptidoform", "peptide", "protein", "fragment"]
AxisColumnType = Literal["string", "integer", "number", "boolean"]

PEPTIDE_LEVELS: frozenset[QuantificationLevel] = frozenset(
    {"ion", "fragment", "peptidoform", "peptide"}
)
"""Quantification levels whose features carry a peptide sequence."""
EncodingMode = Literal["numeric", "factor"]
DuplicateMode = Literal["error", "aggregate", "keep_first", "keep_all_as_raw_table"]
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


class RuleCompositionError(ValueError):
    """A source document cannot produce a complete effective rule."""

    def __init__(self, message: str, path: tuple[str, ...]) -> None:
        super().__init__(message)
        self.path = path


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
    """Declared axis columns: required ``select``, best-effort ``optional_select``, ``compute``.

    ``select`` sources must be present in the input — they gate rule recognition. Columns a
    vendor emits only for some configurations belong in ``optional_select``: captured when the
    source is present, silently skipped when absent, and never part of recognition. This
    mirrors ``Layer.required = false`` so one vendor document can declare the full superset of
    an export's volatile columns without rejecting the exports that omit some.
    """

    select: dict[str, str] = Field(default_factory=dict)
    optional_select: dict[str, str] = Field(default_factory=dict)
    types: dict[str, AxisColumnType] = Field(default_factory=dict)
    compute: list[ColumnCompute] = Field(default_factory=list)

    def type_for(self, name: str) -> AxisColumnType:
        """Return one selected column's declared type or the string default."""
        return self.types.get(name, "string")

    @property
    def names(self) -> list[str]:
        return list(
            dict.fromkeys(
                [
                    *self.select,
                    *self.optional_select,
                    *(column.name for column in self.compute),
                ]
            )
        )

    @model_validator(mode="after")
    def _select_groups_are_disjoint(self) -> ColumnGroup:
        both = sorted(set(self.select) & set(self.optional_select))
        if both:
            raise ValueError(f"column name(s) declared in both select and optional_select: {both}")
        return self

    @model_validator(mode="after")
    def _types_name_selected_columns(self) -> ColumnGroup:
        selected = set(self.select) | set(self.optional_select)
        unknown = sorted(set(self.types) - selected)
        if unknown:
            raise ValueError(f"types must name selected columns; unknown: {unknown}")
        return self


class Columns(_Strict):
    obs: ColumnGroup
    var: ColumnGroup


class ColumnRoles(_Strict):
    """Semantic locations needed by downstream canonical-data consumers."""

    protein_assignment: str | None = Field(default=None, min_length=1)
    fasta_accessions: str | None = Field(default=None, min_length=1)

    def declared(self) -> dict[str, str]:
        """Return semantic role names mapped to their declared ``var`` columns."""
        return {
            name: column
            for name in ("protein_assignment", "fasta_accessions")
            if (column := getattr(self, name)) is not None
        }


class NoValuePattern(_Strict):
    """A numeric layer contains directly parseable scalar values."""

    mode: Literal["none"] = "none"


class RegexValuePattern(_Strict):
    """Extract one numeric capture group from each structured layer value."""

    mode: Literal["regex"] = "regex"
    pattern: str = Field(min_length=1)

    @model_validator(mode="after")
    def _validate_pattern(self) -> RegexValuePattern:
        try:
            compiled = re.compile(self.pattern)
        except re.error as exc:
            raise ValueError(f"value_pattern is not a valid regex: {exc}") from exc
        if compiled.groups != 1:
            raise ValueError(
                f"value_pattern must have exactly one capture group, found {compiled.groups}."
            )
        return self


type ValuePattern = Annotated[
    NoValuePattern | RegexValuePattern,
    Field(discriminator="mode"),
]

NO_VALUE_PATTERN = NoValuePattern()


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

    A regex ``value_pattern`` handles vendor columns whose cells are structured strings
    rather than bare numbers (PEAKS ``AScore`` is ``site:modification:score``). Its
    ``pattern`` has exactly one capture group, applied per cell before numeric coercion;
    the first match wins and non-matching cells become NaN. ``mode = "none"`` explicitly
    declares directly parseable scalar values.
    """

    name: str
    source: str
    encoding_mode: EncodingMode = "numeric"
    categories: dict[str, int] = Field(default_factory=dict)
    missing_values: list[float] = Field(default_factory=list)
    value_pattern: ValuePattern = NO_VALUE_PATTERN
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
        if self.encoding_mode == "factor" and isinstance(self.value_pattern, RegexValuePattern):
            raise ValueError(
                f"Layer {self.name!r}: 'value_pattern' is only valid for numeric layers."
            )
        return self


class SampleNameCleanup(_Strict):
    pattern: str = Field(min_length=1)


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
    column_roles: ColumnRoles = Field(default_factory=ColumnRoles)
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
    def _axis_keys_are_not_optional(self) -> ParseRule:
        """An axis key must always exist; a feature/observation key cannot be best-effort."""
        optional_obs = [
            key for key in self.axis.obs_keys if key in self.columns.obs.optional_select
        ]
        optional_var = [
            key for key in self.axis.var_keys if key in self.columns.var.optional_select
        ]
        if optional_obs:
            raise ValueError(f"axis.obs_keys must not name optional_select columns: {optional_obs}")
        if optional_var:
            raise ValueError(f"axis.var_keys must not name optional_select columns: {optional_var}")
        return self

    @model_validator(mode="after")
    def _wide_obs_has_no_optional_select(self) -> ParseRule:
        """``converters.wide`` accepts only the ``<sample>`` placeholder on the wide obs axis."""
        if self.input_shape == "wide" and self.columns.obs.optional_select:
            raise ValueError(
                "columns.obs.optional_select is not valid for wide rules; the wide observation "
                "axis accepts only the '<sample>' placeholder"
            )
        return self

    @model_validator(mode="after")
    def _column_roles_are_declared(self) -> ParseRule:
        declared_columns = set(self.columns.var.names)
        for role, column in self.column_roles.declared().items():
            if column not in declared_columns:
                raise ValueError(
                    f"column_roles.{role} must name a declared var column; got {column!r}"
                )
        return self

    @model_validator(mode="after")
    def _fragments_only_for_fragment_level(self) -> ParseRule:
        if self.fragments is not None and self.quantification_level != "fragment":
            raise ValueError("[fragments] is only valid for quantification_level='fragment'.")
        return self

    @model_validator(mode="after")
    def _computed_column_consistency(self) -> ParseRule:
        _validate_computed_columns(self)
        if self.columns.obs.compute:
            raise ValueError("computed columns are currently supported only for columns.var.")
        return self

    @model_validator(mode="after")
    def _derived_columns_are_not_selected(self) -> ParseRule:
        if self.modifications is None:
            return self
        derived = {self.modifications.output_column, "stripped_sequence"}
        selected_sources = [
            *self.columns.obs.select.values(),
            *self.columns.obs.optional_select.values(),
            *self.columns.var.select.values(),
            *self.columns.var.optional_select.values(),
        ]
        selected = {source for source in selected_sources if source in derived}
        if selected:
            raise ValueError(
                "APB-derived modification columns must be declared in "
                f"columns.var.compute, not select: {sorted(selected)}"
            )
        return self


def _validate_computed_columns(rule: ParseRule) -> None:
    available = set(rule.columns.var.select) | set(rule.columns.var.optional_select)
    if rule.fragments is not None:
        # Fragment expansion injects this source before computed columns materialize.
        available.add(rule.fragments.label_output)
    for column in rule.columns.var.compute:
        _validate_computed_column(rule, column, available)
        available.add(column.name)


def _validate_computed_column(
    rule: ParseRule,
    column: ColumnCompute,
    available: set[str],
) -> None:
    missing_sources = [source for source in column.from_ if source not in available]
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
        return
    if column.how in {"proforma_sequence", "stripped_sequence"}:
        _validate_sequence_compute(rule, column)
        return
    if column.how == "proforma_ion":
        _validate_ion_compute(rule, column)
        return
    _validate_fragment_compute(rule, column)


def _validate_sequence_compute(rule: ParseRule, column: ColumnCompute) -> None:
    if rule.modifications is None:
        raise ValueError(f"how={column.how!r} requires a [modifications] block.")
    if len(column.from_) != 1:
        raise ValueError(f"how={column.how!r} requires exactly one source column.")


def _validate_ion_compute(rule: ParseRule, column: ColumnCompute) -> None:
    # At fragment level ProForma_ion is an intermediate for ProForma_fragment.
    if rule.quantification_level not in {"ion", "fragment"}:
        raise ValueError("how='proforma_ion' is valid only for ion or fragment rules.")
    if len(column.from_) != 2:
        raise ValueError("how='proforma_ion' requires exactly two source columns.")
    charge_column = column.from_[1]
    if rule.columns.var.type_for(charge_column) != "integer":
        raise ValueError(
            "how='proforma_ion' requires its charge source to declare type='integer'; "
            f"got {charge_column!r}"
        )
    if rule.quantification_level == "ion" and column.name not in rule.axis.var_keys:
        raise ValueError("computed ProForma ion columns must be used in axis.var_keys.")


def _validate_fragment_compute(rule: ParseRule, column: ColumnCompute) -> None:
    if rule.quantification_level != "fragment":
        raise ValueError("how='proforma_fragment' is valid only for fragment rules.")
    if len(column.from_) != 2:
        raise ValueError("how='proforma_fragment' requires exactly two source columns.")
    if column.name not in rule.axis.var_keys:
        raise ValueError("computed ProForma fragment columns must be used in axis.var_keys.")


class PartialAxis(_Strict):
    """Axis fragment that becomes complete after merging a document level with its base."""

    obs_keys: list[str] | None = Field(default=None, min_length=1)
    var_keys: list[str] | None = Field(default=None, min_length=1)
    x_layer: str | None = None
    duplicates: Duplicates | None = None


class PartialColumnGroup(_Strict):
    """Column-group fragment with presence preserved for deterministic merging."""

    select: dict[str, str] | None = None
    optional_select: dict[str, str] | None = None
    types: dict[str, AxisColumnType] | None = None
    compute: list[ColumnCompute] | None = None


class PartialColumns(_Strict):
    """Partial obs/var column declarations used inside a source document."""

    obs: PartialColumnGroup | None = None
    var: PartialColumnGroup | None = None


class RuleFragment(_Strict):
    """Strict partial ParseRule body used by a document base."""

    input_shape: InputShape | None = None
    axis: PartialAxis | None = None
    columns: PartialColumns | None = None
    column_roles: ColumnRoles | None = None
    layers: list[Layer] | None = None
    sample_name_cleanup: SampleNameCleanup | None = None
    modifications: Modifications | None = None
    fragments: Fragments | None = None


class SearchParameterCondition(Parameters):
    """Typed equality conditions over explicitly declared search-parameter fields."""

    def matches(self, search_parameters: Parameters) -> bool:
        """Return whether every explicitly declared condition value matches."""
        fields = self.model_fields_set
        expected = self.model_dump(include=fields, mode="python")
        observed = search_parameters.model_dump(include=fields, mode="python")
        return expected == observed

    def agrees_with(self, other: SearchParameterCondition) -> bool:
        """Return whether overlapping equality conditions agree."""
        shared_fields = self.model_fields_set & other.model_fields_set
        own_values = self.model_dump(include=shared_fields, mode="python")
        other_values = other.model_dump(include=shared_fields, mode="python")
        return own_values == other_values


class SearchParameterOverride(_Strict):
    """One search-parameter condition and its level-local axis override."""

    when_search_parameters: SearchParameterCondition = Field(json_schema_extra={"minProperties": 1})
    axis: PartialAxis

    @model_validator(mode="after")
    def _fragments_are_not_empty(self) -> SearchParameterOverride:
        if not self.when_search_parameters.model_fields_set:
            raise ValueError("search-parameter condition must declare at least one field")
        if not self.axis.model_fields_set:
            raise ValueError("search-parameter override axis must declare at least one field")
        return self

    def matches(self, search_parameters: Parameters) -> bool:
        """Return whether all normalized parameter conditions match."""
        return self.when_search_parameters.matches(search_parameters)


class LevelRuleFragment(RuleFragment):
    """Strict level body with non-recursive search-parameter axis overrides.

    ``requires_search_parameters`` gates the level's *availability* rather than its content:
    a tool whose quantification settings decide the level of its own output table (Sage's
    ``combine_charge_states`` collapses charge states, making one and the same ``lfq.tsv``
    schema ion- or peptidoform-level) cannot be resolved by version or headers, because
    neither differs. A level declaring conditions here is offered only when the parsed
    search parameters satisfy every one of them. Levels with no conditions are always
    available, so this changes nothing for documents that do not use it.
    """

    search_parameter_overrides: list[SearchParameterOverride] = Field(default_factory=list)
    requires_search_parameters: SearchParameterCondition = Field(
        default_factory=SearchParameterCondition
    )

    def is_unconditionally_available(self) -> bool:
        """Return whether this level requires no parsed search-parameter evidence."""
        return not self.requires_search_parameters.model_fields_set

    def is_available_for(self, search_parameters: Parameters) -> bool:
        """Return whether this level is offered for concrete search parameters."""
        return self.is_unconditionally_available() or self.requires_search_parameters.matches(
            search_parameters
        )


class ParseRuleDocument(_Strict):
    """One self-contained software-version document with shared and level rules."""

    schema_version: str
    file_version: str
    software_name: str
    software_version: str
    base: RuleFragment
    levels: dict[QuantificationLevel, LevelRuleFragment] = Field(min_length=1)

    def effective_rule(
        self,
        level: QuantificationLevel,
    ) -> ParseRule:
        """Compose one level without search-parameter overrides."""
        if level not in self.levels:
            raise KeyError(level)
        return self._materialize_rule(level, [])

    def parameterized_effective_rule(
        self,
        level: QuantificationLevel,
        search_parameters: Parameters,
    ) -> ParseRule:
        """Compose one level with every matching search-parameter override."""
        if level not in self.levels:
            raise KeyError(level)
        level_fragment = self.levels[level]
        matching_overrides = [
            override
            for override in level_fragment.search_parameter_overrides
            if override.matches(search_parameters)
        ]
        return self._materialize_rule(level, matching_overrides)

    def level_is_unconditionally_available(
        self,
        level: QuantificationLevel,
    ) -> bool:
        """Return whether ``level`` exists and requires no parameter evidence."""
        if level not in self.levels:
            return False
        return self.levels[level].is_unconditionally_available()

    def level_is_available_for(
        self,
        level: QuantificationLevel,
        search_parameters: Parameters,
    ) -> bool:
        """Return whether ``level`` exists and its parameter gate is satisfied."""
        if level not in self.levels:
            return False
        return self.levels[level].is_available_for(search_parameters)

    def unconditionally_available_levels(self) -> list[QuantificationLevel]:
        """Return ungated levels in stable source order."""
        return [
            level
            for level, fragment in self.levels.items()
            if fragment.is_unconditionally_available()
        ]

    def available_levels_for(
        self,
        search_parameters: Parameters,
    ) -> list[QuantificationLevel]:
        """Return levels whose parameter gate is satisfied, in source order."""
        return [
            level
            for level, fragment in self.levels.items()
            if fragment.is_available_for(search_parameters)
        ]

    def effective_rules(self) -> dict[QuantificationLevel, ParseRule]:
        """Return every declared level without search-parameter overrides."""
        return {level: self.effective_rule(level) for level in self.levels}

    def parameterized_effective_rules(
        self,
        search_parameters: Parameters,
    ) -> dict[QuantificationLevel, ParseRule]:
        """Return every declared level with matching parameter overrides applied.

        Returns every declared level regardless of its availability gate: this is the
        document's full contract, used by validation and rule listing. Resolution filters
        by ``available_levels_for``.
        """
        return {
            level: self.parameterized_effective_rule(level, search_parameters)
            for level in self.levels
        }

    def validate_effective_rule_variants(self, level: QuantificationLevel) -> None:
        """Validate the default and every condition-compatible override combination."""
        if level not in self.levels:
            raise KeyError(level)
        self._materialize_rule(level, [])
        overrides = self.levels[level].search_parameter_overrides
        for size in range(1, len(overrides) + 1):
            for selected in combinations(overrides, size):
                if _search_parameter_conditions_are_compatible(selected):
                    self._materialize_rule(level, selected)

    def _materialize_rule(
        self,
        level: QuantificationLevel,
        overrides: Sequence[SearchParameterOverride],
    ) -> ParseRule:
        """Compose typed fragments into one effective converter contract."""
        level_fragment = self.levels[level]
        sample_name_cleanup = self.base.sample_name_cleanup
        if "sample_name_cleanup" in level_fragment.model_fields_set:
            level_cleanup = level_fragment.sample_name_cleanup
            base_cleanup = self.base.sample_name_cleanup
            sample_name_cleanup = level_cleanup
            if level_cleanup is not None and base_cleanup is not None:
                sample_name_cleanup = _compose_sample_name_cleanup(
                    base_cleanup,
                    level_cleanup,
                )

        modifications = self.base.modifications
        if "modifications" in level_fragment.model_fields_set:
            level_modifications = level_fragment.modifications
            base_modifications = self.base.modifications
            modifications = level_modifications
            if level_modifications is not None and base_modifications is not None:
                modifications = _compose_modifications(
                    base_modifications,
                    level_modifications,
                )

        fragments = self.base.fragments
        if "fragments" in level_fragment.model_fields_set:
            level_fragments = level_fragment.fragments
            base_fragments = self.base.fragments
            fragments = level_fragments
            if level_fragments is not None and base_fragments is not None:
                fragments = _compose_fragments(base_fragments, level_fragments)

        return ParseRule(
            schema_version=self.schema_version,
            file_version=self.file_version,
            software_name=self.software_name,
            software_version=self.software_version,
            quantification_level=level,
            input_shape=_compose_input_shape(self.base, level_fragment),
            axis=_compose_axis(self.base, level_fragment, overrides),
            columns=_compose_columns(self.base, level_fragment),
            column_roles=_compose_column_roles(self.base, level_fragment),
            layers=_compose_layers(self.base, level_fragment),
            sample_name_cleanup=sample_name_cleanup,
            modifications=modifications,
            fragments=fragments,
        )


def _search_parameter_conditions_are_compatible(
    overrides: tuple[SearchParameterOverride, ...],
) -> bool:
    """Return whether selected equality conditions can form valid parameters."""
    conditions = [override.when_search_parameters for override in overrides]
    equalities_agree = all(left.agrees_with(right) for left, right in combinations(conditions, 2))
    return (
        equalities_agree
        and _charge_range_is_compatible(conditions)
        and _peptide_length_range_is_compatible(conditions)
        and _precursor_mz_range_is_compatible(conditions)
        and _fragment_mz_range_is_compatible(conditions)
    )


def _compose_input_shape(base: RuleFragment, level: LevelRuleFragment) -> InputShape:
    """Return the declared input shape after the level overrides the base."""
    value = level.input_shape if "input_shape" in level.model_fields_set else base.input_shape
    if value is None:
        raise RuleCompositionError("effective rule requires input_shape", ("input_shape",))
    return value


def _compose_axis(
    base: RuleFragment,
    level: LevelRuleFragment,
    overrides: Sequence[SearchParameterOverride],
) -> Axis:
    """Compose partial axis declarations in source order."""
    fragments: list[PartialAxis] = []
    if base.axis is not None:
        fragments.append(base.axis)
    if "axis" in level.model_fields_set:
        if level.axis is None:
            fragments.clear()
        else:
            fragments.append(level.axis)
    fragments.extend(parameter_override.axis for parameter_override in overrides)
    return Axis(
        obs_keys=_compose_axis_obs_keys(fragments),
        var_keys=_compose_axis_var_keys(fragments),
        x_layer=_compose_axis_x_layer(fragments),
        duplicates=_compose_axis_duplicates(fragments),
    )


def _compose_axis_obs_keys(fragments: Sequence[PartialAxis]) -> list[str]:
    """Return the last declared observation keys, or raise a schema error."""
    for fragment in reversed(fragments):
        if "obs_keys" in fragment.model_fields_set:
            if fragment.obs_keys is None:
                raise RuleCompositionError(
                    "effective rule requires axis.obs_keys", ("axis", "obs_keys")
                )
            return fragment.obs_keys
    raise RuleCompositionError("effective rule requires axis.obs_keys", ("axis", "obs_keys"))


def _compose_axis_var_keys(fragments: Sequence[PartialAxis]) -> list[str]:
    """Return the last declared feature keys, or raise a schema error."""
    for fragment in reversed(fragments):
        if "var_keys" in fragment.model_fields_set:
            if fragment.var_keys is None:
                raise RuleCompositionError(
                    "effective rule requires axis.var_keys", ("axis", "var_keys")
                )
            return fragment.var_keys
    raise RuleCompositionError("effective rule requires axis.var_keys", ("axis", "var_keys"))


def _compose_axis_x_layer(fragments: Sequence[PartialAxis]) -> str:
    """Return the last declared primary layer, or raise a schema error."""
    for fragment in reversed(fragments):
        if "x_layer" in fragment.model_fields_set:
            if fragment.x_layer is None:
                raise RuleCompositionError(
                    "effective rule requires axis.x_layer", ("axis", "x_layer")
                )
            return fragment.x_layer
    raise RuleCompositionError("effective rule requires axis.x_layer", ("axis", "x_layer"))


def _compose_axis_duplicates(fragments: Sequence[PartialAxis]) -> Duplicates:
    """Return the last declared duplicate policy or its schema default."""
    for fragment in reversed(fragments):
        if "duplicates" in fragment.model_fields_set:
            if fragment.duplicates is None:
                raise RuleCompositionError("axis.duplicates cannot be null", ("axis", "duplicates"))
            return fragment.duplicates
    return Duplicates()


def _compose_columns(base: RuleFragment, level: LevelRuleFragment) -> Columns:
    """Compose partial observation and feature column declarations."""
    fragments: list[PartialColumns] = []
    if base.columns is not None:
        fragments.append(base.columns)
    if "columns" in level.model_fields_set:
        if level.columns is None:
            fragments.clear()
        else:
            fragments.append(level.columns)
    return Columns(
        obs=_compose_axis_column_group(fragments, "obs"),
        var=_compose_axis_column_group(fragments, "var"),
    )


def _compose_axis_column_group(
    columns: Sequence[PartialColumns],
    group: Literal["obs", "var"],
) -> ColumnGroup:
    """Compose one required axis column group from typed fragments."""
    fragments: list[PartialColumnGroup] = []
    is_null = False
    for fragment in columns:
        if group not in fragment.model_fields_set:
            continue
        value = fragment.obs if group == "obs" else fragment.var
        if value is None:
            fragments.clear()
            is_null = True
        else:
            if is_null:
                fragments.clear()
            fragments.append(value)
            is_null = False
    if is_null or not fragments:
        raise RuleCompositionError(f"effective rule requires columns.{group}", ("columns", group))
    return _compose_column_group(fragments, group)


def _compose_column_group(
    fragments: Sequence[PartialColumnGroup],
    group: Literal["obs", "var"],
) -> ColumnGroup:
    """Merge typed mapping fields and append typed computed-column declarations."""
    return ColumnGroup(
        select=_compose_selected_columns(fragments, group),
        optional_select=_compose_optional_columns(fragments, group),
        types=_compose_column_types(fragments, group),
        compute=_compose_computed_columns(fragments, group),
    )


def _compose_selected_columns(
    fragments: Sequence[PartialColumnGroup],
    group: Literal["obs", "var"],
) -> dict[str, str]:
    """Compose required source-column mappings."""
    selected: dict[str, str] = {}
    is_null = False
    for fragment in fragments:
        if "select" not in fragment.model_fields_set:
            continue
        if fragment.select is None:
            selected.clear()
            is_null = True
        else:
            if is_null:
                selected.clear()
            selected.update(fragment.select)
            is_null = False
    if is_null:
        raise RuleCompositionError(
            f"columns.{group}.select cannot be null", ("columns", group, "select")
        )
    return selected


def _compose_optional_columns(
    fragments: Sequence[PartialColumnGroup],
    group: Literal["obs", "var"],
) -> dict[str, str]:
    """Compose best-effort source-column mappings."""
    selected: dict[str, str] = {}
    is_null = False
    for fragment in fragments:
        if "optional_select" not in fragment.model_fields_set:
            continue
        if fragment.optional_select is None:
            selected.clear()
            is_null = True
        else:
            if is_null:
                selected.clear()
            selected.update(fragment.optional_select)
            is_null = False
    if is_null:
        raise RuleCompositionError(
            f"columns.{group}.optional_select cannot be null",
            ("columns", group, "optional_select"),
        )
    return selected


def _compose_column_types(
    fragments: Sequence[PartialColumnGroup],
    group: Literal["obs", "var"],
) -> dict[str, AxisColumnType]:
    """Compose logical selected-column types; omitted entries default to string."""
    declared: dict[str, AxisColumnType] = {}
    is_null = False
    for fragment in fragments:
        if "types" not in fragment.model_fields_set:
            continue
        if fragment.types is None:
            declared.clear()
            is_null = True
        else:
            if is_null:
                declared.clear()
            declared.update(fragment.types)
            is_null = False
    if is_null:
        raise RuleCompositionError(
            f"columns.{group}.types cannot be null",
            ("columns", group, "types"),
        )
    return declared


def _compose_computed_columns(
    fragments: Sequence[PartialColumnGroup],
    group: Literal["obs", "var"],
) -> list[ColumnCompute]:
    """Append computed-column declarations in document order."""
    computed: list[ColumnCompute] = []
    is_null = False
    for fragment in fragments:
        if "compute" not in fragment.model_fields_set:
            continue
        if fragment.compute is None:
            computed.clear()
            is_null = True
        else:
            if is_null:
                computed.clear()
            computed.extend(fragment.compute)
            is_null = False
    if is_null:
        raise RuleCompositionError(
            f"columns.{group}.compute cannot be null", ("columns", group, "compute")
        )
    return computed


def _compose_column_roles(base: RuleFragment, level: LevelRuleFragment) -> ColumnRoles:
    """Compose semantic column-role declarations field by field."""
    base_roles = base.column_roles
    roles = base_roles
    if "column_roles" in level.model_fields_set:
        if level.column_roles is None:
            raise RuleCompositionError("column_roles cannot be null", ("column_roles",))
        roles = level.column_roles
    elif "column_roles" in base.model_fields_set and base_roles is None:
        raise RuleCompositionError("column_roles cannot be null", ("column_roles",))
    if roles is None:
        return ColumnRoles()
    if base_roles is None or roles is base_roles:
        return roles

    protein_assignment = base_roles.protein_assignment
    fasta_accessions = base_roles.fasta_accessions
    if "protein_assignment" in roles.model_fields_set:
        protein_assignment = roles.protein_assignment
    if "fasta_accessions" in roles.model_fields_set:
        fasta_accessions = roles.fasta_accessions
    return ColumnRoles(
        protein_assignment=protein_assignment,
        fasta_accessions=fasta_accessions,
    )


def _compose_layers(base: RuleFragment, level: LevelRuleFragment) -> list[Layer]:
    """Append typed layer declarations from the base and selected level."""
    if "layers" not in level.model_fields_set:
        selected = base.layers
    elif level.layers is None:
        selected = None
    elif base.layers is None:
        selected = level.layers
    else:
        selected = [*base.layers, *level.layers]
    if selected is None:
        raise RuleCompositionError("effective rule requires layers", ("layers",))
    if not selected:
        raise RuleCompositionError("effective rule requires at least one layer", ("layers",))
    return selected


def _compose_sample_name_cleanup(
    base: SampleNameCleanup,
    level: SampleNameCleanup,
) -> SampleNameCleanup:
    """Compose two present wide-sample cleanup declarations."""
    if "pattern" in level.model_fields_set:
        return level
    return base


def _compose_modifications(
    base: Modifications,
    level: Modifications,
) -> Modifications:
    """Compose two present modification-parser declarations."""
    if isinstance(base, TokenRegexModifications) and isinstance(level, TokenRegexModifications):
        return _compose_token_regex_modifications(base, level)
    if isinstance(base, SiteListModifications) and isinstance(level, SiteListModifications):
        return _compose_site_list_modifications(base, level)
    raise ValueError("a level modification parser cannot change the parser declared by its base")


def _compose_fragments(base: Fragments, level: Fragments) -> Fragments:
    """Compose two present packed-fragment declarations."""
    if isinstance(base, PositionalFragments) and isinstance(level, PositionalFragments):
        return _compose_positional_fragments(base, level)
    if isinstance(base, ColumnLabeledFragments) and isinstance(level, ColumnLabeledFragments):
        return _compose_column_labeled_fragments(base, level)
    raise ValueError(
        "a level fragment declaration cannot change the label strategy declared by its base"
    )


def _compose_token_regex_modifications(
    base: TokenRegexModifications,
    level: TokenRegexModifications,
) -> TokenRegexModifications:
    """Compose two token-regex parser declarations without untyped dictionaries."""
    return TokenRegexModifications(
        parser="token_regex",
        source_column=level.source_column,
        token_pattern=level.token_pattern,
        token_position=(
            level.token_position
            if "token_position" in level.model_fields_set
            else base.token_position
        ),
        case_sensitive=(
            level.case_sensitive
            if "case_sensitive" in level.model_fields_set
            else base.case_sensitive
        ),
        unknown_policy=(
            level.unknown_policy
            if "unknown_policy" in level.model_fields_set
            else base.unknown_policy
        ),
        output_column=(
            level.output_column if "output_column" in level.model_fields_set else base.output_column
        ),
        map=[*base.map, *level.map],
    )


def _compose_site_list_modifications(
    base: SiteListModifications,
    level: SiteListModifications,
) -> SiteListModifications:
    """Compose two site-list parser declarations without untyped dictionaries."""
    return SiteListModifications(
        parser="site_list",
        sequence_column=level.sequence_column,
        modification_column=level.modification_column,
        site_column=level.site_column,
        delimiter=(level.delimiter if "delimiter" in level.model_fields_set else base.delimiter),
        site_base=level.site_base if "site_base" in level.model_fields_set else base.site_base,
        case_sensitive=(
            level.case_sensitive
            if "case_sensitive" in level.model_fields_set
            else base.case_sensitive
        ),
        unknown_policy=(
            level.unknown_policy
            if "unknown_policy" in level.model_fields_set
            else base.unknown_policy
        ),
        output_column=(
            level.output_column if "output_column" in level.model_fields_set else base.output_column
        ),
        map=[*base.map, *level.map],
    )


def _compose_positional_fragments(
    base: PositionalFragments,
    level: PositionalFragments,
) -> PositionalFragments:
    """Compose two positional-fragment declarations field by field."""
    return PositionalFragments(
        label_strategy="positional",
        value_columns=level.value_columns,
        delimiter=(level.delimiter if "delimiter" in level.model_fields_set else base.delimiter),
        label_output=(
            level.label_output if "label_output" in level.model_fields_set else base.label_output
        ),
    )


def _compose_column_labeled_fragments(
    base: ColumnLabeledFragments,
    level: ColumnLabeledFragments,
) -> ColumnLabeledFragments:
    """Compose two column-labeled fragment declarations field by field."""
    return ColumnLabeledFragments(
        label_strategy="column",
        value_columns=level.value_columns,
        label_column=level.label_column,
        delimiter=(level.delimiter if "delimiter" in level.model_fields_set else base.delimiter),
        label_output=(
            level.label_output if "label_output" in level.model_fields_set else base.label_output
        ),
    )


def _charge_range_is_compatible(conditions: Sequence[SearchParameterCondition]) -> bool:
    minimum = None
    maximum = None
    for condition in conditions:
        if "min_precursor_charge" in condition.model_fields_set:
            minimum = condition.min_precursor_charge
        if "max_precursor_charge" in condition.model_fields_set:
            maximum = condition.max_precursor_charge
    return minimum is None or maximum is None or minimum <= maximum


def _peptide_length_range_is_compatible(
    conditions: Sequence[SearchParameterCondition],
) -> bool:
    minimum = None
    maximum = None
    for condition in conditions:
        if "min_peptide_length" in condition.model_fields_set:
            minimum = condition.min_peptide_length
        if "max_peptide_length" in condition.model_fields_set:
            maximum = condition.max_peptide_length
    return minimum is None or maximum is None or minimum <= maximum


def _precursor_mz_range_is_compatible(conditions: Sequence[SearchParameterCondition]) -> bool:
    minimum = None
    maximum = None
    for condition in conditions:
        if "min_precursor_mz" in condition.model_fields_set:
            minimum = condition.min_precursor_mz
        if "max_precursor_mz" in condition.model_fields_set:
            maximum = condition.max_precursor_mz
    return minimum is None or maximum is None or minimum <= maximum


def _fragment_mz_range_is_compatible(conditions: Sequence[SearchParameterCondition]) -> bool:
    minimum = None
    maximum = None
    for condition in conditions:
        if "min_fragment_mz" in condition.model_fields_set:
            minimum = condition.min_fragment_mz
        if "max_fragment_mz" in condition.model_fields_set:
            maximum = condition.max_fragment_mz
    return minimum is None or maximum is None or minimum <= maximum
