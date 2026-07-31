"""Architectural guards for typed, backend-independent computation."""

from __future__ import annotations

import ast
from collections.abc import Callable
from pathlib import Path
from typing import override

import pytest

PACKAGE = Path(__file__).resolve().parents[1] / "src" / "anndata_proteomics"

# Everything under src/ is a guarded pure computation module unless it is a composition
# root, a storage adapter, or one of the reviewed exemptions below. A new package is
# therefore guarded on creation rather than when somebody remembers to list it.
UNGUARDED_DIRECTORIES = frozenset({"adapters", "scripts"})
UNGUARDED_FILES = frozenset(
    {
        # External input boundaries: vendor tables and vendor parameter files arrive
        # untyped and are narrowed here, so nullable parsed fields are the contract.
        "params/model.py",
        "params/parsers/_common.py",
        "params/parsers/alphadia.py",
        "params/parsers/alphapept.py",
        "params/parsers/diann.py",
        "params/parsers/fragpipe.py",
        "params/parsers/maxquant.py",
        "params/parsers/metamorpheus.py",
        "params/parsers/msaid.py",
        "params/parsers/peaks.py",
        "params/parsers/sage.py",
        "params/parsers/spectronaut.py",
        "params/parsers/wombat.py",
        "params/registry.py",
        "readers/dispatch.py",
        "readers/tabular.py",
        # Validated schema contracts: Pydantic models whose optional fields describe a
        # document that genuinely omits them.
        "modifications/model.py",
        "modifications/schema.py",
        "modifications/sdrf.py",
        "rules/loader.py",
        "rules/registry.py",
        "rules/schema.py",
        # Test-data catalogue lookup shipped with the package.
        "test_data.py",
    }
)


def _guarded_files() -> tuple[Path, ...]:
    """Return every module the pure-computation rules apply to."""
    return tuple(
        path
        for path in sorted(PACKAGE.rglob("*.py"))
        if path.name != "__init__.py"
        and path.relative_to(PACKAGE).parts[0] not in UNGUARDED_DIRECTORIES
        and path.relative_to(PACKAGE).as_posix() not in UNGUARDED_FILES
    )


PURE_FILES = _guarded_files()
WORKFLOW_FILES = tuple(sorted((PACKAGE / "workflows").glob("*.py")))
ADAPTER_FILES = tuple(sorted((PACKAGE / "adapters").rglob("*.py")))

CONTAINER_SLOTS = {
    "X",
    "layers",
    "mod",
    "obs",
    "obsm",
    "uns",
    "var",
    "varm",
    "varp",
}

SCIENTIFIC_CONTAINER_SLOTS = CONTAINER_SLOTS - {"mod", "uns"}
DATAFRAME_MANIPULATION_METHODS = {
    "apply",
    "assign",
    "astype",
    "copy",
    "drop",
    "drop_duplicates",
    "fillna",
    "groupby",
    "join",
    "map",
    "merge",
    "pivot",
    "pivot_table",
    "reindex",
    "rename",
    "replace",
    "reset_index",
    "set_index",
    "sort_index",
    "sort_values",
    "stack",
    "to_numpy",
    "transform",
    "unstack",
}

MIXED_FASTA_ADAPTER_ENTRYPOINTS = {
    "annotate_protein_anndata_from_fasta",
    "annotate_protein_mudata_from_fasta",
    "validate_peptide_anndata_against_fasta",
    "validate_peptide_mudata_against_fasta",
    "validate_peptide_mudata_modality_against_fasta",
}

# One table for the whole package: each adapter and the calculation entry points it must
# never call itself. Folded in from the former tests/test_fasta_architecture.py.
FORBIDDEN_ADAPTER_ENTRY_POINTS: dict[str, frozenset[str]] = {
    "adapters/anndata/fasta.py": frozenset(
        {
            "annotate_proteins_from_fasta",
            "build_feature_mapping",
            "build_mulink_feature_mapping",
            "calculate_feature_mapping_update",
            "combined_peptide_protein_matches",
            "custom_cleavage",
            "match_peptide_collections_to_fasta",
            "match_peptides_to_fasta",
            "merge_owned_feature_mapping",
            "peptide_feature_nodes",
            "peptide_protein_matches",
            "protein_feature_nodes",
            "protein_group_accessions",
            "resolve_cleavage_name",
            "resolve_protein_annotation_input",
            "unavailable_reported_protein_validation",
            "validate_peptide_levels",
            "validate_reported_proteins",
        }
    ),
    "adapters/anndata/proteobench.py": frozenset(
        {
            "align_runs",
            "build_scores",
            "compute_intermediate",
            "score_level",
            "score_levels",
        }
    ),
    "adapters/anndata/annotation.py": frozenset(
        {
            "annotate_observations",
            "annotation_diagnostics",
            "match_sample_annotation",
            "run_sample_annotation",
        }
    ),
    "adapters/anndata/conversion.py": frozenset({"convert_long", "convert_table", "convert_wide"}),
}

NO_EXCEPTION_CONTROL_FLOW = (
    (PACKAGE / "readers" / "tabular.py", "detect_text_delimiter"),
    (PACKAGE / "params" / "model.py", "_coerce_scan_window"),
    (PACKAGE / "params" / "parsers" / "alphadia.py", "_ppm"),
    (PACKAGE / "params" / "parsers" / "diann.py", "_parse_diann_version"),
    (PACKAGE / "params" / "parsers" / "maxquant.py", "_min_peptide_length"),
    (PACKAGE / "params" / "parsers" / "_common.py", "read_text"),
)

type SignatureSlot = tuple[str, str, str, str]
type AnnotationPredicate = Callable[[ast.expr], bool]

# Exact, reviewed external-schema boundaries. Removing a nullable slot or raw-object
# normalizer must also remove it here; adding one requires an explicit architecture review.
AUDITED_OPTIONAL_SIGNATURE_SLOTS: frozenset[SignatureSlot] = frozenset(
    {
        # Pydantic field validators and model normalization.
        (
            "modifications/model.py",
            "ModificationOccurrence._non_negative_index",
            "parameter",
            "value",
        ),
        (
            "modifications/model.py",
            "ModificationOccurrence._non_negative_index",
            "return",
            "return",
        ),
        ("modifications/model.py", "ModificationOccurrence._valid_accession", "parameter", "value"),
        ("modifications/model.py", "ModificationOccurrence._valid_accession", "return", "return"),
        ("modifications/model.py", "SearchedModification._valid_accession", "parameter", "value"),
        ("modifications/model.py", "SearchedModification._valid_accession", "return", "return"),
        ("params/model.py", "MassTolerance.parse", "return", "return"),
        ("params/model.py", "Probability.parse", "return", "return"),
        ("params/model.py", "_coerce_float", "return", "return"),
        ("params/model.py", "_serialize_modifications", "return", "return"),
        ("params/model.py", "_validate_order", "parameter", "maximum"),
        ("params/model.py", "_validate_order", "parameter", "minimum"),
        # Vendor formats genuinely omit these parsed fields.
        ("params/parsers/alphadia.py", "_nested_range", "return", "return"),
        ("params/parsers/diann.py", "_arguments", "return", "return"),
        ("params/parsers/diann.py", "_extract_cfg_float", "return", "return"),
        ("params/parsers/diann.py", "_extract_cfg_int", "return", "return"),
        ("params/parsers/diann.py", "_extract_cfg_text", "return", "return"),
        ("params/parsers/diann.py", "_extract_modifications", "return", "return"),
        ("params/parsers/diann.py", "_extract_with_regex", "return", "return"),
        ("params/parsers/diann.py", "_find_cmdline", "return", "return"),
        ("params/parsers/diann.py", "_flag", "return", "return"),
        ("params/parsers/diann.py", "_predictors_library", "return", "return"),
        ("params/parsers/peaks.py", "_fdr", "return", "return"),
        ("params/parsers/peaks.py", "_mass_tolerance", "return", "return"),
        ("params/parsers/peaks.py", "_value", "return", "return"),
        ("params/parsers/spectronaut.py", "_extract_tolerances", "return", "return"),
        ("params/parsers/spectronaut.py", "_first_regex_value", "return", "return"),
        ("params/parsers/spectronaut.py", "_main_search_block", "return", "return"),
        ("params/parsers/spectronaut.py", "_tolerance_patterns", "return", "return"),
        ("params/parsers/spectronaut.py", "_value", "return", "return"),
        ("params/parsers/spectronaut.py", "_value_regex", "return", "return"),
        # Cyclopts omission is part of the external command-line schema.
        ("scripts/cli.py", "annotate", "parameter", "output"),
        ("scripts/cli.py", "convert", "parameter", "level"),
        ("scripts/cli.py", "proteobench", "parameter", "output"),
        ("scripts/cli.py", "summary_cmd", "parameter", "modality"),
        ("scripts/extract_raw_file_db.py", "_feature_count", "return", "return"),
        ("scripts/extract_raw_file_db.py", "select", "parameter", "module"),
    }
)

AUDITED_OBJECT_SIGNATURE_SLOTS: frozenset[SignatureSlot] = frozenset(
    {
        # Raw backend/vendor scalar narrowing.
        ("adapters/anndata/namespace.py", "parse_namespace", "parameter", "stored"),
        ("adapters/anndata/proteobench.py", "_require_anndata", "parameter", "target"),
        ("converters/_fragments.py", "_split_packed", "parameter", "value"),
        ("converters/assemble.py", "_format_charge", "parameter", "value"),
        ("params/parsers/maxquant.py", "_joined_text", "parameter", "value"),
        ("params/parsers/maxquant.py", "_text", "parameter", "value"),
        ("serialization.py", "to_json_compatible", "parameter", "value"),
        # Pydantic before-validators receive untrusted values by contract.
        ("params/model.py", "MassTolerance.parse", "parameter", "value"),
        ("params/model.py", "Parameters._canonicalize_enzyme", "parameter", "value"),
        ("params/model.py", "Parameters._canonicalize_enzyme", "return", "return"),
        ("params/model.py", "Parameters._coerce_bool", "parameter", "value"),
        ("params/model.py", "Parameters._coerce_bool", "return", "return"),
        ("params/model.py", "Parameters._coerce_modifications", "parameter", "value"),
        ("params/model.py", "Parameters._coerce_modifications", "return", "return"),
        ("params/model.py", "Parameters._coerce_non_negative_float", "parameter", "value"),
        ("params/model.py", "Parameters._coerce_non_negative_float", "return", "return"),
        ("params/model.py", "Parameters._coerce_non_negative_int", "parameter", "value"),
        ("params/model.py", "Parameters._coerce_non_negative_int", "return", "return"),
        ("params/model.py", "Parameters._coerce_positive_int", "parameter", "value"),
        ("params/model.py", "Parameters._coerce_positive_int", "return", "return"),
        ("params/model.py", "Parameters._coerce_probability", "parameter", "value"),
        ("params/model.py", "Parameters._coerce_probability", "return", "return"),
        ("params/model.py", "Parameters._coerce_scan_window", "parameter", "value"),
        ("params/model.py", "Parameters._coerce_scan_window", "return", "return"),
        ("params/model.py", "Parameters._coerce_tolerance", "parameter", "value"),
        ("params/model.py", "Parameters._coerce_tolerance", "return", "return"),
        ("params/model.py", "Parameters._coerce_unparsed", "parameter", "value"),
        ("params/model.py", "Parameters._coerce_unparsed", "return", "return"),
        ("params/model.py", "Parameters._empty_strings_to_none", "parameter", "value"),
        ("params/model.py", "Parameters._empty_strings_to_none", "return", "return"),
        ("params/model.py", "Probability._coerce_value", "parameter", "value"),
        ("params/model.py", "Probability._coerce_value", "return", "return"),
        ("params/model.py", "Probability.parse", "parameter", "value"),
        ("params/model.py", "_coerce_float", "parameter", "value"),
        ("params/model.py", "_is_missing", "parameter", "value"),
        ("params/model.py", "_is_unknown_literal", "parameter", "value"),
        ("params/model.py", "_modification_from_item", "parameter", "item"),
        ("params/model.py", "_series_json_value", "parameter", "value"),
        ("params/model.py", "_to_scalar", "parameter", "value"),
    }
)


def _tree(path: Path) -> ast.Module:
    return ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def _function(path: Path, name: str) -> ast.FunctionDef | ast.AsyncFunctionDef:
    matches = [
        node
        for node in ast.walk(_tree(path))
        if isinstance(node, ast.FunctionDef | ast.AsyncFunctionDef) and node.name == name
    ]
    assert len(matches) == 1
    return matches[0]


def _container_import(node: ast.Import | ast.ImportFrom) -> bool:
    if isinstance(node, ast.Import):
        roots = {alias.name.split(".", 1)[0] for alias in node.names}
    else:
        roots = {node.module.split(".", 1)[0]} if node.module else set()
    return bool(roots & {"anndata", "mudata"})


def _function_annotations(
    node: ast.FunctionDef | ast.AsyncFunctionDef,
) -> list[ast.expr]:
    arguments = [*node.args.posonlyargs, *node.args.args, *node.args.kwonlyargs]
    annotations = [argument.annotation for argument in arguments if argument.annotation]
    if node.args.vararg and node.args.vararg.annotation:
        annotations.append(node.args.vararg.annotation)
    if node.args.kwarg and node.args.kwarg.annotation:
        annotations.append(node.args.kwarg.annotation)
    if node.returns:
        annotations.append(node.returns)
    return annotations


def _is_explicit_optional(annotation: ast.expr) -> bool:
    return any(
        isinstance(candidate, ast.BinOp)
        and isinstance(candidate.op, ast.BitOr)
        and any(
            isinstance(item, ast.Constant) and item.value is None for item in ast.walk(candidate)
        )
        for candidate in ast.walk(annotation)
    )


def _is_generic_object_property_bag(node: ast.AST) -> bool:
    if not isinstance(node, ast.Subscript) or not isinstance(node.value, ast.Name):
        return False
    if node.value.id not in {"dict", "Mapping", "MutableMapping"}:
        return False
    if not isinstance(node.slice, ast.Tuple) or len(node.slice.elts) != 2:
        return False
    key, value = node.slice.elts
    return (
        isinstance(key, ast.Name)
        and key.id == "str"
        and isinstance(value, ast.Name)
        and value.id == "object"
    )


def _is_unparameterized_array(annotation: ast.expr) -> bool:
    """Bare ``np.ndarray``/``NDArray`` is ``ndarray[Any, dtype[Any]]``."""
    if isinstance(annotation, ast.Name):
        return annotation.id == "NDArray"
    return isinstance(annotation, ast.Attribute) and annotation.attr == "ndarray"


def _contains_unparameterized_array(annotation: ast.expr) -> bool:
    """True when an array annotation carries no dtype, making it implicitly ``Any``."""
    subscripted = {
        id(node.value) for node in ast.walk(annotation) if isinstance(node, ast.Subscript)
    }
    return any(
        isinstance(node, ast.Name | ast.Attribute)
        and _is_unparameterized_array(node)
        and id(node) not in subscripted
        for node in ast.walk(annotation)
    )


def _contains_object(annotation: ast.expr) -> bool:
    return any(
        isinstance(candidate, ast.Name) and candidate.id == "object"
        for candidate in ast.walk(annotation)
    )


class _SignatureSlotVisitor(ast.NodeVisitor):
    """Collect matching parameter/return slots with stable qualified names."""

    def __init__(self, relative_path: str, predicate: AnnotationPredicate) -> None:
        self.relative_path = relative_path
        self.predicate = predicate
        self.scope: list[str] = []
        self.slots: set[SignatureSlot] = set()

    @override
    def visit_ClassDef(self, node: ast.ClassDef) -> None:
        self.scope.append(node.name)
        self.generic_visit(node)
        self.scope.pop()

    @override
    def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
        self._visit_function(node)

    @override
    def visit_AsyncFunctionDef(self, node: ast.AsyncFunctionDef) -> None:
        self._visit_function(node)

    def _visit_function(self, node: ast.FunctionDef | ast.AsyncFunctionDef) -> None:
        qualified_name = ".".join((*self.scope, node.name))
        arguments = [*node.args.posonlyargs, *node.args.args, *node.args.kwonlyargs]
        if node.args.vararg is not None:
            arguments.append(node.args.vararg)
        if node.args.kwarg is not None:
            arguments.append(node.args.kwarg)
        for argument in arguments:
            if argument.annotation is not None and self.predicate(argument.annotation):
                self.slots.add((self.relative_path, qualified_name, "parameter", argument.arg))
        if node.returns is not None and self.predicate(node.returns):
            self.slots.add((self.relative_path, qualified_name, "return", "return"))

        self.scope.append(node.name)
        self.generic_visit(node)
        self.scope.pop()


def _production_signature_slots(predicate: AnnotationPredicate) -> frozenset[SignatureSlot]:
    slots: set[SignatureSlot] = set()
    for path in PACKAGE.rglob("*.py"):
        visitor = _SignatureSlotVisitor(path.relative_to(PACKAGE).as_posix(), predicate)
        visitor.visit(_tree(path))
        slots.update(visitor.slots)
    return frozenset(slots)


SCIENTIFIC_CONSTRUCTORS = frozenset(
    {
        "DataFrame",
        "Series",
        "array",
        "asarray",
        "coo_matrix",
        "csc_array",
        "csc_matrix",
        "csr_array",
        "csr_matrix",
    }
)


def _cli_scientific_manipulation_offenders(tree: ast.Module) -> list[str]:
    dataframe_names = {
        argument.arg
        for function in ast.walk(tree)
        if isinstance(function, ast.FunctionDef | ast.AsyncFunctionDef)
        for argument in (*function.args.posonlyargs, *function.args.args, *function.args.kwonlyargs)
        if argument.annotation and ast.unparse(argument.annotation) in {"DataFrame", "pd.DataFrame"}
    }
    dataframe_names.update(
        target.id
        for assignment in ast.walk(tree)
        if isinstance(assignment, ast.Assign)
        and isinstance(assignment.value, ast.Call)
        and isinstance(assignment.value.func, ast.Name)
        and assignment.value.func.id == "read_table"
        for target in assignment.targets
        if isinstance(target, ast.Name)
    )

    offenders: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute) and node.attr in SCIENTIFIC_CONTAINER_SLOTS:
            offenders.append(f"container slot {node.attr}:{node.lineno}")
        if (
            isinstance(node, ast.Subscript)
            and isinstance(node.value, ast.Name)
            and node.value.id in dataframe_names
        ):
            offenders.append(f"DataFrame subscript {node.value.id}:{node.lineno}")
        if not isinstance(node, ast.Call):
            continue
        if isinstance(node.func, ast.Attribute):
            if isinstance(node.func.value, ast.Name) and node.func.value.id in {"np", "pd"}:
                offenders.append(f"scientific call {ast.unparse(node.func)}:{node.lineno}")
            if (
                isinstance(node.func.value, ast.Name)
                and node.func.value.id in dataframe_names
                and node.func.attr in DATAFRAME_MANIPULATION_METHODS
            ):
                offenders.append(f"DataFrame call {ast.unparse(node.func)}:{node.lineno}")
        elif isinstance(node.func, ast.Name) and node.func.id in SCIENTIFIC_CONSTRUCTORS:
            offenders.append(f"scientific constructor {node.func.id}:{node.lineno}")
    return offenders


def test_explicit_optional_detector_reaches_nested_unions() -> None:
    annotation = ast.parse("tuple[str | None, ...]", mode="eval").body
    assert _is_explicit_optional(annotation)


def test_production_optional_signatures_match_reviewed_external_boundaries() -> None:
    observed = _production_signature_slots(_is_explicit_optional)
    assert observed == AUDITED_OPTIONAL_SIGNATURE_SLOTS


def test_production_object_signatures_match_reviewed_raw_boundaries() -> None:
    observed = _production_signature_slots(_contains_object)
    assert observed == AUDITED_OBJECT_SIGNATURE_SLOTS


@pytest.mark.parametrize("path", PURE_FILES, ids=lambda path: path.relative_to(PACKAGE).as_posix())
def test_computation_modules_do_not_import_storage_backends(path: Path) -> None:
    offenders = [
        node.lineno
        for node in ast.walk(_tree(path))
        if isinstance(node, ast.Import | ast.ImportFrom) and _container_import(node)
    ]
    assert offenders == []


def test_storage_backend_imports_are_confined_to_adapter_and_cli() -> None:
    offenders: list[str] = []
    for path in PACKAGE.rglob("*.py"):
        relative = path.relative_to(PACKAGE)
        permitted = relative.parts[:2] == ("adapters", "anndata") or relative == Path(
            "scripts/cli.py"
        )
        if permitted:
            continue
        for node in ast.walk(_tree(path)):
            if isinstance(node, ast.Import | ast.ImportFrom) and _container_import(node):
                offenders.append(f"{relative}:{node.lineno}")
    assert offenders == []


@pytest.mark.parametrize("path", PURE_FILES, ids=lambda path: path.relative_to(PACKAGE).as_posix())
def test_computation_signatures_are_concrete(path: Path) -> None:
    offenders: list[str] = []
    for node in ast.walk(_tree(path)):
        if not isinstance(node, ast.FunctionDef | ast.AsyncFunctionDef):
            continue
        for annotation in _function_annotations(node):
            contains_any = any(
                isinstance(item, ast.Name) and item.id == "Any" for item in ast.walk(annotation)
            )
            if contains_any:
                offenders.append(f"{node.name}:{node.lineno}: Any")
            if _is_explicit_optional(annotation):
                offenders.append(f"{node.name}:{node.lineno}: {ast.unparse(annotation)}")
    assert offenders == []


@pytest.mark.parametrize(
    "path",
    sorted((PACKAGE / "workflows").glob("*.py")),
    ids=lambda path: path.name,
)
def test_workflows_do_not_access_container_slots(path: Path) -> None:
    offenders = [
        f"{node.attr}:{node.lineno}"
        for node in ast.walk(_tree(path))
        if isinstance(node, ast.Attribute) and node.attr in CONTAINER_SLOTS
    ]
    assert offenders == []


def test_cli_composition_does_not_manipulate_scientific_frames_or_matrices() -> None:
    path = PACKAGE / "scripts" / "cli.py"
    assert _cli_scientific_manipulation_offenders(_tree(path)) == []


@pytest.mark.parametrize("path", WORKFLOW_FILES, ids=lambda path: path.name)
def test_workflows_do_not_manipulate_scientific_frames_or_matrices(path: Path) -> None:
    """Orchestration orders calculations; it never reshapes a frame or a matrix itself."""
    assert _cli_scientific_manipulation_offenders(_tree(path)) == []


def test_production_signatures_have_no_unparameterized_arrays() -> None:
    """Bare ``np.ndarray``/``NDArray`` is ``Any`` in both axes; require a dtype."""
    assert _production_signature_slots(_contains_unparameterized_array) == frozenset()


def test_unparameterized_array_detector_distinguishes_dtyped_annotations() -> None:
    bare = ast.parse("np.ndarray", mode="eval").body
    dtyped = ast.parse("NDArray[np.float64]", mode="eval").body
    nested = ast.parse("dict[str, np.ndarray]", mode="eval").body
    assert _contains_unparameterized_array(bare)
    assert not _contains_unparameterized_array(dtyped)
    assert _contains_unparameterized_array(nested)


def test_dependency_direction_keeps_adapters_below_workflows() -> None:
    """Only the composition root may import the storage adapter."""
    offenders: list[str] = []
    for path in PACKAGE.rglob("*.py"):
        relative = path.relative_to(PACKAGE)
        if relative.parts[0] == "adapters" or relative == Path("scripts/cli.py"):
            continue
        for node in ast.walk(_tree(path)):
            if not isinstance(node, ast.Import | ast.ImportFrom):
                continue
            names = (
                [alias.name for alias in node.names]
                if isinstance(node, ast.Import)
                else [node.module or ""]
            )
            if any(name.startswith("anndata_proteomics.adapters") for name in names):
                offenders.append(f"{relative}:{node.lineno}")
    assert offenders == []


def test_namespace_key_is_named_in_exactly_one_module() -> None:
    """One module owns the APB ``uns`` namespace key."""
    offenders = [
        str(path.relative_to(PACKAGE))
        for path in PACKAGE.rglob("*.py")
        if path.relative_to(PACKAGE).parts[0] == "adapters"
        and path.relative_to(PACKAGE) != Path("adapters/anndata/namespace.py")
        and '"anndata_proteomics"' in path.read_text(encoding="utf-8")
    ]
    assert offenders == []


def test_cli_scientific_manipulation_guard_detects_direct_operations() -> None:
    tree = ast.parse(
        """
def invalid(data: DataFrame, target: object) -> None:
    data.copy()
    data["feature"]
    target.X
    pd.merge(data, data)
"""
    )

    offenders = _cli_scientific_manipulation_offenders(tree)

    assert len(offenders) == 4
    assert any(offender.startswith("DataFrame call data.copy") for offender in offenders)
    assert any(offender.startswith("DataFrame subscript data") for offender in offenders)
    assert any(offender.startswith("container slot X") for offender in offenders)
    assert any(offender.startswith("scientific call pd.merge") for offender in offenders)


def test_fasta_adapter_has_no_mixed_calculation_and_persistence_entrypoints() -> None:
    path = PACKAGE / "adapters" / "anndata" / "fasta.py"
    definitions = {
        node.name
        for node in ast.walk(_tree(path))
        if isinstance(node, ast.FunctionDef | ast.AsyncFunctionDef)
    }
    assert definitions.isdisjoint(MIXED_FASTA_ADAPTER_ENTRYPOINTS)


def test_fasta_adapter_has_no_empty_string_level_sentinel() -> None:
    """Absence of a stored level is a presence check, never an empty string."""
    functions = {
        node.name
        for node in ast.walk(_tree(PACKAGE / "adapters" / "anndata" / "fasta.py"))
        if isinstance(node, ast.FunctionDef)
    }
    assert "_namespace_level" not in functions
    assert "_has_namespace_level" in functions
    assert "_require_namespace_level" in functions


def test_production_has_no_typing_any_or_suppressions() -> None:
    offenders: list[str] = []
    for path in PACKAGE.rglob("*.py"):
        source = path.read_text(encoding="utf-8")
        tree = ast.parse(source, filename=str(path))
        if any(isinstance(node, ast.Name) and node.id == "Any" for node in ast.walk(tree)):
            offenders.append(f"{path.relative_to(PACKAGE)}: Any")
        if "# noqa" in source or "type: ignore" in source:
            offenders.append(f"{path.relative_to(PACKAGE)}: suppression")
    assert offenders == []


def test_production_has_no_broad_exception_handlers() -> None:
    offenders: list[str] = []
    for path in PACKAGE.rglob("*.py"):
        for node in ast.walk(_tree(path)):
            if not isinstance(node, ast.ExceptHandler):
                continue
            if node.type is None or (
                isinstance(node.type, ast.Name) and node.type.id in {"BaseException", "Exception"}
            ):
                offenders.append(f"{path.relative_to(PACKAGE)}:{node.lineno}")
    assert offenders == []


@pytest.mark.parametrize(
    ("path", "function_name"),
    NO_EXCEPTION_CONTROL_FLOW,
    ids=lambda value: value.name if isinstance(value, Path) else value,
)
def test_parsing_decisions_do_not_use_exceptions_as_control_flow(
    path: Path,
    function_name: str,
) -> None:
    handlers = [
        node.lineno
        for node in ast.walk(_function(path, function_name))
        if isinstance(node, ast.ExceptHandler)
    ]
    assert handlers == []


@pytest.mark.parametrize(
    "function_name",
    ("extract_gene_name", "parse_header_id"),
)
def test_fasta_absence_results_do_not_use_empty_strings(function_name: str) -> None:
    function = _function(PACKAGE / "fasta" / "annotation.py", function_name)
    empty_strings = [
        node.lineno
        for node in ast.walk(function)
        if isinstance(node, ast.Constant) and node.value == ""
    ]
    assert empty_strings == []


def test_description_primary_metadata_fields_are_concrete_contracts() -> None:
    tree = _tree(PACKAGE / "description.py")
    metadata = next(
        node
        for node in tree.body
        if isinstance(node, ast.ClassDef) and node.name == "DescriptionMetadata"
    )
    annotations = {
        ast.unparse(node.target): ast.unparse(node.annotation)
        for node in metadata.body
        if isinstance(node, ast.AnnAssign)
    }
    assert annotations["quantification_level"] == "DescriptionQuantificationLevel"
    assert annotations["software_name"] == "DescriptionSoftwareName"
    assert annotations["conversion"] == "DescriptionConversionMetadata"


@pytest.mark.parametrize(
    ("adapter", "forbidden"),
    tuple(FORBIDDEN_ADAPTER_ENTRY_POINTS.items()),
    ids=lambda value: value if isinstance(value, str) else "entry-points",
)
def test_adapters_do_not_import_scientific_entry_points(
    adapter: str,
    forbidden: frozenset[str],
) -> None:
    """A storage adapter persists a calculated result; it never runs the calculation."""
    offenders = [
        f"{alias.name}:{node.lineno}"
        for node in ast.walk(_tree(PACKAGE / adapter))
        if isinstance(node, ast.ImportFrom)
        for alias in node.names
        if alias.name in forbidden
    ]
    assert offenders == []


def test_production_has_no_generic_object_property_bags() -> None:
    offenders: list[str] = []
    for path in PACKAGE.rglob("*.py"):
        offenders.extend(
            f"{path.relative_to(PACKAGE)}:{node.lineno}"
            for node in ast.walk(_tree(path))
            if isinstance(node, ast.Subscript) and _is_generic_object_property_bag(node)
        )
    assert offenders == []


def test_production_does_not_trigger_validation_with_a_null_sentinel() -> None:
    offenders: list[str] = []
    for path in PACKAGE.rglob("*.py"):
        for node in ast.walk(_tree(path)):
            if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
                continue
            if node.func.attr != "validate_python" or not node.args:
                continue
            first = node.args[0]
            if isinstance(first, ast.Constant) and first.value is None:
                offenders.append(f"{path.relative_to(PACKAGE)}:{node.lineno}")
    assert offenders == []
