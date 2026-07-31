"""Keep the rules, params, and modifications subpackages in an acyclic order."""

from __future__ import annotations

import ast
from collections.abc import Iterator
from pathlib import Path

from anndata_proteomics.modifications import schema as modification_schema
from anndata_proteomics.rules import schema as rule_schema

_SOURCE_ROOT = Path(__file__).parents[1] / "src" / "anndata_proteomics"
_LAYERS = ("modifications", "params", "rules")


def _imported_modules(source: Path) -> Iterator[str]:
    tree = ast.parse(source.read_text(encoding="utf-8"), filename=str(source))
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            yield from (alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            yield node.module


def test_package_dependency_layers_are_acyclic() -> None:
    """A package may import itself or a lower layer, never a higher layer."""
    violations: list[str] = []
    for layer_index, layer in enumerate(_LAYERS):
        allowed = set(_LAYERS[: layer_index + 1])
        for source in sorted((_SOURCE_ROOT / layer).rglob("*.py")):
            for module in _imported_modules(source):
                prefix = "anndata_proteomics."
                if not module.startswith(prefix):
                    continue
                imported_package = module.removeprefix(prefix).partition(".")[0]
                if imported_package in _LAYERS and imported_package not in allowed:
                    relative_source = source.relative_to(_SOURCE_ROOT)
                    violations.append(f"{relative_source}: {module}")

    assert violations == []


def test_layer_packages_have_empty_initializers() -> None:
    for layer in ("modifications", "params"):
        assert (_SOURCE_ROOT / layer / "__init__.py").read_bytes() == b""


def test_rule_schema_preserves_modification_model_imports() -> None:
    """Existing rules.schema imports remain aliases of the owning models."""
    names = (
        "ModificationMapEntry",
        "Modifications",
        "SiteListModifications",
        "TokenPosition",
        "TokenRegexModifications",
        "UnknownPolicy",
    )
    for name in names:
        assert getattr(rule_schema, name) is getattr(modification_schema, name)
