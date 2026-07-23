"""Low-level discovery primitives shared by rule loading and registry APIs."""

from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
from pathlib import Path

from anndata_proteomics.rules.schema import QuantificationLevel


@dataclass(frozen=True, order=True, slots=True)
class RuleLocator:
    """Operational address of one level inside a software-version document."""

    path: Path
    level: QuantificationLevel


def packaged_rules_root() -> Path:
    """Return the filesystem path to the packaged parsing-rule documents."""
    traversable = resources.files("anndata_proteomics") / "parsing_rules"
    return Path(str(traversable))


def iter_packaged_documents() -> tuple[Path, ...]:
    """Return self-contained JSON documents in stable path order."""
    root = packaged_rules_root()
    paths = set(root.glob("*/rules.json")) | set(root.glob("*/v*/rules.json"))
    return tuple(sorted(paths))


def document_vendor(path: Path | str) -> str:
    """Return the vendor folder slug for a packaged document path."""
    resolved = Path(path).resolve()
    root = packaged_rules_root().resolve()
    try:
        return resolved.relative_to(root).parts[0]
    except (ValueError, IndexError) as exc:
        raise ValueError(f"not a packaged rule document: {resolved}") from exc


def document_paths_for_software(software: str) -> tuple[Path, ...]:
    """Return every packaged version document for a vendor slug."""
    vendor = packaged_rules_root() / software
    if not vendor.is_dir():
        return ()
    paths = set(vendor.glob("rules.json")) | set(vendor.glob("v*/rules.json"))
    return tuple(sorted(paths))
