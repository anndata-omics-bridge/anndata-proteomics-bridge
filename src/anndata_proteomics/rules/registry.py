"""Locate packaged parsing-rule TOMLs, resolving version subfolders by software version.

DIA-NN report columns change across versions, so DIA-NN rules whose columns are
version-dependent live in version subfolders (``diann/v1/``, ``diann/v2/`` … , finer
``diann/v1_9/`` when needed); version-agnostic levels stay at the vendor root. The folder name
is the selector: ``resolve_rule_path`` maps a parsed software version to the most-specific
covering folder. Vendors with a single format keep their flat ``<vendor>/parse_*.toml``.
"""

from __future__ import annotations

import re
from collections.abc import Iterator
from importlib import resources
from pathlib import Path

from packaging.version import InvalidVersion, Version


class RuleNotFound(LookupError):
    """Raised when (software, level, version) resolves to no packaged TOML."""


def packaged_rules_root() -> Path:
    """Filesystem path to the parsing_rules/ directory inside the package."""
    traversable = resources.files("anndata_proteomics") / "parsing_rules"
    return Path(str(traversable))


def iter_packaged_rules() -> Iterator[Path]:
    """Every packaged parse_*.toml — flat (``<vendor>/``) and version-foldered (``<vendor>/v*/``)."""
    root = packaged_rules_root()
    paths = set(root.glob("*/parse_*.toml")) | set(root.glob("*/v*/parse_*.toml"))
    yield from sorted(paths)


_VERSION_FOLDER_RE = re.compile(r"^v(\d+(?:_\d+)*)$")


def _folder_version(name: str) -> tuple[int, ...] | None:
    """``(1, 9)`` for ``v1_9``, ``(2,)`` for ``v2``; ``None`` if not a version folder."""
    match = _VERSION_FOLDER_RE.match(name)
    return tuple(int(part) for part in match.group(1).split("_")) if match else None


def _release_tuple(version: str | None) -> tuple[int, ...] | None:
    """``(1, 9, 2)`` from ``'1.9.2 Academia '``; ``None`` if empty/unparseable."""
    if not version:
        return None
    try:
        return Version(str(version).split()[0]).release
    except (InvalidVersion, IndexError):
        return None


def resolve_rule_path(software: str, level: str, version: str | None = None) -> Path | None:
    """Path to the packaged rule for (software, level) at a given software version.

    Picks the most-specific ``v*`` subfolder whose version is a prefix of ``version`` and that
    contains ``parse_<software>_<level>.toml``; otherwise the vendor-root file
    (``parse_<software>_<level>.toml`` or a legacy ``parse_<software>_<level>_<n>.toml``).
    Returns ``None`` when nothing matches (e.g. a version-specific level the file's version does
    not provide). ``version=None`` skips the folders and uses the root file only.
    """
    vendor = packaged_rules_root() / software
    if not vendor.exists():
        return None
    filename = f"parse_{software}_{level}.toml"

    release = _release_tuple(version)
    if release is not None:
        matched = [
            (fv, sub)
            for sub in vendor.iterdir()
            if sub.is_dir()
            and (fv := _folder_version(sub.name)) is not None
            and release[: len(fv)] == fv
        ]
        for _fv, sub in sorted(matched, key=lambda item: len(item[0]), reverse=True):
            candidate = sub / filename
            if candidate.exists():
                return candidate

    root_file = vendor / filename
    if root_file.exists():
        return root_file
    legacy = sorted(vendor.glob(f"parse_{software}_{level}_*.toml"))
    return legacy[0] if legacy else None


def find_rule(software: str, level: str, version: str | None = None) -> Path:
    """``resolve_rule_path(...)`` or raise :class:`RuleNotFound`."""
    path = resolve_rule_path(software, level, version)
    if path is None:
        raise RuleNotFound(
            f"no packaged rule for software={software!r} level={level!r} version={version!r}"
        )
    return path
