"""Load and validate a parsing-rule TOML into a ParseRule."""

from __future__ import annotations

import re
import tomllib
from pathlib import Path

from pydantic import ValidationError

from anndata_proteomics.rules.registry import _VERSION_FOLDER_RE, find_rule, resolve_rule_path
from anndata_proteomics.rules.schema import ParseRule


def _is_array_of_tables(value: list) -> bool:
    """True if a TOML array holds tables (``[[...]]``), e.g. layers/compute/map entries.

    A non-empty array whose every element is a table — distinguishes ``[[layers]]`` (append on
    merge) from scalar arrays like ``obs_keys`` / ``var_keys`` (replace on merge).
    """
    return bool(value) and all(isinstance(item, dict) for item in value)


def _merge_rule_dicts(base: dict, leaf: dict) -> dict:
    """Deep-merge a leaf rule dict onto its vendor ``base`` (leaf wins).

    - scalars: the leaf value replaces the base value;
    - tables (dicts, e.g. ``[columns.var.select]``, ``[axis]``): deep-merged, leaf keys win;
    - arrays of tables (``[[layers]]``, ``[[columns.var.compute]]``, ``[[modifications.map]]``):
      base entries first, then the leaf's appended — this preserves compute dependency order;
    - scalar arrays (``obs_keys``, ``var_keys``): the leaf value replaces the base value.
    """
    merged = dict(base)
    for key, leaf_val in leaf.items():
        base_val = merged.get(key)
        if isinstance(base_val, dict) and isinstance(leaf_val, dict):
            merged[key] = _merge_rule_dicts(base_val, leaf_val)
        elif isinstance(base_val, list) and isinstance(leaf_val, list):
            if _is_array_of_tables(base_val) or _is_array_of_tables(leaf_val):
                merged[key] = base_val + leaf_val
            else:
                merged[key] = leaf_val
        else:
            merged[key] = leaf_val
    return merged


def _vendor_base_path(leaf_path: Path) -> Path | None:
    """The vendor base file ``<vendor>/<vendor>.toml`` for a leaf rule, or None.

    Convention-based inheritance: a leaf at ``<vendor>/parse_*.toml`` or
    ``<vendor>/v*/parse_*.toml`` is merged onto its vendor base. The vendor directory is the
    leaf's parent (skipping a ``v*`` version folder). Returns None when no base file exists
    (the leaf then loads standalone) or when the path *is* the base.

    Note: a vendor directory must not be named like a version folder (``v1``, ``v2_3``, …); such
    a name would be mistaken for a version subfolder and its base skipped. No packaged vendor is.
    """
    parent = leaf_path.parent
    vendor_dir = parent.parent if _VERSION_FOLDER_RE.match(parent.name) else parent
    base = vendor_dir / f"{vendor_dir.name}.toml"
    return base if base.exists() and base != leaf_path else None


def load_rule(path: Path | str) -> ParseRule:
    """Load a TOML file and validate it as a ParseRule.

    If a vendor base file (``<vendor>/<vendor>.toml``) sits alongside the rule, its shared
    blocks are merged in first (see :func:`_merge_rule_dicts`); this is the single choke point,
    so ``recognize`` / ``validate`` / ``convert`` all see the merged rule.

    Raises FileNotFoundError if the path doesn't exist.
    On pydantic validation failure, the file path is attached as an exception note
    so it shows up in the traceback alongside pydantic's field-level message.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    data = tomllib.loads(p.read_text())
    base_path = _vendor_base_path(p)
    if base_path is not None:
        data = _merge_rule_dicts(tomllib.loads(base_path.read_text()), data)
    try:
        return ParseRule.model_validate(data)
    except ValidationError as e:
        e.add_note(f"in {p}")
        raise


def load_packaged_rule(
    software: str, quantification_level: str, version: str | None = None
) -> ParseRule:
    """Load the packaged rule for (software, level) at ``version`` (None → version-agnostic root)."""
    rule = load_rule(find_rule(software, quantification_level, version))
    if version is not None and not _software_version_matches(rule, version):
        raise ValueError(
            f"rule software_version={rule.software_version!r} does not match "
            f"parsed version={version!r}"
        )
    return rule


def resolve_rule_for_version(
    software: str, quantification_level: str, version: str | None
) -> ParseRule | None:
    """The rule whose version subfolder covers ``version``, or ``None`` if no variant applies.

    DIA-NN report columns vary by version; this picks the right variant by the software version
    (parsed from the param file), via ``registry.resolve_rule_path``. ``None`` means the level
    is not available for that version (e.g. fragment on DIA-NN 2.x).
    """
    path = resolve_rule_path(software, quantification_level, version)
    if path is None:
        return None
    rule = load_rule(path)
    if version is not None and not _software_version_matches(rule, version):
        return None
    return rule


def _software_version_matches(rule: ParseRule, version: str) -> bool:
    """True when ``rule.software_version`` regex matches a parsed parameter version."""
    try:
        return re.search(rule.software_version, version) is not None
    except re.error as exc:
        raise ValueError(f"invalid software_version regex {rule.software_version!r}") from exc
