"""Typed access to the APB-owned AnnData ``uns`` namespace.

This is the only production module that names the namespace key. Readers come in two
shapes: :func:`read_namespace` decodes the whole payload for description and summary
views, while :func:`read_namespace_text` reads one stored JSON string without decoding
everything beside it. Writers use :func:`update_namespace`, which merges shallowly so
values a call does not touch are persisted exactly as they were stored.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from anndata import AnnData
from mudata import MuData

from anndata_proteomics.serialization import JsonObject, JsonValue, to_json_compatible

NAMESPACE = "anndata_proteomics"


@dataclass(frozen=True, slots=True)
class MissingNamespaceText:
    """No text payload is stored under one APB namespace key."""


MISSING_NAMESPACE_TEXT = MissingNamespaceText()


def read_namespace(target: AnnData | MuData) -> JsonObject:
    """Read and decode the complete APB namespace as a recursive JSON object."""
    return parse_namespace(target.uns.get(NAMESPACE))


def parse_namespace(stored: object) -> JsonObject:
    """Validate and decode one raw namespace value at a storage boundary."""
    if stored is None:
        return {}
    if not isinstance(stored, Mapping):
        raise TypeError(f"uns[{NAMESPACE!r}] must be a mapping")
    for key in stored:
        if not isinstance(key, str):
            raise TypeError(f"uns[{NAMESPACE!r}] keys must be strings")
    decoded = to_json_compatible(stored)
    if not isinstance(decoded, dict):
        raise TypeError(f"uns[{NAMESPACE!r}] must decode to a JSON object")
    return decoded


def has_namespace_key(target: AnnData | MuData, key: str) -> bool:
    """Return whether one APB namespace key is present, decoding nothing."""
    stored = target.uns.get(NAMESPACE)
    if stored is None:
        return False
    if not isinstance(stored, Mapping):
        raise TypeError(f"uns[{NAMESPACE!r}] must be a mapping")
    return key in stored


def read_namespace_text(target: AnnData | MuData, key: str) -> str | MissingNamespaceText:
    """Read one stored JSON string without decoding the rest of the namespace."""
    stored = target.uns.get(NAMESPACE)
    if stored is None:
        return MISSING_NAMESPACE_TEXT
    if not isinstance(stored, Mapping):
        raise TypeError(f"uns[{NAMESPACE!r}] must be a mapping")
    if key not in stored:
        return MISSING_NAMESPACE_TEXT
    value = stored[key]
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if not isinstance(value, str):
        raise TypeError(
            f"uns[{NAMESPACE!r}][{key!r}] must be a JSON string; got {type(value).__name__}"
        )
    return value


def write_namespace(target: AnnData | MuData, namespace: JsonObject) -> None:
    """Replace the APB namespace with an already validated payload."""
    target.uns[NAMESPACE] = namespace


def update_namespace(target: AnnData | MuData, updates: Mapping[str, JsonValue]) -> None:
    """Merge keys into the APB namespace, leaving untouched keys in place.

    The surviving keys are decoded through :func:`parse_namespace`, so the persisted
    payload stays uniformly JSON-typed rather than mixing decoded and raw backend values.
    """
    namespace = read_namespace(target)
    namespace.update(updates)
    write_namespace(target, namespace)
