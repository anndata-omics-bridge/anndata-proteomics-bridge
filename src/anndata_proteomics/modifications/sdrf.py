"""SDRF-Proteomics modification rendering."""

from __future__ import annotations

from anndata_proteomics.modifications.model import ModType, SearchedModification


def to_sdrf_value(mod: SearchedModification) -> str:
    """Render a :class:`SearchedModification` as an SDRF key=value string.

    Output uses the canonical SDRF-Proteomics order:
    ``NT=<name>;AC=<accession>;MT=<fixed|variable>;TA=<target>;PP=<position>``.

    Empty fields are omitted (the SDRF spec accepts subsets when values are
    unknown).
    """
    parts: list[str] = [f"NT={mod.name}"]
    if mod.accession:
        parts.append(f"AC={mod.accession}")
    if mod.mod_type is not ModType.unknown:
        parts.append(f"MT={mod.mod_type.value}")
    if mod.target:
        parts.append(f"TA={mod.target}")
    if mod.position:
        parts.append(f"PP={mod.position}")
    return ";".join(parts)


def from_sdrf_value(value: str) -> SearchedModification:
    """Parse an SDRF ``comment[modification parameters]`` value.

    Order-insensitive. Unknown keys are silently dropped.
    """
    fields: dict[str, str] = {}
    for token in value.split(";"):
        if "=" not in token:
            continue
        key, raw = token.split("=", 1)
        fields[key.strip().upper()] = raw.strip()

    if "NT" not in fields:
        raise ValueError(f"SDRF modification missing NT field: {value!r}")

    mt = ModType.unknown
    if "MT" in fields:
        try:
            mt = ModType(fields["MT"].lower())
        except ValueError as exc:
            raise ValueError(f"unknown MT value: {fields['MT']!r}") from exc

    return SearchedModification(
        name=fields["NT"],
        accession=fields.get("AC"),
        mod_type=mt,
        target=fields.get("TA"),
        position=fields.get("PP") or "Anywhere",
    )
