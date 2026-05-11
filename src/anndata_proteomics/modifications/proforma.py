"""ProForma sequence rendering from modification occurrences."""

from __future__ import annotations

from anndata_proteomics.modifications.model import ModificationOccurrence


def render_proforma(
    stripped: str,
    occurrences: list[ModificationOccurrence],
    unknown_tokens: dict[int, str] | None = None,
) -> str:
    """Build a ProForma 2.0 string from a stripped sequence + modifications.

    Parameters
    ----------
    stripped
        Unmodified amino acid sequence.
    occurrences
        Localized modifications. ``sequence_index`` is 0-based into
        ``stripped``; ``position`` may be ``"N-term"`` / ``"C-term"`` for
        terminal modifications (then ``sequence_index`` is ignored).
    unknown_tokens
        Optional mapping ``{sequence_index: original_vendor_token}`` for
        unresolved tokens. Index ``-1`` denotes N-term, ``len(stripped)``
        denotes C-term. Rendered verbatim inside brackets.

    Notes
    -----
    Mods at the same residue are concatenated (``M[Oxidation][Acetyl]``).
    Preferred label per occurrence: accession when present
    (``[UNIMOD:35]``), else name (``[Oxidation]``).
    """
    unknown_tokens = unknown_tokens or {}
    nterm: list[str] = []
    cterm: list[str] = []
    by_residue: dict[int, list[str]] = {}

    for occ in occurrences:
        tag = occ.accession or occ.name
        if occ.position == "N-term":
            nterm.append(tag)
        elif occ.position == "C-term":
            cterm.append(tag)
        elif occ.sequence_index is not None:
            by_residue.setdefault(occ.sequence_index, []).append(tag)

    for idx, token in unknown_tokens.items():
        if idx == -1:
            nterm.append(token)
        elif idx == len(stripped):
            cterm.append(token)
        else:
            by_residue.setdefault(idx, []).append(token)

    out: list[str] = []
    if nterm:
        out.append("[" + "][".join(nterm) + "]-")
    for i, residue in enumerate(stripped):
        out.append(residue)
        if i in by_residue:
            out.append("[" + "][".join(by_residue[i]) + "]")
    if cterm:
        out.append("-[" + "][".join(cterm) + "]")
    return "".join(out)
