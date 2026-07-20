---
name: proforma
description: Use when working in APB on ProForma peptide, peptidoform, ion, or fragment identifiers; vendor modified-sequence normalization; parsing-rule TOML modifications or fragments blocks; Unimod token mapping; or code under anndata_proteomics/modifications.
---

# ProForma

## Purpose

Use this skill in APB when adding, fixing, or reviewing ProForma-related behavior:
vendor modified-sequence columns, peptidoform and ion identifiers, fragment feature
keys, `[modifications]` TOML blocks, `[fragments]` TOML blocks, and Unimod registry
mappings.

APB currently implements a small, practical ProForma normalization subset. Do not
silently widen it into a general ProForma parser or add a new public API surface.

## First Read

Before editing ProForma behavior, read the local source of truth:

- `AGENTS.md` for APB-specific rules.
- `docs/toml_schema.md`, especially "Var-axis naming convention", "Modifications",
  and "Fragments".
- `src/anndata_proteomics/rules/schema.py` for schema invariants.
- `src/anndata_proteomics/modifications/apply_rules.py` for token extraction,
  token-position handling, unknown-token policy, and mapping behavior.
- `src/anndata_proteomics/modifications/pipeline.py` for TOML-to-runtime conversion.
- `src/anndata_proteomics/modifications/proforma.py` for final rendering.
- `src/anndata_proteomics/modifications/unimod_registry.toml` before adding a
  new accession to a parsing rule.
- `tests/test_packaged_modifications.py`, `tests/test_rule_modifications.py`, and
  converter/schema tests touched by the change.

External reference: Pyteomics `pyteomics.proforma` documents ProForma parsing and
formatting, including `parse()`, `to_proforma()`, `ProForma.parse()`, chimeric
peptidoforms, charge state/adduct support, and controlled-vocabulary caching:
https://pyteomics.readthedocs.io/en/latest/api/proforma.html

Treat Pyteomics as a reference or validation aid unless the user explicitly asks to
add it as a dependency. APB does not currently depend on `pyteomics`, and adding it
would change dependency and runtime behavior.

## Pyteomics Validation

Use Pyteomics as an optional test oracle only when it directly checks APB's intended
surface:

- Good fit: smoke-validate generated `ProForma_peptidoform` strings and
  `ProForma_ion` strings with `pyteomics.proforma.ProForma.parse()`.
- Good fit: compare APB rendering of simple localized Unimod modifications against
  `pyteomics.proforma.to_proforma()` in an optional test.
- Possible fit: force controlled-vocabulary resolution by checking parsed tag `mass`
  values, not merely by parsing. Pyteomics may parse an unknown `UNIMOD:*` token
  without raising until resolution/mass access.
- Poor fit: validating APB `ProForma_fragment` values. APB currently uses
  `{peptidoform}/{charge}/{fragment_label}`, while Pyteomics parses the
  `{peptidoform}/{charge}` ion grammar and rejects the extra fragment-label segment.
- Poor fit: validating allowed residue/terminus placement. Pyteomics can parse
  syntactically valid but biologically mismatched placements such as an oxidation tag
  on a non-M residue; APB's target checks live in `apply_rules.py` plus
  `unimod_registry.toml`.

If adding this to tests, keep it behind an optional dependency such as
`pyteomics[proforma]`, not in APB's runtime dependencies. The ProForma extra brings
the Pyteomics ProForma support stack; plain `pyteomics` may be insufficient for
Unimod-backed parsing in a clean environment.

## APB Mental Model

APB uses schema-pinned names for ProForma-derived feature identifiers:

- `ProForma_peptide`: bare peptide sequence from `how = "stripped_sequence"`, for
  example `PEPTIDE`.
- `ProForma_peptidoform`: modified peptide sequence from
  `how = "proforma_sequence"`, for example `PEPM[UNIMOD:35]TIDE`.
- `ProForma_ion`: peptidoform plus charge from `how = "proforma_ion"`, for example
  `PEPM[UNIMOD:35]TIDE/2`.
- `ProForma_fragment`: ion plus fragment label from `how = "proforma_fragment"`,
  for example `PEPM[UNIMOD:35]TIDE/2/b4-unknown^1`.

Protein-level rules use protein identifiers and normally do not compute ProForma
columns.

APB docs currently describe the target as ProForma 2.0. Pyteomics documents a
ProForma v2.1 implementation. Do not change APB version claims or semantics without
checking the local docs and tests.

## Modification Workflow

When adding or fixing vendor modified-sequence conversion:

1. Start from real vendor examples and identify the exact source column containing
   modification tokens.
2. If the vendor has shared base TOML inheritance, put shared `[modifications]`
   blocks in the vendor base file. Put only level-specific behavior in leaf rules.
3. Keep vendor columns under `[columns.var.select]`. Put APB-derived columns under
   `[[columns.var.compute]]`; never select `proforma_sequence` or
   `stripped_sequence` as if they came from the input table.
4. Use only `parser = "token_regex"` unless the runtime implementation in
   `modifications/` is deliberately extended.
5. Make `token_pattern` capture the vendor token in its first capture group.
6. Set `token_position` from the real syntax, then verify the runtime behavior in
   `apply_rules.py`. The schema has several literals, but implementation details
   live in the runtime parser.
7. Map each vendor token to an accession with `[[modifications.map]]`. Canonical
   name, target, position, and mass delta come from `unimod_registry.toml`, not from
   per-vendor TOML.
8. Add missing accessions to `unimod_registry.toml` before referencing them.
9. Add or update packaged-rule tests with concrete vendor input and expected
   ProForma output.

Minimal TOML shape:

```toml
[[columns.var.compute]]
name = "ProForma_peptidoform"
from = ["Modified_Sequence"]
how = "proforma_sequence"

[[columns.var.compute]]
name = "ProForma_peptide"
from = ["Modified_Sequence"]
how = "stripped_sequence"

[[columns.var.compute]]
name = "ProForma_ion"
from = ["ProForma_peptidoform", "Precursor_Charge"]
how = "proforma_ion"

[modifications]
source_column = "Modified.Sequence"
parser = "token_regex"
token_pattern = "\\(([^()]*)\\)"
token_position = "after_residue"
unknown_policy = "preserve"
output_column = "proforma_sequence"

[[modifications.map]]
token = "UniMod:35"
accession = "UNIMOD:35"
```

## Fragment Workflow

For fragment-level conversion:

- Use `quantification_level = "fragment"`.
- Use `axis.var_keys = ["ProForma_fragment"]`.
- Compute `ProForma_peptidoform`, then intermediate `ProForma_ion`, then
  `ProForma_fragment`.
- Add a `[fragments]` block only when the vendor packs fragment values in parallel
  delimiter-separated columns that must be exploded before conversion.
- Use `[fragments].label_output` as the second source for
  `how = "proforma_fragment"`.
- Remember that `ProForma_ion` is an intermediate at fragment level and should not be
  the var key there.

## Guardrails

- Fix the root cause in the upstream APB file: schema, TOML, registry, parser, or
  renderer. Do not add wrapper normalization, try/except fallbacks, or skip logic
  unless the user explicitly asks for a temporary workaround.
- Keep column names and compute names exactly consistent with `schema.py`.
- Do not add public methods or a new dependency just to handle one vendor token.
- Do not fake advanced ProForma features such as ambiguous localization, labile
  modifications, chimeric peptidoforms, adducts, global modifications, glycans, or
  cross-links with ad hoc string manipulation. If the task needs those semantics,
  first decide whether APB should depend on a real ProForma parser such as Pyteomics.

## Validation

Run focused tests from the APB repo root after edits:

```bash
pytest tests/test_packaged_modifications.py tests/test_rule_modifications.py
```

Also run the specific converter, schema, or rule-loading tests touched by the change.
If a dependency/environment issue prevents testing, report the exact command and
failure.
