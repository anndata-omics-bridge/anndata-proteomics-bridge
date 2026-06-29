# Review: parsing-rule TOMLs and modification handling

Date: 2026-06-27 (findings re-verified against the working tree on the same day)

## Status: ALL FINDINGS IMPLEMENTED (2026-06-27)

- **Finding 3** — removed the unimplemented `already_proforma` / `separate_mod_column` parser
  modes (schema, pipeline, `docs/toml_schema.md`, `docs/parsing_architecture.md` + regenerated
  `.html`, regenerated `parse_rule.schema.json`); only `token_regex` remains.
- **Finding 4 + 5** — `UnimodEntry.target` / `MapEntry.target` are now `list[str]`;
  `_target_matches` is list-aware (residue + terminus + empty guard); `UNIMOD:21` → `[S,T,Y]`;
  added `UNIMOD:121` (GG on K, 114.04293). Duplicate-accession guard kept; `sdrf.py` / `model.py`
  intentionally unchanged (those targets are the params/SDRF side, not the registry).
- **Finding 1** — convention-based base/leaf inheritance: merge engine in `rules.loader.load_rule`
  (single choke point); vendor bases `diann/diann.toml` + `spectronaut/spectronaut.toml`; the 7
  DIA-NN/Spectronaut leaves stripped of shared blocks. Base files are excluded from the
  `parse_*.toml` glob (still 11 packaged rules).

**Verification:** full suite green (364 passed / 4 skipped, ruff clean); a before/after fingerprint
proved converted output (obs/var/X/layers + MuData) is **byte-identical** for both vendors across
all resolved levels. An adversarial-review pass (18 issues → fixes) added: an `apply_modifications`
guard so protein levels inherit `[modifications]` but skip applying it when no compute consumes it
(restores exact pre-refactor protein behavior, closes a recognize/KeyError gap and the per-row perf
cost), git-tracked the base files, regenerated the rendered HTML, and strengthened the
duplicate-accession / merge-override / parser-mode / fragment-inheritance tests.

Scope:

- APB parsing rules under `src/anndata_proteomics/parsing_rules/`
- APB modification schema/runtime: `src/anndata_proteomics/rules/schema.py`,
  `src/anndata_proteomics/modifications/` (`pipeline.py`, `apply_rules.py`, `unimod_registry.py`,
  `unimod_registry.toml`)
- ProteoBench as upstream behaviour evidence (`origin/main`, `origin/intermediate_format_interface`)

Every APB finding below was checked against the current code; file:line anchors are given so the
next person can confirm without re-deriving. None of these overlap the recent `apb convert` / `apb
list` / README reconciliation — modifications are a separate subsystem — but findings 1 and 3 have
**documentation knock-ons** (README + `docs/toml_schema.md`) that are called out so the docs don't
drift the way the CLI docs had.

## Findings

### 1. Vendor `[modifications]` blocks are duplicated verbatim across level TOMLs — VERIFIED

The `[modifications]` blocks are byte-for-byte identical between a vendor's levels:

- DIA-NN: [`diann/parse_diann_ion.toml`](../src/anndata_proteomics/parsing_rules/diann/parse_diann_ion.toml) ≡ [`diann/v1/parse_diann_fragment.toml`](../src/anndata_proteomics/parsing_rules/diann/v1/parse_diann_fragment.toml)
- Spectronaut: [`spectronaut/parse_spectronaut_ion_1.toml`](../src/anndata_proteomics/parsing_rules/spectronaut/parse_spectronaut_ion_1.toml) ≡ [`spectronaut/parse_spectronaut_fragment.toml`](../src/anndata_proteomics/parsing_rules/spectronaut/parse_spectronaut_fragment.toml)

(Confirmed by diffing the extracted blocks — identical in both pairs.) Eight rules carry a
`[modifications]` block today; all eight use `parser = "token_regex"`.

The duplicated content is *not* the canonical Unimod metadata — that is already centralised in
[`modifications/unimod_registry.toml`](../src/anndata_proteomics/modifications/unimod_registry.toml),
and per-rule entries reference it by accession only (see `ModificationMapEntry`,
[schema.py:102](../src/anndata_proteomics/rules/schema.py#L102)). What duplicates is the **vendor
parser profile**: `source_column`, `token_pattern`, `token_position`, `case_sensitive`,
`unknown_policy`, and the vendor-token→accession map.

Recommendation: a **convention-based base/leaf hierarchy** per vendor — the duplicated blocks move
into a vendor base file and the level TOMLs inherit them. This dedups not just `[modifications]` but
every block shared across a vendor's levels (obs.select, shared var.select, shared computes). See
**"Plan: convention-based rule inheritance"** at the end for the full design. (Supersedes an earlier
idea of a profile registry keyed by id; a base/leaf hierarchy covers more duplication with **no
invented TOML keys** — TOML has no `extends`/include of its own, so inheritance is inferred from
directory layout, not declared in the files.)

### 2. ProteoBench duplicates the same parser blocks — useful as evidence, not as a model

On `origin/intermediate_format_interface`, ProteoBench repeats identical modification-parser blocks
across many per-module/per-tool TOMLs (DIA-NN, Spectronaut, PEAKS, FragPipe each in several files).
This confirms the *behaviour* APB must match, but ProteoBench's copy-paste factoring is **not** the
structure APB should adopt — finding 1's base/leaf hierarchy is the APB-side answer. (Exact per-tool
file counts depend on the ProteoBench branch/commit and are not load-bearing for this conclusion.)

### 3. Schema + docs advertise parser modes the runtime silently no-ops — VERIFIED

The schema defines a three-way discriminated union on `parser`
([schema.py:154](../src/anndata_proteomics/rules/schema.py#L154)):

- `token_regex` (`TokenRegexModifications`)
- `already_proforma` (`AlreadyProformaModifications`)
- `separate_mod_column` (`SeparateModColumnModifications`)

But the runtime implements only the first —
[`pipeline.apply_modifications`](../src/anndata_proteomics/modifications/pipeline.py) opens with
`if mods.parser != "token_regex": return df`, a **silent no-op**. A TOML using either other mode
passes `apb validate` but produces no `proforma_sequence` / `stripped_sequence` columns, so it fails
later when `var_keys` / computes reference the missing column. The schema docstring and the pipeline
docstring both admit token_regex is "the only parser with a runtime implementation today" — but the
contract is still validate-OK / convert-broken.

This is made worse by **docs**: [docs/toml_schema.md:263-285](../docs/toml_schema.md#L263) documents
`already_proforma` and `separate_mod_column` as if they were usable, with no "not implemented" note.

No shipped rule uses the unimplemented modes (all 8 are `token_regex`), so the cheapest correct fix
is to **remove `already_proforma` and `separate_mod_column` from the schema and the docs until a rule
needs them** — reintroduce alongside the runtime + a test when one lands. (Alternative: implement
them now; ProteoBench's `proforma.py` already has helpers for ProForma-like inputs, e.g. I2MassChroQ,
and for sequence-plus-separate-mod-column inputs, if a near-term rule justifies it.) A silent no-op
is the one option to rule out.

### 4. The Unimod registry cannot represent one accession on multiple residues — VERIFIED

[`unimod_registry.py`](../src/anndata_proteomics/modifications/unimod_registry.py) loads the TOML into
`{accession: UnimodEntry}` and **raises on a duplicate accession** (`load_registry`,
[unimod_registry.py:44](../src/anndata_proteomics/modifications/unimod_registry.py#L44)); `UnimodEntry.target`
is a single `str`. So an accession is locked to exactly one residue target.

Today the registry holds 6 entries — `UNIMOD:1, 4, 21, 27, 28, 35`. `UNIMOD:21` (Phospho) targets
**`S` only**, and `UNIMOD:121` (GG) is absent — even though APB's *parameter* parser and ProteoBench
upstream both recognise phospho on S/T/Y and GG on K. Because `apply_rules` matches a token by its
adjacent residue against the entry's single `target`
([apply_rules.py `_target_matches`](../src/anndata_proteomics/modifications/apply_rules.py)), a
`UniMod:21` token on T or Y resolves to nothing and is preserved as an unknown token.

Recommendation: before adding phospho/GG to any parsing rule, change the registry model so one
accession can carry multiple residue targets — either `target: list[str]` (e.g. `["S","T","Y"]`) or
multiple records per accession resolved by token context. Then add `UNIMOD:121` and extend tests to
cover `UNIMOD:21` on S, T, and Y.

### 5. DIA-NN / Spectronaut ion+fragment token maps are not currently missing a ProteoBench token — VERIFIED

The DIA-NN and Spectronaut ion/fragment maps cover the same three common mods — Acetyl
(`UNIMOD:1`), Carbamidomethyl (`UNIMOD:4`), Oxidation (`UNIMOD:35`) — which is exactly what the
current ProteoBench DIA-NN/Spectronaut ion parse-settings cover. So these specific token maps are
**not** behind ProteoBench today. The drift risk lives in the shared registry/runtime/schema
(findings 3-4), not in these two maps. (FragPipe additionally uses `UNIMOD:27`/`UNIMOD:28`.)

## Suggested implementation order

Ordered cheapest-correctness-first, dependencies before dependents:

1. **Close the validate-OK / convert-broken gap (finding 3).** Remove `already_proforma` and
   `separate_mod_column` from the schema **and** from `docs/toml_schema.md` until the runtime + a
   real rule need them. Safe now — no rule uses them. Re-add as a unit: schema + `pipeline.py` +
   test, together.
2. **Make the registry multi-target (finding 4).** `target: list[str]` (or multi-record). Keep the
   "raise on true duplicate" guard but allow multiple residues per accession.
3. **Add `UNIMOD:121` and S/T/Y phospho coverage**, with tests for `UNIMOD:21` on S, T, Y. Depends
   on step 2.
4. **De-duplicate the level TOMLs (finding 1)** via the convention-based base/leaf hierarchy — see
   the plan below. Independent of steps 1-3; can be done first if dedup is the priority. Doc
   knock-ons: README §"Adding a new conversion" / §"ProForma sequences & modifications",
   `docs/toml_schema.md`, and the `[columns.*.select]` note in `AGENTS.md`.
5. **Keep level TOMLs focused** on quantification-level structure once the base exists: axes,
   selected metadata, layers, fragments, level-specific computes — no copied shared blocks.

Do not adopt ProteoBench's per-module copy-paste factoring (finding 2).

## Plan: convention-based rule inheritance (addresses Finding 1)

**Decided design — no invented TOML keys.** TOML has no `extends`/include of its own; inheritance is
inferred from directory layout by the loader, so every TOML file stays 100% standard.

- **Base file:** `parsing_rules/<vendor>/<vendor>.toml` (e.g. `diann/diann.toml`) holds content
  shared by all of that vendor's levels. It does **not** match the `parse_*.toml` glob in
  `iter_packaged_rules`, so it is never treated as a rule, listed by `apb list`, or converted
  directly — it is only ever merged into a leaf.
- **Leaf rules** keep today's names and locations: `parse_<vendor>_<level>.toml` at the vendor root
  (= all versions) or in `v*/` (= version-specific). Version-folder resolution
  (`resolve_rule_path`) is unchanged.
- **Inheritance by convention:** loading a leaf implicitly merges it onto its vendor base
  (`<vendor>/<vendor>.toml`, found by walking up from the leaf, skipping any `v*/` folder). No base
  file present → the leaf loads standalone (backward compatible). No key in any file.

**Merge semantics** (child = leaf overrides parent = base):

- scalars (`x_layer`, `quantification_level`, `software_version`, …): child replaces base;
- tables (`[columns.obs.select]`, `[columns.var.select]`, `[axis]`, `[axis.duplicates]`):
  deep-merge, child keys win;
- arrays of tables (`[[columns.var.compute]]`, `[[layers]]`, `[[modifications.map]]`): **base
  entries first, child appended** — preserves compute dependency order (base defines
  `ProForma_peptidoform`; leaf appends `ProForma_ion` which references it);
- sub-objects (`[modifications]`, `[fragments]`): inherited whole unless the leaf defines its own;
- append-only for now — no remove/replace directive (no leaf needs to drop a base entry today; add
  one later if that changes).

**Single choke point:** `rules/loader.load_rule` resolves the base, merges, then validates the
**merged** dict as the existing `ParseRule`. `recognize`, `apb validate`, and `convert` all consume
merged rules unchanged; `parse_rule.schema.json` (the merged shape) is unchanged.

**Base/leaf split** (DIA-NN, from the diffs):

- `diann/diann.toml` (base): `software_name`, `input_shape`, `schema_version`, `[axis].obs_keys` +
  `[axis.duplicates]`, `[columns.obs.select]`, the shared `[columns.var.select]` columns, the shared
  computes (`ProForma_peptidoform`, `ProForma_peptide`), `[modifications]`.
- each leaf: `quantification_level`, `software_version`, `file_version`, `[axis].var_keys` +
  `x_layer`, level-specific `[columns.var.select]` columns, the level compute (`ProForma_ion` /
  `ProForma_fragment`), `[[layers]]`, `[fragments]`.

**Scope:** DIA-NN and Spectronaut only — the vendors with ≥2 leaves, where duplication exists.
Single-rule vendors (MaxQuant, FragPipe, PEAKS, WOMBAT) stay single files; the optional-base
convention means no base is created for them (splitting would add a file without removing any
duplication). **Two-level** (vendor base → leaf); add a mid-level base only if a vendor ever needs
version-wide-but-not-level-specific overrides.

**Phases:**

1. Loader merge engine + unit tests on fixtures (merge a synthetic base+leaf; assert the documented
   scalar/table/array/sub-object semantics). No real-file moves yet.
2. Confirm `apb validate` routes through the merging `load_rule` (a leaf that relies on the base for
   required fields like `layers` must still validate); fix if it parses raw TOML.
3. Create `diann/diann.toml` + `spectronaut/spectronaut.toml`; strip the shared blocks out of their
   leaves. Prove no behavioural change: convert the existing DIA-NN/Spectronaut test data before and
   after and diff the resulting `.h5ad`/`.h5mu` (plus the full pytest run).
4. Docs: document the convention in `docs/toml_schema.md`, README §"Adding a new conversion", and
   the `[columns.*.select]` note in `AGENTS.md`. No `parse_rule.schema.json` change.
