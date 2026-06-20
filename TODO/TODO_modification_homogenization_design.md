# Plan — Centralized (registry-backed) Modification Homogenization for Parameter Parsing

**Status:** DEFERRED (2026-06-20). The *functional* homogenization is already done (the
ProteoBench-sync integration makes every vendor emit correct ProForma-style mod strings; code-review
Action 15 deduped the shared mechanics), and the current state is green and ProteoBench-matched.
This plan would add the **registry-backed architecture** — resolve every vendor mod to a *typed*
`SearchedModification` (accession + target + mass_delta from `unimod_registry`), collapsing the
per-vendor `MASS_TO_MOD` / alias dicts into one source of truth and enabling accession-backed SDRF
export.

**Why deferred:** nothing consumes that metadata yet. SDRF export (`to_sdrf_value`) and accession
resolution are wired only for the modified-*sequence* path (`apply_rules`/`pipeline`), **not** the
parameter-file path — so this would build a new module + resolver + per-vendor adapters + TOML alias
tables + registry expansion + a 10-vendor rewrite to populate fields no caller reads. Revisit when
SDRF export of search *parameters* is a real deliverable.

**Two caveats for whoever picks it up (see analysis below):**
1. The lower-risk first step is **choke-point enrichment** in `_coerce_modifications` (fill
   `accession`/`target`/`mass_delta` from the registry *after* the string is built, leaving the
   rendered `name`/`source` untouched) — not the full per-vendor rewrite. The byte-for-byte mod
   tests then guarantee parity because the rendered string never changes.
2. Collapsing FragPipe/Sage `MASS_TO_MOD` into the registry **changes rendered names**
   (`Pyro-glu` → `Glu->pyro-Glu`; `GG` and `Label:*` are absent from the 6-entry registry), so it
   needs vendor-name overrides to keep the ProteoBench match. "Delete the mass dict" is not free.

Companion to the archived ProteoBench-sync plan and `TODO_code_review_june.md` (Actions 4 + 15).
The design below is retained as the reference for when this is revisited.

## Current state (2026-06-20) — the baseline this plan builds on

**Already in place (do _not_ redo):**

- **Functional homogenization works.** Every vendor's `extract_params` emits `fixed_mods` /
  `variable_mods` as comma-joined ProForma-style strings (`C[Carbamidomethyl], M[Oxidation]`),
  matching ProteoBench. The per-vendor mod fields are now compared in `tests/test_params_*.py`
  (e.g. `test_params_sage`, `test_params_diann`) — the old "tests skip mod fields" gap is **closed**.
- **Shared mechanics extracted (Action 15).** `params/_common.homogenize_paren_mods(mod, mapping)`
  (MaxQuant + Spectronaut) and `lookup_mass_mod(mass, mapping, *, tol=)` (FragPipe + Sage) already
  centralize the *tokenizing* and *mass-match* mechanics. MetaMorpheus (` on `) and WOMBAT (` of `)
  keep their own single-use homogenizers.
- **Typed tolerances (Action 4) done** — the ANTI-002 follow-up below is resolved.

**Still raw / unbuilt (the gap this plan closes):**

- **No registry resolution.** `Parameters._coerce_modifications` splits the homogenized string and
  stores `SearchedModification(name="C[Carbamidomethyl]", source=...)` — `accession`, `target`,
  `mass_delta` stay `None`. `SearchedModification` *has* those fields (name/accession/mod_type/
  target/position/mass_delta/source); the param path simply never fills them, so there is no
  accession-backed metadata for SDRF/QC.
- **Per-vendor data dicts persist.** `MODIFICATION_MAPPING` (diann/peaks/alphapept/msaid/...),
  `_MASS_TO_MOD` (fragpipe) and `MASS_TO_MOD_MAPPING` (sage) still hold mass/name data the registry
  could own. `unimod_registry.toml` has only **6 entries** (Acetyl, Carbamidomethyl, Phospho,
  Glu->pyro-Glu, Gln->pyro-Glu, Oxidation) and exposes `resolve(accession)` only — no `by_name` /
  `by_code` / `by_mass`.
- The token `name` is the whole `residue[Name]` string, not the bare mod name.

So this plan = **resolve, don't re-derive**: feed the already-tokenized vendor mods through a
registry-backed resolver to produce *typed* `SearchedModification`s, and collapse the per-vendor
mass/name dicts into the registry + small alias tables.

## Problem

Both APB **and** ProteoBench `main` now emit a canonical **`{residue|terminal}[{ModName}]`**
form, comma-joined, for every vendor (the ProteoBench-sync integration added the per-vendor
mappers, so the strings now match):

| Vendor input | ProteoBench output |
|---|---|
| MaxQuant `Oxidation (M),Acetyl (Protein N-term)` | `M[Oxidation], Protein N-term[Acetyl]` |
| MaxQuant `Phospho (STY)` | `S[Phospho], T[Phospho], Y[Phospho]` |
| FragPipe mass `57.02146` on `C` | `C[Carbamidomethyl]` |
| Sage `{"C": 57.0215}` | `C[Carbamidomethyl]` |
| AlphaPept `cC`, `oxM` | `C[Carbamidomethyl]`, `M[Oxidation]` |

The rendered strings now match `main`. **What's still wrong is downstream:**
`_coerce_modifications` comma-splits that string into `SearchedModification(name=str(item))`, so
`name` is the whole `residue[Name]` token and `accession` / `target` / `mass_delta` are never
resolved. This plan replaces that untyped split with one registry-backed normalizer + small
per-vendor tables, reusing APB's existing modification infrastructure rather than keeping
ProteoBench's per-file dicts.

## Core insight

There are three vendor **input classes**, but one **target representation** and one
**source of truth**:

- **Source of truth (reuse, don't duplicate):** `modifications.unimod_registry` already
  carries `name`, `mass_delta`, `target`, `position`, `accession` per accession. This
  single registry subsumes ProteoBench's three scattered dicts:
  - `MODIFICATION_MAPPING` (name fallback) → registry name + alias table
  - `MASS_TO_MOD` (FragPipe) / `MASS_TO_MOD_MAPPING` (Sage) → registry `mass_delta` index
  - `RESIDUE_MAP` → registry `target`

- **Three input classes → three resolution strategies → one renderer:**

  | Class | Vendors | Strategy |
  |---|---|---|
  | name + residue (syntactic) | MaxQuant, Spectronaut, MSAID, MetaMorpheus, Wombat | reshape `Name (res)` / `Name of res` / `Name on res`, expand multi-residue, keep name verbatim |
  | code-mapped | AlphaPept, DIA-NN short codes | alias `code → accession` |
  | mass-based | FragPipe, Sage | `mass_delta → accession` (±tolerance) |

  All three converge on a typed `SearchedModification`, rendered as `residue[Name]`.

## Architecture

New module **`modifications/param_mods.py`** (sits beside the existing sequence-mod infra;
the existing `apply_rules`/`pipeline`/`proforma` handle *result-table modified sequences*,
this handles *parameter-file mod lists* — same target form, different input, so a sibling
module, not a fork).

```
                       ┌─────────────────────────────────────────┐
  raw vendor string ──▶│  per-vendor ADAPTER  (the only per-vendor │
  "Oxidation (M),..."  │  code: tokenize → RawModToken[])          │
                       └───────────────────┬──────────────────────┘
                                           ▼
                       RawModToken(name|code|mass, residue, terminal)
                                           ▼
                       ┌─────────────────────────────────────────┐
                       │  ParamModResolver  (registry-backed)      │
                       │   .by_name(alias) .by_code(c) .by_mass(d) │  ◀── unimod_registry (REUSE)
                       └───────────────────┬──────────────────────┘
                                           ▼
                       SearchedModification(name, accession, target, mass_delta, source)
                                           ▼
                       render_param_mods(...) ──▶ "M[Oxidation], Protein N-term[Acetyl]"
```

### Pieces

1. **`ParamModResolver`** (registry-backed; one instance, cached):
   ```python
   class ParamModResolver:
       def by_accession(self, acc: str) -> UnimodEntry | None: ...
       def by_name(self, alias: str) -> UnimodEntry | None:   # alias table → accession
       def by_code(self, code: str) -> UnimodEntry | None:    # vendor code → accession
       def by_mass(self, delta: float, tol: float = 1e-3) -> UnimodEntry | None:
   ```
   `by_mass` builds a one-time sorted mass index over `load_registry()` (replaces
   FragPipe's `MASS_TO_MOD` + `MASS_TOLERANCE = 0.001` and Sage's `MASS_TO_MOD_MAPPING`).

2. **Per-vendor adapter** — the *only* per-vendor code. A tiny callable per vendor:
   ```python
   def maxquant_tokens(raw: str) -> list[RawModToken]: ...   # "Name (residues)", expand STY, terminals
   def fragpipe_tokens(raw: str) -> list[RawModToken]: ...   # mass,residue,active,sites table
   def sage_tokens(raw: dict) -> list[RawModToken]: ...      # {residue: delta | [deltas]}
   # Spectronaut, MSAID reuse maxquant_tokens; Wombat/MetaMorpheus = "of"/"on" variants
   ```
   `RawModToken = (key: str|float, residue: str|None, terminal: str|None, kind: Literal["name","code","mass"])`.

3. **Per-vendor mapping tables** — declarative, in TOML (per AGENTS.md: "migrate
   modification mapping rules into APB parsing TOMLs/schema"). Two small tables:
   - `code → accession` (AlphaPept `cC`/`oxM`/`a<^`; DIA-NN short forms)
   - `vendor name alias → accession` (only where a vendor name ≠ registry name)
   Mass-based vendors need **no table** — the registry `mass_delta` index covers them.

4. **Orchestrator:**
   ```python
   def homogenize_param_mods(raw, *, vendor: str, mod_type: ModType) -> list[SearchedModification]:
       tokens = _ADAPTERS[vendor](raw)
       return [_resolve(tok, mod_type) for tok in tokens]
   ```

5. **Renderer** (small, dedicated — `proforma.render_proforma` solves the *different*
   full-peptide problem and is not reused here):
   ```python
   def render_param_mods(mods: list[SearchedModification]) -> str:
       # "M[Oxidation], Protein N-term[Acetyl]"  — joins target[name], comma+space
   ```

### Wiring

- Each vendor's `extract_params` calls `homogenize_param_mods(raw, vendor=..., mod_type=...)`
  (the vendor is known there) and passes `list[SearchedModification]` to `Parameters(...)`.
- `_coerce_modifications` keeps the string-split branch as a **typed fallback only**
  (already-homogenized `SearchedModification`s pass through unchanged).
- Serialization (`to_series` / `_serialize_modifications`) renders via `render_param_mods`
  so the CSV matches ProteoBench byte-for-byte.

## The key decision: verbatim-match vs registry-canonicalization

ProteoBench keeps the **vendor's** mod name verbatim (`C[Carbamidomethyl]`), reshaping only
syntax. APB's registry has *canonical* names + accessions. These can disagree.

- **Option A — sync-first (verbatim):** render the vendor name as-is; use the registry
  only to *derive* names for mass/code vendors. Byte-for-byte match with ProteoBench CSVs.
- **Option B — canonical-first:** always map to the registry's canonical name/accession.
  More correct, enables accession-backed SDRF export, but may diverge from ProteoBench
  strings where vendor name ≠ registry name.

**Recommendation — layered (A's string, B's metadata).** Resolve to a `SearchedModification`
that carries **both**: `name` = vendor-verbatim (drives the ProteoBench-compatible render),
`accession`/`target`/`mass_delta` = resolved from the registry (metadata for SDRF/QC),
`source` = the raw token. Default render uses `name` → stays in sync now; accession is
enrichment that never changes the rendered string. This reuses `SearchedModification`'s
existing fields exactly as designed and avoids choosing between sync and correctness.

## Reuse ledger (per AGENTS.md "name the canonical tool first")

| Need | Reuse | New |
|---|---|---|
| canonical name/mass/target source | `modifications.unimod_registry` (+ extend TOML beyond 8 entries) | mass index over it |
| typed mod record | `SearchedModification`, `ModType` | — |
| mass→name lookup | registry `mass_delta` | `ParamModResolver.by_mass` |
| name/code aliases | — | per-vendor TOML alias tables |
| tokenizing vendor strings | — | per-vendor adapters |
| `residue[Name]` rendering | — (not `render_proforma`; different problem) | `render_param_mods` (small) |

No per-vendor `MASS_TO_MOD` / `MODIFICATION_MAPPING` / `RESIDUE_MAP` dicts are copied from
ProteoBench — all of that data collapses into the one registry + two small alias tables.

**Already partly there (Action 15):** `_common.homogenize_paren_mods` and `lookup_mass_mod` are
the tokenizing / mass-match primitives this plan builds the per-vendor adapters and `by_mass` on
top of. Today they return *strings* over per-vendor dicts; this plan repoints them at the registry
and at typed `SearchedModification` output (the adapters return `RawModToken`s instead of joined
strings, and `by_mass` indexes `load_registry()` instead of `_MASS_TO_MOD`).

## Sequencing (from the 2026-06-20 baseline)

AGENTS.md current scope is ion-level **DIA-NN, MaxQuant, Spectronaut** — do these first. One
commit per step, tree green at each; `extract_params` signatures stay stable. Because the rendered
strings already match ProteoBench today, **every step is regression-guarded by the existing
mod-field comparisons** — they must stay byte-for-byte identical while the *types* underneath
gain accession/target/mass_delta.

1. **Registry + resolver + renderer.** Extend `unimod_registry.toml` past the current 6 entries to
   cover what the vendor dicts reference: add **GG** (`UNIMOD:121`), the SILAC **`Label:*`** set
   FragPipe carries (`Label:2H(4)`, `Label:13C(6)`, `Label:13C(6)15N(2)`, `Label:13C(6)15N(4)`),
   and confirm the two pyro-Glu accessions already present cover FragPipe's `-17.0265`/`-18.0106`.
   Add `ParamModResolver` with `by_accession` (wraps the existing `resolve`), `by_name`, `by_code`,
   and `by_mass` (one-time sorted mass index over `load_registry()`, reusing `lookup_mass_mod`'s
   `1e-3` tolerance). Add `render_param_mods(mods) -> "residue[Name], ..."`.
2. **Typed param path.** Add `homogenize_param_mods(raw, *, vendor, mod_type) ->
   list[SearchedModification]`. `_coerce_modifications` already passes `SearchedModification`
   instances through unchanged, so vendors can hand it typed lists; serialize via
   `render_param_mods` (or keep `source` = rendered token) so the CSV round-trip stays identical.
3. **MaxQuant + Spectronaut.** Refactor `_common.homogenize_paren_mods` into a `maxquant_tokens`
   adapter returning `RawModToken[]` (name + residue/terminal, expand `STY`); resolve `by_name`.
   Spectronaut + MSAID reuse it.
4. **DIA-NN.** code/name alias table → resolve `by_code`/`by_name`. **3 in-scope vendors now typed.**
5. **Mass-based: FragPipe, Sage.** Replace `_MASS_TO_MOD` / `MASS_TO_MOD_MAPPING` + the
   `lookup_mass_mod` wrappers with `by_mass` over the registry; delete the per-vendor mass dicts.
6. **Remaining:** AlphaPept (codes), WOMBAT (` of `), MetaMorpheus (` on `).
7. **Verify:** the per-vendor mod-field comparisons already exist (regression guard for the
   strings); **add** assertions that the resolved `accession` / `mass_delta` are now populated on
   the in-scope vendors — the new capability this plan adds.

## Out of scope / follow-ups

- This addresses parameter-file mod *lists*. The result-table modified-*sequence*
  homogenization (`apply_rules`/`pipeline`) is already handled and unchanged.
- ANTI-002 (typed `MassTolerance`) is **done** (code-review Action 4): vendors with numeric values
  build `MassTolerance(...)` directly; this plan does the analogous thing for modifications.
- Track ProteoBench PR #1012 (`parse_input()` intermediate-format API) — the longer-term
  direction may let APB consume ProteoBench's homogenized output directly instead of
  re-deriving it; this normalizer is the bridge until then.
