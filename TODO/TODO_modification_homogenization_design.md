# Design Sketch — Centralized Modification Homogenization for Parameter Parsing

**Status:** proposal / sketch (no code yet). Companion to the ProteoBench-sync finding and `TODO_code_review_june.md` (ANTI-002).

## Problem

APB's vendor parsers hand raw `fixed_mods` / `variable_mods` strings straight into
`Parameters`, where [`_coerce_modifications`](../src/anndata_proteomics/params/model.py)
merely comma-splits them into `SearchedModification(name=str(item))`. ProteoBench `main`
(modification-homogenization work, merged Jun 2026) now emits a canonical
**`{residue|terminal}[{ModName}]`** form, comma-joined, for every vendor:

| Vendor input | ProteoBench output |
|---|---|
| MaxQuant `Oxidation (M),Acetyl (Protein N-term)` | `M[Oxidation], Protein N-term[Acetyl]` |
| MaxQuant `Phospho (STY)` | `S[Phospho], T[Phospho], Y[Phospho]` |
| FragPipe mass `57.02146` on `C` | `C[Carbamidomethyl]` |
| Sage `{"C": 57.0215}` | `C[Carbamidomethyl]` |
| AlphaPept `cC`, `oxM` | `C[Carbamidomethyl]`, `M[Oxidation]` |

Result: every APB vendor's modification fields mismatch `main`. This sketch replaces the
raw passthrough with one normalizer + small per-vendor tables, reusing APB's existing
modification infrastructure rather than copying ProteoBench's per-file dicts.

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

## Sequencing (scope-aware)

AGENTS.md current scope is ion-level **DIA-NN, MaxQuant, Spectronaut** — do these first.

1. `ParamModResolver` (+ mass index) and `render_param_mods`; extend `unimod_registry.toml`
   to cover the common mods (Carbamidomethyl, Oxidation, Acetyl, Phospho, pyro-Glu, …).
2. `maxquant_tokens` (also serves Spectronaut + MSAID) → wire MaxQuant + Spectronaut.
3. DIA-NN code/name alias table → wire DIA-NN. **At this point the 3 in-scope vendors sync.**
4. Mass-based: `fragpipe_tokens`, `sage_tokens` (resolve via `by_mass`).
5. Remaining: AlphaPept (codes), Wombat (`"of"`), MetaMorpheus (`"on"`).
6. **Test:** compare `render_param_mods(...)` against ProteoBench's `test/params/*.csv`
   expected `fixed_mods`/`variable_mods` columns — closes the current false-confidence gap
   (APB tests skip mod fields today, e.g. Sage).

## Out of scope / follow-ups

- This addresses parameter-file mod *lists*. The result-table modified-*sequence*
  homogenization (`apply_rules`/`pipeline`) is already handled and unchanged.
- ANTI-002 (typed `MassTolerance` re-serialization adding `.0`) is a separate sync fix.
- Track ProteoBench PR #1012 (`parse_input()` intermediate-format API) — the longer-term
  direction may let APB consume ProteoBench's homogenized output directly instead of
  re-deriving it; this normalizer is the bridge until then.
