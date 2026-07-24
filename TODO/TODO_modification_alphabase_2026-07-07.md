# TODO (handoff) — Bring alphabase's modification assets into APB, and unify the modification TODOs

Date: 2026-07-07
Status: HANDOFF — scope + references only. Design and implementation happen in a **separate session**.

## Task

APB's [alphabase review](../../REVIEW_what_we_learn_from_alphapep_base.md) found three mature,
reusable modification assets in MannLabs `alphabase` (review items **1, 3, 4**). Pull them into
APB's modification handling, and **unify** this with the two existing modification TODOs so we end
up with one registry-backed design instead of three overlapping ones.

Do **not** start the unification here — this note only records what the future session must cover
and where to look.

## In scope — the three alphabase assets to adopt

1. **Curated `modification_mappings` tables (item 1).** alphabase ships a large, battle-tested
   vendor-token → canonical-mod table (MaxQuant token zoo `S(Phospho (STY))`/`pS`/…; the full
   MSFragger mass table; Dimethyl / TMT / iTRAQ / mTRAQ / DiLeu plex families). APB currently
   re-curates a handful of tokens per parsing-rule JSON file by hand. Treat alphabase's table as the
   expansion source for APB's registry + per-vendor alias tables.
   - Source: `related_work/alphabase/alphabase/constants/const_files/psm_reader.yaml` (`modification_mappings:`).

2. **Mass-token matching with tolerance (item 3).** alphabase matches mass modifications within
   `mod_mass_tol` (0.1 Da) after truncating to 4 dp (`mass_mapped_mods`). APB's parsing-rule path
   matches **exact strings** (e.g. FragPipe `token = "57.0215"`, and even mixes 4 dp and 6 dp in one
   file — brittle to precision drift across versions/configs).
   - Note: APB **already** has a tolerant mass matcher on the *params* side — `lookup_mass_mod(mass, mapping, *, tol=)` in `params/_common` — to reuse rather than reinvent.
   - APB side: `apb/src/anndata_proteomics/parsing_rules/fragpipe/rules.json`, level `ion`.

3. **`fixed_C57`-style fixed-modification handling (item 4).** When a fixed modification
   (classically Carbamidomethyl@C) is not written into the vendor's modified sequence, APB has **no
   re-insertion mechanism**, so a fixed-mod search can silently under-report it. alphabase carries a
   `fixed_C57` flag for exactly this. Confirm against a fixed-Carbamidomethyl MaxQuant fixture, then
   design the APB equivalent.
   - Source: `related_work/alphabase/alphabase/psm_reader/maxquant_reader.py` (`parse_mod_seq(..., fixed_C57)`).
   - APB side: `apb/src/anndata_proteomics/parsing_rules/maxquant/rules.json`, level `ion`.

## Fold in / unify with (the reason this is one design, not three)

- **[TODO_modification_homogenization_design.md](TODO_modification_homogenization_design.md)** —
  the registry-backed design for the **parameter-file** mod path (DEFERRED). Its spine is exactly
  what items 1/3/4 need: `unimod_registry.toml` as the single source of truth, "resolve, don't
  re-derive", collapse per-vendor mass/name dicts into registry + small alias tables. Items 1/3/4
  land mainly on the **parsing-rule / modified-sequence** path, but they share that registry + alias
  + mass-tolerance machinery — hence unify.
- **[TODO_proforma.md](Archive/TODO_proforma.md)** — the ProForma encoding note and the peptidoform-vs-
  parameter split. Its "Mapping tables in APB" section (parsing-rule `[[modifications.map]]` →
  `unimod_registry.toml`) is the exact surface alphabase's tables would expand, and its
  `unknown_policy` (preserve/drop/error) governs unmatched tokens.

**Unifying question for the design session:** one registry + alias/mass-index layer feeding *both*
the parameter path (homogenization TODO) and the parsing-rule path (ProForma TODO + items 1/3/4) —
without changing rendered strings that are currently ProteoBench-matched byte-for-byte.

## Reference documents

- [REVIEW_what_we_learn_from_alphapep_base.md](../../REVIEW_what_we_learn_from_alphapep_base.md) — items 1, 3, 4 with evidence.
- [TODO_modification_homogenization_design.md](TODO_modification_homogenization_design.md) — params-path registry design.
- [TODO_proforma.md](Archive/TODO_proforma.md) — ProForma peptidoform encoding + APB mapping tables.
- `related_work/alphabase/alphabase/constants/const_files/psm_reader.yaml` — `modification_mappings`, `mass_mapped_mods`, `mod_mass_tol`, `fixed_C57`.
- `related_work/alphabase/alphabase/psm_reader/modification_mapper.py`, `maxquant_reader.py` — how alphabase applies them.
- APB: `src/anndata_proteomics/modifications/unimod_registry.toml`, `params/_common.py` (`lookup_mass_mod`), the per-vendor `parse_*` JSON files.
