# Plan: Integrate ProteoBench Param-Parsing Changes Into APB

**Status:** ✅ COMPLETE (2026-06-19). All 10 vendors integrated; full suite 237 passed / 0 failed / 0 skipped; `generate_report.py` ok=6. Fixtures copied into `tests/params/` (78 files), tests repointed and self-contained. Refactor phase (TODO_code_review_june.md) can now proceed on top. Implementation notes appended at end.
**Sequencing decision (already made):** merge ProteoBench changes FIRST, then run the
[TODO_code_review_june.md](TODO_code_review_june.md) refactors on top. Rationale: the
review items are behavior-preserving and touch the exact same files; refactoring first
would force a painful manual reconciliation against restructured code. Integrating first
lets us lock the combined behavior with the (refreshed) fixtures, then refactor with a
green test suite as the safety net.

## Why now

APB's `params/` package was ported from ProteoBench `proteobench/io/params/` on
2026-05-12 as a **standalone re-architecture** (pydantic `Parameters`/`MassTolerance` +
`_common.py`, no `proteobench` imports), verified byte-for-byte against ProteoBench
fixtures (see [Archive/PLAN_20260511_parameters-and-modifications.md](Archive/PLAN_20260511_parameters-and-modifications.md)).

Since that fork, ProteoBench `main` gained **~920 lines** of param-parsing changes
(`parameter-homogenization-and-filtering` merged 2026-06-04, `some-param-updates` merged
2026-06-16). APB's parsers are now based on a **pre-homogenization** snapshot and have
drifted from current ProteoBench behavior and fixtures.

Delta measured: `git diff dfb755e3..origin/main -- proteobench/io/params/` (ProteoBench).

## Scope

In-scope (vendors APB ships), with ProteoBench source-side line delta since fork:

| APB module | ProteoBench source | Δ lines | Notes |
|---|---|---:|---|
| `diann.py` | `diann.py` | +77 | regexes, tolerance, mod mapping, defaults |
| `fragpipe.py` | `fragger.py` | +95 | new `fragpipe_v23_noMBR` path, FDR/MBR |
| `maxquant.py` | `maxquant.py` | +62 | tolerance, version-gated mods |
| `metamorpheus.py` | `metamorpheus.py` | +102 | largest per-vendor change |
| `sage.py` | `sage.py` | +48 | |
| `peaks.py` | `peaks.py` | +22 | |
| `alphapept.py` | `alphapept.py` | +10 | |
| `msaid.py` | `msaid.py` | +10 | |
| `spectronaut.py` | `spectronaut.py` | +6 | |
| `wombat.py` | `wombat.py` | +6 | |
| `model.py` (validators) | `__init__.py` normalize layer | +293 | reconcile into pydantic, see below |

Out of scope: de novo parsers (casanovo/instanovo/contranovo/adanovo/etc.) and
quant tools APB does not ship (alphadia, i2masschroq, msangel, proline, quantms, maxdia)
unless explicitly added.

## Architecture mapping (key decision)

ProteoBench's June work added a **runtime normalization layer** on its dataclass in
`io/params/__init__.py`:
`_MISSING_SENTINELS`, `_ENZYME_MAP`, `_AUTO_CALIBRATION_SENTINELS`/`_LABEL`,
`_TOLERANCE_FIELDS`/`_FLOAT_FIELDS`/`_INT_FIELDS`, and `fill_none()` / `normalize()` /
`normalize_dataframe_columns()`.

APB already re-expresses most of this as **pydantic validators** in `params/model.py`
(`Probability`, `MassTolerance.parse`, charge/length/mz before-validators, `_is_missing`,
`_validate_ranges`, etc.).

**Decision: fold ProteoBench's *rules and data* into APB's existing pydantic layer; do
NOT port `normalize()`/`fill_none()` as a parallel system.** Porting the dataclass layer
wholesale would duplicate validation APB already performs (violates the repo's
reuse-before-duplicate and root-cause rules). Concretely, reconcile rule-by-rule:

Confirmed gaps / reconciliation points (more may surface during execution):
1. **Enzyme canonicalization** — ProteoBench added `_ENZYME_MAP`; APB has a bare
   `enzyme: str | None` with no canonicalization. → Add an enzyme-canonicalization
   helper (in `_common.py`) + a `model.py` before-validator, seeded from `_ENZYME_MAP`.
2. **FDR as percentage** — ProteoBench treats numeric FDR `>= 1` as a percentage
   (`/100`). APB's `Probability(value: float = Field(ge=0, le=1))` *rejects* numeric `1.0`.
   → Extend the `Probability._coerce_value` before-validator to divide numeric `>= 1` by
   100, matching ProteoBench, but only if a current fixture actually exercises it (confirm
   first — do not add speculative coercion).
3. **Auto-calibration tolerance sentinels** — ProteoBench maps `"dynamic"`, `"0 ppm"`,
   `"[-0.0 ppm, 0.0 ppm]"`, etc. → `"Automatic calibration"`. Confirm APB's
   `MassTolerance` "automatic" mode recognizes the same expanded sentinel set.
4. **Missing sentinels** — diff ProteoBench `_MISSING_SENTINELS` against APB `_is_missing`
   and align.

## Verification strategy

Fixtures changed too: `git diff dfb755e3..origin/main -- test/params/` = 59 files,
+1162/-228, incl. a new `fragpipe_v23_noMBR` fixture pair and edits to fragpipe/sage/
metamorpheus/wombat/mqpar expected CSVs.

- **Refresh APB's copied fixtures** from current ProteoBench `test/params/` (in-scope
  vendors only) before/with each vendor's port, so APB validates against current expected
  output, not the May snapshot.
- Per the May-11 plan's contract, APB parser serialized output must match the refreshed
  ProteoBench fixtures (byte-for-byte where practical).
- Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_params_*.py`
  plus `tests/test_rule_validate.py` and `tests/test_converters_e2e.py`.
- After pipeline-touching changes, run `uv run python tools/generate_report.py`.

## Segment plan (one commit per segment; tree green at each)

- **S0 — Baseline.** Run APB param tests as-is; record current pass/fail against the
  *existing* (stale) fixtures. Establish the diff baseline.
- **S1 — Model/normalization reconciliation.** Apply the 4 gap items above in
  `model.py`/`_common.py`. Unit tests for enzyme map, FDR%, auto-cal sentinels.
- **S2 — Refresh fixtures + port DIA-NN.** (current-scope vendor #1) Refresh DIA-NN
  fixtures, port the +77 of parser logic, green.
- **S3 — MaxQuant.** (current-scope #2) Refresh fixtures, port +62, green.
- **S4 — Spectronaut.** (current-scope #3) Refresh fixtures, port +6, green.
- **S5 — FragPipe.** Port +95 incl. new `fragpipe_v23_noMBR` fixture/path, green.
- **S6 — Remaining shipped vendors.** metamorpheus (+102), sage (+48), peaks (+22),
  alphapept (+10), msaid (+10), wombat (+6) — one commit each, each fixture-verified.
- **S7 — Full suite + report.** Run everything incl. `generate_report.py`; confirm no
  regressions. Hand off to the code-review refactors.

If a segment is larger than expected, split it; do not merge segments.

## Constraints

- APB stays free of `proteobench` imports (standalone, per the original migration
  contract). ProteoBench is **not** modified by this work.
- Keep `extract_params` signatures and public API unchanged (the code-review refactor
  depends on stable surfaces).
- Root-cause only: port logic into the correct APB module; no normalize-wrapper/tryCatch
  bandaids.

## Resolved decisions (2026-06-19)

1. **Vendor scope** — **all 10 shipped vendors** (DIA-NN, MaxQuant, Spectronaut,
   FragPipe, metamorpheus, sage, peaks, alphapept, msaid, wombat). S5–S6 are in this pass.
2. **Fixture refresh mechanism** — **copy ProteoBench `test/params/` fixtures into the APB
   tree** (snapshot). APB tests stay self-contained; re-sync on future ProteoBench changes.

## Implementation notes (2026-06-19)

**S1 — model.py (shared, done once):**
- Added `_ENZYME_MAP` (ported from ProteoBench) + a `_canonicalize_enzyme` before-validator on `enzyme`; removed `enzyme` from the generic empty-string validator.
- Extended `MassTolerance.parse` to treat ProteoBench's auto-calibration sentinels (`0 ppm`, `[-0.0 ppm, 0.0 ppm]`, numeric `0`, `0`) as automatic, collapsing the label to the canonical `"Automatic calibration"`.
- Extended `_MISSING_STRINGS` with `placeholder` and `-`.
- Changed `_coerce_probability` FDR threshold from `> 1` to `>= 1` (so a `1.0%` FDR → `0.01`), matching ProteoBench's `normalize()`. (Confirmed needed by the PEAKS 1% fixture, not speculative.)
- Symmetric-safety: both test sides (`from_series→to_series` and `extract_params→to_series`) pass through these validators, so model-level normalization only fixes raw↔canonical mismatches.

**S2–S6 — per-vendor parsers (10 vendors):** added each vendor's modification mapping producing ProForma-style strings (`C[Carbamidomethyl]`, `M[Oxidation]`, …) applied to `fixed_mods`/`variable_mods`; ported new regexes/defaults/enzyme-specific logic (e.g. DIA-NN `normalize_enzyme`, default `Trypsin/P`, `max_mods` default 0; FragPipe mass-token map + new `fragpipe_v23_noMBR` case; PEAKS `%`-stripping). No `fill_none()`/`normalize()` ported (pydantic model owns it); no `proteobench` imports. Sage's `restrict`-based enzyme branch kept (still needed for the `restrict="P"` → `Trypsin` case).

**Test changes:** stale enzyme assertions updated (`KR` → `Trypsin/P`) in `test_params_model.py` and `test_params_sage.py` (now that the model canonicalizes unconditionally); sage/maxquant equivalence tests aligned to the model round-trip pattern; FragPipe v23 case added.

**Fixtures:** 78 input+expected files copied from `ProteoBench/test/params/` into `tests/params/`; `PROTEOBENCH_PARAMS` repointed to `Path(__file__).resolve().parent / "params"` in all 11 referencing test files. Tests no longer depend on the sibling ProteoBench checkout. (`tools/generate_report.py` still reads param files from the ProteoBench checkout — it is a dev/demo tool, out of scope for test self-containment.)
