# apb parsing-rule coverage gaps — surfaced by an apb_studio corpus run

**Date:** 2026-07-02 · **Source:** `apb_studio` full-corpus `make run` over the ProteoBench
`json_dir` cache (2026-07-01). **Owner:** apb (parsing rules). **Status:** open — scoped, not started.

---

## TL;DR

`apb convert` fails on ~9 ProteoBench submissions with a **clean `ValueError` — "no rule covers
software version …"**. These are **not crashes**; apb is correctly reporting that its packaged parsing
rules don't cover those vendor **versions**. The work is to **add / widen parsing-rule variants** for
the newer FragPipe and PEAKS versions (and to investigate one PEAKS file whose version can't be
parsed). WOMBAT and the DIA-NN crash from the same run are already handled elsewhere (see
[§ Out of scope](#out-of-scope)).

The failure is raised in
[converters/pipeline.py](../src/anndata_proteomics/converters/pipeline.py#L68-L88) `select_rule` →
[rules/loader.py](../src/anndata_proteomics/rules/loader.py#L103-L126) `resolve_rule_for_version`.

---

## How rule↔version matching works (the thing to change)

- `apb convert` reads the **version** from `--params`
  ([`get_software_version` → `parse_params(...).software_version`](../src/anndata_proteomics/converters/pipeline.py#L63)).
- `select_rule(slug, level, version, headers)` calls `resolve_rule_for_version`, which picks the rule
  whose **`software_version` field (a regex)** matches the parsed version
  ([`_software_version_matches` = `re.search(rule.software_version, version)`](../src/anndata_proteomics/rules/loader.py#L121-L126)).
- Packaged rules pin that regex **narrowly**. E.g.
  [`parsing_rules/fragpipe/parse_fragpipe_ion_1.toml`](../src/anndata_proteomics/parsing_rules/fragpipe/parse_fragpipe_ion_1.toml)
  has `software_version = "^22\\.1-build02$"` — so `22.0` and `23.0` don't match → **no rule covers**.
- If a rule *is* selected but the file's columns don't fit, `select_rule` raises **"file columns don't
  match the rule …"** (no silent fallback).

So each gap is fixed by either **widening the `software_version` regex** of an existing rule (when the
columns are unchanged across that version) **or adding a new rule variant** (a `parse_<vendor>_<level>_2.toml`,
or a version subfolder like DIA-NN's `diann/v1`, `diann/v2`) when the columns changed. See
[docs/toml_schema.md](../docs/toml_schema.md) and the base/leaf convention in [AGENTS.md](../AGENTS.md).

---

## The gaps (from the 2026-07-01 run)

| Vendor / level | Version(s) reported | Count | Packaged rule today | Likely action |
|----------------|---------------------|-------|---------------------|---------------|
| **FragPipe / ion** | `22.0`, `23.0` | 3 | `parse_fragpipe_ion_1.toml` pinned `^22\.1-build02$` | verify columns for 22.0/23.0 vs the rule → widen the regex if unchanged, else add a `_2` variant |
| **PEAKS / ion** | `13 20250515`, `13 20250520` | 5 | `parse_peaks_ion_1.toml` | same: check PEAKS 13 columns vs the rule → widen / add variant |
| **PEAKS / ion** | `None` — "file columns don't match … version None" | 1 | — | version couldn't be parsed from `--params`; **investigate the param file** (parse gap) *and/or* the columns |

Representative inputs (paths relative to `apb/`):

```
# FragPipe 23.0
test_data_download/json_dir/Results_quant_ion_DDA_Astral/0c36dceb48ce85c9a56fb3d30c4f65ace4ed0aaf/input_file.tsv
# FragPipe 22.0
test_data_download/json_dir/Results_quant_ion_DDA/45486140efcbe205e2485f1ef4d668ec3d79fb99/input_file.txt
# PEAKS 13 (20250515/20250520)
test_data_download/json_dir/Results_quant_ion_DDA/9d1361331b165d6cc779ccf614419eb77057f573/input_file.csv
# PEAKS — version None / columns mismatch
test_data_download/json_dir/Results_quant_ion_DIA_AIF/b5fddd9b5d27918e8d31ec07bcf599cbd214027a/input_file.txt
```
The co-located `param_0..*` file in each dir is the `--params` (it supplies the version).

---

## Reproduce (one dataset)

```bash
cd apb && source .venv/bin/activate
D=test_data_download/json_dir/Results_quant_ion_DDA_Astral/0c36dceb48ce85c9a56fb3d30c4f65ace4ed0aaf
apb convert "$D/input_file.tsv" --software fragpipe --level ion \
    --params "$D/"param_0..* --output /tmp/probe.h5ad
# → ValueError: fragpipe ion: no rule covers software version '23.0'
```
To see the columns you're matching against: `head -1 "$D/input_file.tsv"` and compare with the
`[columns.*.select]` in the current `parse_fragpipe_ion_1.toml`.

---

## Suggested approach

1. For each (vendor, version): read the version's header columns and diff against the existing rule's
   `[columns.*.select]` / `[axis]`.
2. **Columns unchanged** → widen the rule's `software_version` regex to also match the new version
   (e.g. `^(22\.|23\.)` or a broader pattern). Keep it a regex; prefer an explicit set over `.*`.
3. **Columns changed** → add a new variant rule (`parse_<vendor>_<level>_2.toml`, or a version
   subfolder) with the new columns and a `software_version` regex for the new range; keep the old
   variant for the old version.
4. For the **PEAKS "version None"** case: check whether the param parser fails to extract a PEAKS
   version (a `params/parsers/peaks.py` gap) or whether it's simply a different column layout. Fix the
   parser and/or add the rule variant.
5. Add a test per new version under `tests/` (mirror `test_converters_e2e.py` / the packaged-rule
   tests), skipping when the ProteoBench data isn't present (existing pattern).

**Acceptance:** `apb convert <that input> --software <v> --level ion --params <p>` succeeds for each
listed version and writes the `<level>.h5ad`; `pytest tests/` green; `ruff` clean.

---

## Out of scope (already handled — context only)

- **DIA-NN param crash → degrade (fixed 2026-07-01, apb).** A FragPipe `.workflow` file mis-attached
  to a DIA-NN submission made the DIA-NN param parser raise `packaging.version.InvalidVersion` on an
  empty version. Fixed: `params/parsers/diann.py` tolerates a missing/garbage version and raises a
  clean `ParamsError` (new, `params/model.py`) for non-DIA-NN files; `converters/assemble.py`
  **degrades** — convert still produces the quant data, records
  `uns['anndata_proteomics']['search_parameters_error']`, and logs a warning. This is why some
  datasets now *convert with a params warning* instead of crashing.
- **WOMBAT `no rule covers version '0.9.11'` at `ion` → an apb_studio scaffold bug, NOT an apb gap.**
  apb intentionally supports WOMBAT only at **peptidoform** (and that path works). `make scaffold` was
  mislabeling WOMBAT datasets as `level: ion`; apb_studio now assigns the vendor-native level
  (WOMBAT → peptidoform). Re-running `make scaffold` there removes this failure — no apb change needed.
- **Run resilience / visibility (apb_studio).** The corpus runner now uses `--keep-going`, tees each
  rule to a per-dataset `<artifact>.log`, and surfaces "convert failed: <apb error>" per dataset — so
  these coverage gaps show as a clear per-dataset note rather than a wall of console errors.
