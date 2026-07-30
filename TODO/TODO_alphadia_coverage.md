# AlphaDIA ion coverage across the DIA modules

> Convert the 11 cached AlphaDIA submissions. Unlike MaxQuant/Sage/AlphaPept this needs three new
> pieces, not one rule document: a modification parser for alphabase-style parallel
> `mods`/`mod_sites` columns, an AlphaDIA parameter-log parser, and three shape-distinct rule
> documents.

**Date:** 2026-07-30 · **Source:** header and content survey of the cached AlphaDIA fixtures.
**Owner:** apb (modifications, params, parsing rules). **Status:** planned — awaiting approval on
the modification-parser question below.

---

## Verified ground truth

11 cached submissions, 24 catalogued (13 more downloadable), across four DIA modules and **six
versions in three distinct shapes**:

| Shape | Versions | Format | Cols | Layout |
|---|---|---|---:|---|
| A | 1.10.3, 1.10.4-dev0 | `.tsv` | 13 | **wide** — one bare run column per sample |
| B | 1.12.1, 1.12.2 | `.tsv` | 45 / 92 | **long** — `run` + `intensity` |
| C | 2.1.0, 2.1.1 | `.parquet` | 33 | **long** — dotted namespaces, `raw.name` + `precursor.intensity` |

Per module: dia_astral 5 (1.10.3, 1.10.4-dev0, 1.12.1, 2.1.0, 2.1.1), dia_diapasef 3 (1.10.3,
1.10.4-dev0, 1.12.1), dia_zenotof 2 (1.12.1, 1.12.2), dia_singlecell 1 (2.1.0).

### The blocker: modifications are parallel columns, not inline tokens

AlphaDIA carries alphabase-style modifications in **two parallel columns** beside a bare sequence:

```
sequence     mods                            mod_sites
TCSSFIAAMER  Oxidation@M;Carbamidomethyl@C   9;2
TCSSFIAAMER  Carbamidomethyl@C               2
RPIAHLPCPGK  Carbamidomethyl@C               8
```

`mod_sites` are 1-based residue positions, paired index-wise with `mods`, and **not sorted** —
`9;2` pairs `Oxidation@M`→9 and `Carbamidomethyl@C`→2. Observed tokens across the corpus:
`Carbamidomethyl@C`, `Oxidation@M`, `Acetyl@Protein_N-term`, and repeats of the same token
(`Carbamidomethyl@C;Carbamidomethyl@C` for two cysteines).

APB's only modification runtime is `token_regex`, which needs one `source_column` with inline
tokens. It cannot read a name list against a position list.

**This is not optional.** Rows 1 and 2 above are the *same* sequence and charge but *different*
peptidoforms. Without a parser, both compute `ProForma_ion = TCSSFIAAMER/2` and collapse into one
feature — silently summing an oxidised and a non-oxidised precursor. Dropping modifications here
corrupts the matrix rather than merely losing an annotation, so no partial-credit shortcut exists.

`mod_seq_charge_hash` is present and *is* a correct unique key, but an opaque integer hash is
useless as `var_names` for ProteoBench scoring or any downstream report, and it breaks APB's
`ProForma_ion` convention.

### No parameter parser exists

`param_0..txt` is an AlphaDIA **run log**, not a config file: ANSI-coloured, timestamped, with an
indented config tree and `[user defined, default: X]` annotations. Everything needed is greppable
after stripping ANSI codes:

```
version: 1.10.3
├──enzyme: trypsin/p [user defined, default: trypsin]
├──fixed_modifications: Carbamidomethyl@C
├──variable_modifications: Oxidation@M;Acetyl@Protein_N-term
├──max_var_mod_num: 1 [user defined, default: 2]
├──missed_cleavages: 1
├──target_ms1_tolerance: 10 [user defined, default: 5]
├──target_ms2_tolerance: 15 [user defined, default: 10]
├──fdr: 0.01
```

`params/registry.py` has no `alphadia` entry, so capability discovery stops at "no parameter parser"
before any rule is consulted.

### Shape A needs a negative-lookahead sample regex

The wide files' sample columns are bare run names with no suffix to anchor on, interleaved with
per-feature columns — the exact failure mode recorded for PEAKS in
[TODO_fragpipe_peaks_version_coverage.md](Archive/TODO_fragpipe_peaks_version_coverage.md). **Verified** on
all four wide fixtures: a regex excluding the known metadata names splits them correctly, 6 samples
and 7 metadata columns every time:

```
^(?P<sample>(?!(?:mod_seq_charge_hash|genes|decoy|mods|mod_sites|sequence|charge|proteins|pg)$).+)$
```

The diaPASEF wide runs carry a trailing acquisition id (`..._Sample_Alpha_01_11494`) that the module
annotation's `raw_file` lacks, so shape A also needs `sample_name_cleanup` with `_\d+$`. The Astral
run names are unaffected by that pattern.

Shape A also repeats each feature (98 570 rows for 26 166 unique `mod_seq_charge_hash`), so it needs
an explicit `duplicates` policy rather than `error`.

---

## Plan

### Phase 1 — `site_list` modification parser

`rules/schema.py` gains a second modification parser beside `TokenRegexModifications`, selected by
the existing `parser` field as a discriminated union:

```json
{
  "parser": "site_list",
  "sequence_column": "sequence",
  "modification_column": "mods",
  "site_column": "mod_sites",
  "delimiter": ";",
  "site_base": 1,
  "map": [{"token": "Carbamidomethyl@C", "accession": "UNIMOD:4"}]
}
```

`modifications/apply_rules.py` gains a matching runtime that pairs the two lists index-wise, sorts
by site, resolves tokens through the same Unimod registry, and renders ProForma through the existing
`render_proforma`. `modifications/pipeline.apply_modifications` dispatches on `parser`. The public
output contract is unchanged: `proforma_sequence`, `stripped_sequence`, `unknown_mod_tokens`.

**Deliberately narrow.** This adds a parser *kind*; it does not touch the vendor-token → accession
mapping tables. Adopting alphabase's curated `modification_mappings` stays in
[TODO_modification_alphabase_2026-07-07.md](TODO_modification_alphabase_2026-07-07.md), whose stated
plan is one unified registry-backed design. Nothing here forecloses that.

Registry additions: `Acetyl@Protein_N-term` maps to `UNIMOD:1`, whose registry `position` is
`N-term` — confirm the site value AlphaDIA writes for a protein N-terminal mod (expected `0`).

### Phase 2 — AlphaDIA parameter parser

`params/parsers/alphadia.py` + a `params/registry.py` entry + a committed fixture in `tests/params/`.
Strip ANSI, strip the `[user defined, default: …]` suffix, read the fields above into `Parameters`
(`software_version`, `enzyme`, `allowed_miscleavages`, `fixed_mods`, `variable_mods`, `max_mods`,
`ident_fdr_psm`, `precursor_mass_tolerance`, `fragment_mass_tolerance` — both ppm).

### Phase 3 — three rule documents

`parsing_rules/alphadia/{v1_10,v1_12,v2}/rules.json`, following the existing `diann/v1`, `diann/v2`
precedent. All three declare the **ion** level only:

- **v1_10** — wide, `^1\.10\.`, the lookahead sample regex, `sample_name_cleanup` `_\d+$`.
- **v1_12** — long, `^1\.12\.`, `run` + `intensity`. 1.12.1 (45 cols) and 1.12.2 (92 cols) are not
  subsets of each other, so the volatile extras use the `optional_select` and optional-layer
  mechanisms already in place.
- **v2** — long, `^2\.`, dotted source names sanitised to underscores per
  [conventions.md](../../anndata_omics_bridge/docs/conventions.md); obs key `raw.name`,
  x_layer `precursor.intensity`.

### Out of scope

- **`peptide.intensity` / `pg.intensity` in shape C.** The 2.1.x table carries precursor, peptide
  and protein-group intensities side by side, so peptide and protein levels are technically
  reachable. Both are rollups of the precursor quantity, and
  [apb-toml-level-design](../skills) says standalone levels follow real quantitative layers rather
  than derived rollups — `pg.intensity` is directLFQ output and arguably qualifies. Worth a separate
  decision; not needed for the ProteoBench ion modules.
- The 13 catalogued-but-not-downloaded submissions.
- alphabase `modification_mappings` adoption (its own TODO, above).

---

## Verification

1. All 11 fixtures report `supported` with `('mudata', 'ion')`.
2. `obs_names` equal each module annotation's `raw_file` values — including the diaPASEF wide
   fixtures, which exercise `sample_name_cleanup`.
3. A peptidoform round-trip: `TCSSFIAAMER` + `Oxidation@M;Carbamidomethyl@C` + `9;2` renders
   `TC[UNIMOD:4]SSFIAAM[UNIMOD:35]ER`, and the oxidised and non-oxidised forms stay distinct
   features.
4. No regression: the other 76 supported fixtures unchanged, 0 failed.
5. `uv run pre-commit run --hook-stage pre-push --all-files` in both repos; `CHANGES.md` entries.
