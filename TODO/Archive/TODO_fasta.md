# TODO: FASTA-derived annotations → protein `var` table

Extend the `var` table of the **protein-level** AnnData with annotations derived
from FASTA file(s), the way prolfquapp's `get_annot_from_fasta()` /
`build_protein_annot()` annotate the protein row table. Merged **only** into the
protein layer.

## Decisions (locked in with the user)

- **Invocation:** a new **`apb fasta`** subcommand + a public function, both
  taking an AnnData/MuData and one-or-more FASTA paths. **No annotation-TOML
  schema changes.** (Subcommand named `fasta` to match the `apb convert /
  annotate / fasta` family; this is the old draft's `annotate-var`.)
- **Protein-group join:** leading accession only (prolfquapp-faithful — mirrors
  `tidyr::separate(protein_Id, c("cleanID", NA))`).
- **Cleavage is NOT assumed tryptic.** We already parse the digestion enzyme
  from the vendor parameter file; the peptide count must use the *actual*
  cleavage rule, not a hardcoded `[KR](?!P)`. (99% of the time it is trypsin —
  but Lys-C / Glu-C / etc. happen, so we read it rather than assume it.) See §0.
- **Reuse, don't re-implement:** build on the existing
  `fasta.annotation.fasta_to_dataframe()` (already a faithful port of
  `get_annot_from_fasta()`) and `annotation._sanitize.sanitize_columns`.

## What already exists

- `src/anndata_proteomics/fasta/annotation.py` — `fasta_to_dataframe()` produces
  `fasta.id`, `fasta.header`, `proteinname`, `gene_name` (gated on >1 match),
  `protein_length`, `nr_tryptic_peptides`, optional `sequence`.
  ⚠️ `count_tryptic_peptides()` **hardcodes trypsin** (`[KR](?!P|$)`, min 7 /
  max 30). This is the thing §0 generalizes.
- `src/anndata_proteomics/annotation/` — the merge module, currently
  `annotate_obs` only. Its schema docstring explicitly reserves a sibling `[var]`
  slot for feature annotations.
- **Cleavage params are already on the object.** `converters/assemble.py` calls
  `params.anndata_io.write_search_parameters`, storing the parsed `Parameters`
  under `uns['anndata_proteomics']['search_parameters']` (when a params file was
  given at convert time). `params.anndata_io.read_search_parameters(adata)`
  returns a `Parameters | None`. The model carries `enzyme` (canonicalised:
  `Trypsin`, `Trypsin/P`, `Lys-C`, `Arg-C`, `Asp-N`, `Glu-C`, `Chymotrypsin`, …),
  `allowed_miscleavages`, `min_peptide_length`, `max_peptide_length`.
- Protein level: `parse_diann_protein.toml` / `parse_spectronaut_protein.toml`
  → `var_names = Protein.Group`; var columns `Protein_Group / Protein_Ids /
  Protein_Names / Genes`. In a MuData this is the `protein` modality (`prt:`
  prefix on `var_names`).
- `uns['anndata_proteomics']['quantification_level']` carries the level (set in
  `converters/assemble.py`).

## 0. Generalize peptide counting to the real cleavage rule

`fasta/annotation.py`:

- Add an **enzyme → cleavage-rule** map keyed by the canonical names the params
  model already emits (reuse those names so the two cannot drift). Each entry is
  a regex of cut positions; start with the common set, default to trypsin:
  - `Trypsin` → `[KR](?!P)` · `Trypsin/P` → `[KR]` · `Lys-C` → `K(?!P)` ·
    `Arg-C` → `R(?!P)` · `Glu-C` → `[DE](?!P)` · `Chymotrypsin` → `[FYW](?!P)` ·
    `Asp-N` → cut *before* `D` (N-terminal) — handle the before-cut case.
  - Unknown / `None` enzyme → **warn once, fall back to Trypsin**.
- Generalize `count_tryptic_peptides(seq, *, cleavage=…, min_length, max_length)`
  to take the cleavage rule (keep a `cleavage=Trypsin` default so existing
  callers/tests are unchanged). `fasta_to_dataframe(...)` gains a `cleavage`
  (and optional `min_length`/`max_length`) pass-through.
- **Column name:** keep `nr_tryptic_peptides` for prolfquapp compatibility, but
  it now reflects the configured enzyme; the enzyme used is recorded in
  provenance (§1.5). *(Open: rename to `nr_peptides` if downstream doesn't need
  the prolfquapp name — flag for the user.)*

## 1. New function — `annotation/var_fasta.py`

```python
def annotate_var_from_fasta(
    obj,                          # protein-level AnnData, or MuData (targets mod["protein"])
    fasta_sources,                # path | list[path] | FastaSource(s)
    *,
    match_on="Protein_Group",     # var column carrying the group; "index" => var_names
    is_uniprot=True,
    decoy_pattern="^REV_|^rev_",
    cleavage=None,                # None => read enzyme from uns search_parameters; else override
    min_length=None,              # None => params.min_peptide_length, else default 7
    max_length=None,              # None => params.max_peptide_length, else default 30
    include_sequence=False,
    columns=None,                 # optional subset of derived columns to merge
) -> obj
```

Behaviour, mirroring `annotate_obs` semantics for consistency:

1. **Resolve target.** MuData → `mod["protein"]` (raise if absent). AnnData →
   require `uns['anndata_proteomics']['quantification_level'] == "protein"`, else
   raise. Enforces "protein layer only."
2. **Resolve cleavage rule.** If `cleavage`/`min_length`/`max_length` are not
   passed, call `read_search_parameters(target)` and take `enzyme`,
   `min_peptide_length`, `max_peptide_length` from it; if no params stored,
   warn and default to Trypsin / 7 / 30. (No trypsin assumption when we have the
   real enzyme.)
3. **Build FASTA frame** via `fasta_to_dataframe(..., cleavage=…, min_length=…,
   max_length=…)`, indexed by `proteinname`.
4. **Compute var-side join key** from `match_on` (default `Protein_Group`, or
   `var_names` if `"index"`): strip a `prt:` modality prefix → split on `;`, take
   first token → if `db|ACC|NAME`, take the middle accession. Yields the same
   form as the FASTA `proteinname` (matches prolfquapp's cleanID).
5. **Left-join** the derived columns (`fasta.id`, `fasta.header`,
   `protein_length`, `nr_tryptic_peptides`, `gene_name` when present, `sequence`
   if requested) onto `var`, aligned by key. Sanitize added names; **raise** on
   overlap with existing var columns, **raise** on zero matches, **warn** on
   partial match. (`fasta.id` / `fasta.header` keep their dots — sanitisation
   preserves dots.)
6. **Provenance:** append a `var_annotations_json` entry under the protein
   AnnData's `uns['anndata_proteomics']` (parallel to obs's
   `obs_annotations_json`), recording the FASTA source(s) **and the enzyme /
   min / max used** for the peptide count.

## 2. New CLI subcommand — `fasta` in `scripts/cli.py`

```
apb fasta <data.h5ad|.h5mu> <fasta...> \
    [--output] [--match-on] [--no-is-uniprot] [--decoy-pattern] \
    [--cleavage ENZYME] [--min-length N] [--max-length N]
```

Loads via `load_converted_result`, calls the function, writes back
(`.h5ad` / `.h5mu`); `--output` defaults to `<stem>.annotated<suffix>`
(non-destructive) — same pattern as the existing `annotate` command. The
`--cleavage/--min-length/--max-length` flags **override** the enzyme read from
`uns` (for objects converted without a params file).

## 3. Tests — `tests/test_annotation_var_fasta.py`

Reuse the existing `PROLFQUAPP_FIXTURE`:

- protein AnnData merge (values + leading-accession join across `;`-groups,
  `sp|..|..` ids, `prt:` prefix);
- **cleavage:** enzyme read from `uns` drives the count (Trypsin vs Trypsin/P vs
  Lys-C give different `nr_tryptic_peptides`); no-params object falls back to
  trypsin with a warning; `--cleavage` override wins;
- MuData targets `mod["protein"]` only and round-trips through `.h5mu`;
- non-protein AnnData (e.g. ion-level) raises;
- zero-match raises / partial-match warns;
- re-run column collision raises;
- provenance (incl. enzyme used) recorded in `uns`;
- one CLI smoke test.

## 4. Docs

- `docs/ARCHITECTURE.md` — annotation module now covers obs + var-from-fasta;
  note the cleavage rule comes from the stored search parameters.
- One-line update to the annotation schema docstring's "var reserved" comment
  pointing at the new function.

## Files touched

- `src/anndata_proteomics/fasta/annotation.py` (generalize cleavage; enzyme map)
- `src/anndata_proteomics/annotation/var_fasta.py` (new)
- `src/anndata_proteomics/scripts/cli.py`
- `tests/test_annotation_var_fasta.py` (new)
- `docs/ARCHITECTURE.md`

## Verify

```bash
pytest tests/test_annotation_var_fasta.py
```
