# TODO: Multi-level quantification and MuData

## Status (2026-06-21): core question answered — feasibility proven on DIA-NN

The basic question ("each level as its own AnnData, MuData as a thin container") is
answered **yes** and demonstrated end-to-end for DIA-NN. Shipped:

- **Five DIA-NN level rules from one `report.tsv`**:
  `parse_diann_{ion,peptidoform,peptide,protein,fragment}_1.toml`.
  - protein uses **native** `PG.MaxLFQ`/`PG.Normalised`/`PG.Quantity`
    (`duplicates="keep_first"` — they are pre-aggregated and repeated per
    `(Run, Protein.Group)`).
  - peptidoform/peptide quant is an **APB-derived rollup** = summed
    `Precursor.Quantity` (`duplicates="aggregate"`), NOT MaxLFQ. Open question
    below: is a sum acceptable, or should a MaxLFQ-style rollup be implemented?
- **Fragment level** added as a fifth `QuantificationLevel` (decision: name is
  **`fragment`**, not `fragment_ion`). DIA-NN packs fragments as parallel
  `;`-delimited lists in each precursor row, so a new `[fragments]` TOML block +
  `converters/_fragments.explode_fragments` fan them out (pandas multi-column
  `explode`) before the normal pivot. New compute mode `proforma_fragment`
  (`ProForma_fragment = "{peptidoform}/{charge}/{fragment_label}"`).
- **MuData proof** (test-only, no public API per step 7): `tests/test_mudata_levels.py`
  builds all five AnnData, wraps them in `MuData(axis=0)`, and verifies the FK links
  and `.h5mu` round-trip. `tests/test_diann_levels.py` checks the per-level hierarchy.

Resolved design points (deviating from / refining this doc):
- **var_names must be prefixed per level** (`frg:/ion:/pfm:/pep:/prt:`). Empirically,
  axis=0 MuData tolerates colliding var_names but then *silently empties the merged
  `.var`*; peptide and peptidoform collide for unmodified peptides, so prefixing is
  mandatory, not optional.
- **FK link columns carry the PREFIXED parent id** (e.g. `ion.var["peptidoform_fk"] =
  "pfm:" + ProForma_peptidoform`). This refines the doc's "keep the bare identifier":
  a bare FK can't match a prefixed parent index. The bare id stays available as the
  non-key `.var` columns the rules already compute.
- **`recognize()` can't pick a level** for a multi-level vendor (all DIA-NN levels match
  the same headers) — it returns `None`; the level is selected explicitly via
  `load_packaged_rule(software, level)`.

Known limitation / follow-up:
- **Fragment level is memory-heavy at full scale**: explode multiplies rows ~12x; a full
  6-run DIA-NN report builds ~827k features. `convert_long` now scatters directly into the
  dense matrix (no `pivot_table`) and the fragment path trims unused columns before
  exploding, bringing the peak from ~13.5 GB down to ~6.5 GB. The matrix is ~90% dense, so
  this is largely irreducible without chunking — a chunked-streaming explode+scatter (never
  materialising the full ~5M-row exploded frame) is the next step if full-scale fragment is
  needed. Tests run on a row-capped subset; in practice, convert fragment per run / filtered.
- DIA-NN fragment columns vary by version (some exports lack `Fragment.Info` or carry a
  reduced `Fragment.Quant.*` set), so the fragment rule does not fit every DIA-NN file.

Still open (see sections below): protein↔peptide ambiguity / relation table, whether to
add the other vendors' levels, and the peptide/peptidoform rollup semantics.

## Question

The current project has TOML parse rules for one quantification level at a time,
currently mostly ion-level output. Many proteomics tools can report related
quantities at several levels:

- fragment ion
- ion / precursor
- peptidoform: peptide sequence plus modification state
- peptide: stripped peptide sequence
- protein / protein group

The question is not yet "how do we implement a full MuData exporter?" The first
question is more basic:

Can we represent each level as its own AnnData object, and can MuData be used
as a thin container around those AnnData objects?

## Current project fit

The existing design already points in this direction:

- one TOML rule maps one vendor table shape to one AnnData object
- the rule filename already includes the level:
  `parse_<software>_<level>_<file_version>.toml`
- `ParseRule.quantification_level` already supports `ion`, `peptidoform`,
  `peptide`, and `protein`
- `registry.find_rule(software, quantification_level, file_version)` already
  resolves a level-specific rule

So the pragmatic next step is not one large multi-level TOML. It is one TOML per
level:

```text
parsing_rules/diann/
  parse_diann_ion_1.toml
  parse_diann_peptidoform_1.toml
  parse_diann_peptide_1.toml
  parse_diann_protein_1.toml
```

Add a fragment-level rule only after checking the actual DIA-NN fragment output
columns and deciding whether the level should be called `fragment` or
`fragment_ion`. Unlike the other levels, fragment is not a pure-TOML addition:
it requires a schema change to extend the `QuantificationLevel` literal in
`rules/schema.py`, and a fragment key would need a new compute mode
(`proforma_ion` is hard-coded to ion-level rules only).

## AnnData level model

Each level should remain a normal AnnData:

- `.obs`: samples/runs
- `.var`: features for exactly one level
- `.X`: primary quantitative matrix
- `.layers`: alternative measurements for the same sample x feature matrix

This matters because "multiple levels" are not AnnData layers. Layers are only
valid when the feature axis is identical. For example, `Precursor_Normalised`,
`Precursor_Quantity`, `Q_Value`, and `RT` can be layers for the same ion-level
matrix. Protein, peptide, peptidoform, ion, and fragment matrices have
different `.var` axes, so they should be separate AnnData objects.

## Linking levels

Keep links pragmatic: use columns in each level's `.var` table. The link column
names are not a convention to maintain by hand — the rule schema enforces them.
`[[columns.var.compute]]` binds `how = "stripped_sequence"` to the name
`ProForma_peptide`, `how = "proforma_sequence"` to `ProForma_peptidoform`, and
`how = "proforma_ion"` to `ProForma_ion` (`_PROFORMA_COMPUTE_NAME` in
`rules/schema.py`). Any other name is rejected at load time, so the same
identifier means the same thing at every level for free.

Example:

- ion-level `.var` primary key (its `var_names`): `ProForma_ion`
- ion-level `.var` also carries the column `ProForma_peptidoform`
- peptidoform-level `.var` primary key (its `var_names`): `ProForma_peptidoform`

Then `ion.var["ProForma_peptidoform"]` points into `peptidoform.var_names`. Note
the asymmetry: the same identifier is a **column** at the child level and the
**index** (`var_names`) at the parent level. The link is a foreign key from a
child column into the parent index, and it works only when the parent index is
unique.

This gives a natural many-to-one relationship:

```text
many ions -> one peptidoform
many peptidoforms -> one peptide
many peptides -> one or more proteins
```

The important requirement is that the target column is unique in the target
level. For example:

- `peptidoform.var["ProForma_peptidoform"]` should be unique
- `peptide.var["ProForma_peptide"]` should be unique
- protein identifiers need more thought because protein groups and shared
  peptides can create one-to-many or many-to-many relationships

Do not create a separate link table in `uns` yet. The first storage mechanism
should be the `.var` columns that the TOML rules already define.

## Protein ambiguity

The main hard case is not ion -> peptidoform or peptidoform -> peptide. Those
are usually straightforward if identifiers are normalized.

The hard case is peptide -> protein:

- a proteotypic peptide maps to one protein or protein group
- a shared peptide can map to several proteins or groups
- vendor output may encode this as a semicolon-separated list

This breaks the foreign-key model that works for the lower levels. A
semicolon-separated list is not a scalar key and the target is not unique, so
`peptide -> protein` cannot be a single `.var` column pointing into a unique
parent index the way `ion -> peptidoform -> peptide` can. Protein linking is a
separate mechanism, and it is the case most likely to require the relation table
deferred below. Keeping the vendor columns is still useful — but as feature
metadata, not as a usable link.

For the first iteration, keep the vendor-derived protein mapping as `.var`
columns such as:

```text
Protein_Group
Protein_Ids
Protein_Names
Genes
```

Then decide later whether these columns need normalization into a relation
table. The relation table should be a derived view, not the first storage
mechanism.

## MuData interpretation

MuData can hold one AnnData per level in `.mod`:

```python
from mudata import MuData

mdata = MuData(
    {
        "ion": ion_adata,
        "peptidoform": peptidoform_adata,
        "peptide": peptide_adata,
        "protein": protein_adata,
    },
    axis=0,
)
```

With `axis=0`, observations are shared. For this project that means samples or
runs are the shared axis. This fits if all level-specific AnnData objects use
the same `.obs_names`.

The `.var` issue:

- each modality keeps its own level-specific `.var`
- `mdata.mod["ion"].var` is the ion feature table
- `mdata.mod["peptidoform"].var` is the peptidoform feature table
- MuData also has a global `.var`, but for this use case the important feature
  metadata should stay in the modality-specific `.var` tables

So different `.var` tables per level are expected and are not a blocker. The one
thing that must be handled is global `.var_names`: MuData concatenates variable
names across modalities, and the peptide and peptidoform levels will collide.
For an unmodified peptide the stripped sequence equals the ProForma sequence, so
`ProForma_peptide` and `ProForma_peptidoform` hold the same string and produce
duplicate global `var_names`. Disambiguating `var_names` per level (for example
a short prefix such as `pep:` / `pfm:`) is therefore required, not optional. The
link columns inside `.var` keep the bare identifier so the foreign-key join
above still works.

### Ignore the global `mdata.var`

MuData maintains a global `mdata.var`, but it is not a usable cross-level
registry and must not be treated as one.

- It is an auto-generated projection of the modality `.var` tables, rebuilt from
  them on every write. It cannot be turned off: `set_options` controls only the
  display style, and `update()` takes no arguments. Clearing its columns in
  memory holds, but the `.h5mu` writer regenerates them on the next save.
- Its rows are every feature across all levels, so `len(mdata.var) == mdata.n_vars`
  is the sum of all per-level feature counts, never just one level.
- Columns shared by *every* modality are pulled up unprefixed; columns present
  in only some modalities are modality-prefixed and `NA`-filled elsewhere (for
  example `ion:ProForma_peptidoform`). The clean column set one might hope for
  does not materialise, and the link pointers in particular never align: a
  parent key is an index at its own level and a column only at child levels.

The rule is simple: never read or write links through `mdata.var`. The
authoritative feature metadata and the links both live in `mdata.mod[level].var`,
and that is the only thing to read. The global frame is a harmless byproduct.

If a single flat cross-level registry is ever wanted, it is a derived view
computed on demand from the modality `.var` tables — the same `uns` relation
table deferred below — not the MuData-managed global `.var`.

## QFeatures prior art

QFeatures is useful as a conceptual reference because it was designed for
quantitative mass-spectrometry data across multiple assay levels.

In QFeatures:

- each assay is a quantitative matrix for one level
- assays can represent PSMs, precursors, peptides, proteins, or protein groups
- `AssayLinks` records relationships between assays
- an `AssayLinks` object contains individual `AssayLink` instances
- these links are used to describe the hierarchy between assays, for example
  peptide-level data derived from lower-level features

The useful idea for this project is the hierarchy, not necessarily the storage
implementation. For AnnData/MuData, the stored relationship is the TOML-defined
`.var` columns in each modality. If an AssayLinks-style cross-level registry is
ever needed, it is a derived view computed from those columns, not the
MuData-managed global `.var`.

## Proposed near-term direction

1. Keep one TOML per software and quantification level.
2. For each level, define the primary feature identifier column in `.var`.
3. Use TOML-defined `.var` column names to point from lower levels to higher
   levels, e.g. `ProForma_ion`, `ProForma_peptidoform`, `ProForma_peptide`,
   `Protein_Group`, `Protein_Ids`, `Protein_Names`, `Genes`.
4. Check uniqueness of target identifiers within each target level.
5. Disambiguate `var_names` per level so peptide and peptidoform features do not
   collide in the global variable axis.
6. Treat MuData as an optional thin wrapper around several AnnData objects.
   Ignore its auto-generated global `.var`; always read `mdata.mod[level].var`.
7. Do not add an `uns` link table, family manifest, CLI, or new public API until
   the level concepts are clearer on real DIA-NN examples.

## Concrete checks for DIA-NN

Before implementation, inspect real DIA-NN outputs and answer:

- Which columns define the ion/precursor feature?
- Which columns define the peptidoform feature?
- Which column should be the peptide key?
- Which protein/protein-group columns are available?
- Are the proposed target keys unique at their level?
- How are shared peptides encoded?
- Are protein-level quantities in the same report or a separate output?

## Open questions

- Which TOML-defined columns should be present for each level-specific `.var`
  table?
- Is `peptidoform` enough, or do modification-site-level quantities need their
  own level?
- Should `fragment` mean fragment ion, or should the level be named
  `fragment_ion`?
- How much protein ambiguity should be normalized during conversion, and how
  much should remain as vendor metadata already captured by the TOML-defined
  `.var` columns?

## References

- MuData documentation: https://mudata.readthedocs.io/stable/io/mudata.html
- MuData quickstart: https://mudata.readthedocs.io/stable/notebooks/quickstart_mudata.html
- QFeatures overview: https://rformassspectrometry.github.io/QFeatures/
- QFeatures AssayLinks: https://rdrr.io/bioc/QFeatures/man/AssayLinks.html
