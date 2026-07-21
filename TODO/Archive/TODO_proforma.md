# ARCHIVED: How does ProForma encode modifications in a peptidoform sequence?

Date: 2026-07-06

Status: answered; archived 2026-07-20. Residual modification work is tracked by the active
alphabase and registry-homogenization TODOs.

## Question

How does ProForma standardize modification encoding in a peptidoform sequence,
and what does APB need to implement or validate for vendor modified-sequence
columns?

## Short answer

Yes: ProForma standardizes how modifications are represented. It keeps the
amino-acid sequence readable from N-terminus to C-terminus, and modifications
are written as bracketed tags attached to the modified position.

The standardization has two parts:

1. placement syntax: where the bracket goes in the sequence;
2. modification identity: what goes inside the bracket.

ProForma does not maintain one APB-local list of all modification names. It
standardizes how to refer to modification terms from controlled vocabularies and
ontologies, plus how to encode mass shifts and richer annotations.

Examples:

| Meaning | ProForma |
|---|---|
| unmodified peptide | `PEPTIDE` |
| oxidized methionine by name | `PEPM[Oxidation]TIDE` |
| oxidized methionine by Unimod accession | `PEPM[UNIMOD:35]TIDE` |
| phosphorylated serine by accession | `PEPS[UNIMOD:21]TIDE` |
| N-terminal acetylation | `[Acetyl]-PEPTIDE` |
| N-terminal acetylation by accession | `[UNIMOD:1]-PEPTIDE` |
| C-terminal modification | `PEPTIDE-[Methyl]` |
| mass shift instead of named term | `PEPM[+15.9949]TIDE` |
| joint accession and name tag | `PEPS[UNIMOD:21|Phospho]TIDE` |
| precursor ion charge | `PEPM[UNIMOD:35]TIDE/2` |

The minimal model is:

```text
residue-localized:  AA[modification]
N-terminal:         [modification]-PEPTIDE
C-terminal:         PEPTIDE-[modification]
charged ion:        peptidoform/charge
```

The bracket content can be standardized in several ProForma-supported ways:

- controlled-vocabulary accession: `[UNIMOD:35]`, `[MOD:00046]`,
  `[RESID:AA0037]`;
- controlled-vocabulary name: `[Oxidation]`, `[Phospho]`;
- mass shift: `[+15.9949]`, `[-79.9663]`;
- formula or glycan tag: `[Formula:C12H20O2]`, `[Glycan:HexNAc1Hex2]`;
- multi-valued tag: `[UNIMOD:21|Phospho]`.

For APB, prefer registry-backed accessions such as `UNIMOD:35` where we can
resolve them, because accessions are more stable than names. Names are still
legal ProForma, but they are weaker as a data interchange target.

## Recommendation: normalize peptidoforms to accessions

For peptidoform sequences, APB should normalize to accession-backed ProForma
when the vendor token can be resolved.

Example:

```text
vendor token:     UniMod:35 / Oxidation / +15.9949 / oxM
APB peptidoform:  PEPM[UNIMOD:35]TIDE
```

This gives the feature identity a stable, machine-readable modification term.
The residue/terminal location is encoded by the ProForma placement itself. The
name, target, position, and mass delta come from APB's registry, not from each
vendor TOML.

If a token cannot be resolved, the existing `[modifications].unknown_policy`
decides what happens:

- `preserve`: keep the token in the rendered output;
- `drop`: remove the token;
- `error`: fail conversion.

The default should remain conservative: resolve known tokens to accessions and
make unknown-token behavior explicit in the parsing rule.

## Mapping tables in APB

APB has two mapping layers for result-table peptidoforms:

1. Vendor token to accession, in parsing-rule TOMLs:

   ```toml
   [[modifications.map]]
   token = "UniMod:35"
   accession = "UNIMOD:35"
   ```

   or for mass-delta vendors:

   ```toml
   [[modifications.map]]
   token = "15.9949"
   accession = "UNIMOD:35"
   ```

2. Accession to canonical metadata, in
   `src/anndata_proteomics/modifications/unimod_registry.toml`:

   ```toml
   [[entries]]
   accession = "UNIMOD:35"
   name = "Oxidation"
   target = ["M"]
   position = "Anywhere"
   mass_delta = 15.9949
   ```

Current registry entries are intentionally small and fixture-driven:

- `UNIMOD:1` Acetyl, N-term;
- `UNIMOD:4` Carbamidomethyl, C;
- `UNIMOD:21` Phospho, S/T/Y;
- `UNIMOD:27` Glu->pyro-Glu, N-term E;
- `UNIMOD:28` Gln->pyro-Glu, N-term Q;
- `UNIMOD:35` Oxidation, M;
- `UNIMOD:121` GG, K.

So yes, there is a mapping table, but it is not a complete Unimod dump. APB
adds mappings when supported vendor fixtures require them.

## Two different modification contexts

Do not collapse the peptidoform and parameter cases into one representation.
They answer different questions.

### 1. Peptidoform modifications

Question: what exact modified peptide feature was observed or quantified?

Representation:

```text
PEPM[UNIMOD:35]TIDE/2
```

Requirements:

- localized on the peptide sequence;
- compact enough to be a feature identifier;
- stable across vendors;
- suitable for `ProForma_peptidoform` and `ProForma_ion`;
- unknown-token policy is local to the result-table parser.

Use ProForma accession tags here.

### 2. Search parameter modifications

Question: what modifications were allowed in the database search?

Representation should be typed search metadata, not only a sequence string.
For SDRF this becomes `comment[modification parameters]`, for example:

```text
NT=Oxidation;AC=UNIMOD:35;MT=variable;TA=M;PP=Anywhere
NT=Carbamidomethyl;AC=UNIMOD:4;MT=fixed;TA=C;PP=Anywhere
```

Requirements:

- carries fixed vs variable (`MT`);
- carries target residue or terminus (`TA`);
- carries position/specificity (`PP`);
- carries accession (`AC`) and name (`NT`);
- is not localized to a particular observed peptide sequence;
- can be repeated once per searched modification in SDRF.

APB already has the right model shape for this:
`SearchedModification(name, accession, mod_type, target, position, mass_delta,
source)`. The current gap is that parameter parsers still mostly emit
ProteoBench-compatible ProForma-like strings such as `M[Oxidation]`, and the
model fallback stores those strings without resolving `accession`, `target`, or
`mass_delta`.

Recommendation for parameters:

- keep ProteoBench-compatible rendering where needed for legacy comparison;
- enrich the underlying `SearchedModification` with registry-backed accession,
  target, position, and mass delta;
- render SDRF from the typed object, not by reverse-parsing peptidoform strings;
- keep vendor aliases/code/mass mappings declarative where practical.

## What exactly is standardized?

ProForma standardizes the grammar and the allowed kinds of modification
descriptors. It supports public controlled vocabularies and ontologies such as:

- Unimod: `UNIMOD:35`, `UNIMOD:21`;
- PSI-MOD: `MOD:00046`;
- RESID: `RESID:AA0037`;
- XL-MOD for cross-linkers;
- GNO for glycans.

That means APB should not invent a private modification notation. APB should
translate vendor-specific tokens into standard ProForma descriptors.

Important distinction:

- `M[UNIMOD:35]` is a standard ProForma representation of oxidation on M.
- `M[Oxidation]` is also valid ProForma if the term resolves as intended.
- `M[+15.9949]` is valid ProForma for a mass shift, but less semantically rich.
- A vendor string such as `M(UniMod:35)`, `Oxidation (M)`, or `oxM` is not
  necessarily ProForma. APB must normalize it into ProForma.

So the APB task is not to decide a new encoding. The APB task is to map each
vendor token onto the existing ProForma encoding, preferably accession-backed.

## What this means for APB

APB's hard problem is not inventing the ProForma modification notation. The hard
problem is converting vendor-specific modified-sequence syntax into standard
ProForma tags.

For each supported vendor, APB needs to answer:

- Which source column contains the vendor modified sequence?
- Which syntax marks a residue-localized modification?
- Which syntax marks an N-terminal or C-terminal modification?
- Does the vendor report a name, accession, short code, numeric mass, formula,
  glycan, or custom token?
- Which vendor tokens map to `unimod_registry.toml` entries?
- What should happen when APB sees an unknown token: preserve, drop, or error?

The output should keep APB's current computed identifiers:

- `ProForma_peptide`: stripped peptide sequence.
- `ProForma_peptidoform`: peptide plus ProForma modification tags.
- `ProForma_ion`: `ProForma_peptidoform` plus `/charge`.
- `ProForma_fragment`: APB-specific fragment identifier,
  `{ProForma_ion}/{fragment_label}`. This is not plain ProForma.

## Scope APB should support first

Support and document the subset APB already needs for current vendor result
tables:

- unmodified peptide sequences;
- residue-localized modifications;
- N-terminal modifications;
- C-terminal modifications;
- multiple localized modifications on one peptidoform;
- positive integer precursor charges.

Keep these out of scope unless a real APB fixture requires them:

- labile modifications with `{...}`;
- global modifications with `<...>`;
- ambiguous or probabilistic localization;
- cross-links;
- chimeric peptidoforms;
- charge adduct expressions;
- full proteoform/top-down syntax;
- formulas, glycans, and arbitrary CV namespaces beyond APB's registry.

## Pyteomics recommendation

Use `pyteomics[proforma]` for optional tests and validation, not as APB's runtime
vendor parser.

Good uses:

- parse APB-generated `ProForma_peptidoform` values;
- parse APB-generated `ProForma_ion` values;
- catch malformed brackets, terminal syntax, and charge formatting;
- force semantic CV/mass resolution in focused tests by accessing parsed `.mass`.

Bad uses:

- do not feed vendor modified-sequence strings directly to Pyteomics before APB
  has normalized them;
- do not treat parse success as proof that a modification is biologically valid
  on a specific residue or terminus;
- do not validate `ProForma_fragment` with Pyteomics, because APB's fragment
  suffix is outside plain ProForma.

## mzSpecLib / mzspeclib-py recommendation

Use mzspeclib-py as prior art for vendor rewrite patterns and fixtures, not as a
runtime dependency.

Useful things to inspect:

- DIA-NN and Spectronaut modified peptide rewrites;
- MSP terminal-modification and alias handling;
- BiblioSpec / EncyclopeDIA modified sequence assumptions;
- controlled vocabulary names used for ProForma sequence and ion identifiers.

Do not import mzspeclib-py into APB unless there is a separate design decision:
it is a spectral-library implementation with a different data model.

## Work packages

1. Build a vendor fixture table with raw modified sequence, stripped peptide,
   expected `ProForma_peptidoform`, charge, expected `ProForma_ion`, and
   fragment label where applicable.
2. Extend `unimod_registry.toml` only from fixture evidence, and prefer
   accession-backed mappings where possible.
3. Add tests for residue, N-term, C-term, ion charge, and fragment identifier
   generation.
4. Add optional Pyteomics validation for the standard ProForma subset.
5. Document APB's supported subset and explicitly list unsupported ProForma
   features.

## Acceptance criteria

- The TODO question is answered in APB docs or tests with concrete examples.
- Current APB column names remain unchanged.
- Unknown-token behavior is covered for preserve/drop/error policies.
- Pyteomics tests, if added, validate only standard ProForma strings.
- `ProForma_fragment` is tested as APB-specific, not as a ProForma sequence.
- `uv run --extra dev pytest tests/` passes.
- `git diff --check` passes.

## References

- Pyteomics ProForma API:
  <https://pyteomics.readthedocs.io/en/latest/api/proforma.html>
- HUPO-PSI ProForma:
  <https://www.psidev.info/proforma>
- HUPO-PSI ProForma repository:
  <https://github.com/HUPO-PSI/ProForma>
- HUPO-PSI mzSpecLib:
  <https://www.psidev.info/mzspeclib>
- mzspeclib-py reference implementation:
  <https://github.com/HUPO-PSI/mzspeclib-py>
