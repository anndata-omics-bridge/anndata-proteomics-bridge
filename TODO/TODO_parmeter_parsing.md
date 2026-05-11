# TODO: Parameter File Parsing Migration

## Goal

Move reusable proteomics parameter-file parsing from ProteoBench into APB
(`anndata_proteomics_bridge`) so APB becomes the upstream owner of vendor
parameter parsers and ProteoBench consumes APB instead of carrying duplicate
generic parsing logic.

This is a planning TODO only. Do not start implementation until the migration
plan is approved.

## Current State

ProteoBench currently owns the reusable parser implementations under:

```text
ProteoBench/proteobench/io/params/
```

Relevant parser modules include:

- `alphadia.py`
- `alphapept.py`
- `diann.py`
- `fragger.py`
- `i2masschroq.py`
- `maxquant.py`
- `msangel.py`
- `peaks.py`
- `proline.py`
- `quantms.py`
- `sage.py`
- `spectronaut.py`
- `wombat.py`
- de novo parsers: `adanovo.py`, `casanovo.py`, `deepnovo.py`,
  `instanovo.py`, `pihelixnovo.py`, `piprimenovo.py`, `pointnovo.py`

The parser contract is currently centered on `ProteoBenchParameters`, with
expected outputs tested by CSV fixtures in:

```text
ProteoBench/test/test_parse_params_*.py
ProteoBench/test/params/
```

ProteoBench also has a parameter overview matrix:

```text
ProteoBench/docs/parsing_overview.tsv
```

## Target Ownership

- APB owns reusable vendor parameter parsers.
- ProteoBench owns submission workflows, UI, benchmarking behavior, and
  ProteoBench-specific presentation of parsed parameters.
- ProteoBench should import APB parsers through a compatibility layer during the
  migration, then delete duplicated parser implementations after behavior is
  proven equivalent.

## Proposed APB Package Shape

Suggested APB module layout:

```text
src/anndata_proteomics/params/
  __init__.py
  model.py
  registry.py
  alphadia.py
  alphapept.py
  diann.py
  fragger.py
  i2masschroq.py
  maxquant.py
  msangel.py
  peaks.py
  proline.py
  quantms.py
  sage.py
  spectronaut.py
  wombat.py
```

De novo parsers should be decided explicitly: either migrate in the same
package under `params/denovo/`, or defer them until quant parameter parsing is
stable.

## Parameter Model

Define an APB-owned typed parameter model before moving parser logic.

Candidate fields are the current ProteoBench parameter fields, including:

- software name and version
- search engine and version
- enzyme
- missed cleavages
- fixed modifications
- variable modifications
- precursor and fragment mass tolerances
- peptide length bounds
- precursor charge bounds
- FDR settings
- MBR / reanalysis setting
- quantification method
- protein inference
- scan window

Open model decisions:

- Use Pydantic or dataclasses. Prefer Pydantic if strict validation and stable
  serialization matter; prefer dataclasses if the model should stay lightweight.
- Decide whether unknown or unsupported parameters are stored as `None`,
  omitted, or tracked in an explicit `unparsed` / `warnings` field.
- Preserve current ProteoBench CSV-compatible serialization during migration.
- Avoid adding public API beyond the parser entry points needed by ProteoBench
  and APB tests.

## Modification Normalization

Parameter parsing and modification normalization are coupled, but should not be
collapsed into one unstructured string format.

ProteoBench currently parses many vendor-specific representations of fixed and
variable modifications. APB should introduce a normalized internal modification
model before trying to unify vendor parser output.

Candidate internal representation:

- `name`: vendor or canonical modification name
- `accession`: controlled vocabulary accession when known, e.g. `UNIMOD:35`
- `target`: amino acid residue, terminus, or broader target
- `mod_type`: `fixed` or `variable`
- `position`: position/localization qualifier when available
- `mass_delta`: optional numeric mass shift
- `source`: original vendor string / file field
- `proforma`: optional ProForma rendering when the modification can be expressed
  safely

Important design point: ProForma is useful for modified peptide/proteoform
sequences. SDRF-Proteomics does not encode searched modifications as ProForma
strings; it encodes them in `comment[modification parameters]` as structured
key/value terms with `NT`, `AC`, `MT`, `TA`, and `PP` fields. Therefore APB
should be able to export SDRF-style modification parameter values separately
from ProForma peptide strings.

Reference examples from SDRF-Proteomics:

```text
NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=fixed;PP=Anywhere
NT=Oxidation;AC=UNIMOD:35;TA=M;MT=variable;PP=Anywhere
```

Relevant SDRF-Proteomics sources:

- Specification v1.1.0: https://sdrf.quantms.org/specification.html
- SDRF terms reference: https://sdrf.quantms.org/sdrf-terms.html
- SDRF GitHub repository: https://github.com/bigbio/proteomics-sample-metadata

Relevant ProForma sources:

- HUPO-PSI ProForma: https://www.psidev.info/proforma
- ProForma 2.0 paper/preprint: https://arxiv.org/abs/2109.11352
- CTDP ProForma resource: https://ctdp.org/resources/proforma/

## Migration Steps

1. Inventory every ProteoBench parser, fixture, expected CSV, and documented
   parameter in `docs/parsing_overview.tsv`.
2. Define the APB parameter model and modification model.
3. Port one low-risk parser first, preferably `sage.py`, with its fixtures and
   tests. Keep expected output byte-for-byte equivalent where practical.
4. Add a parser registry mapping software names and file signatures to parser
   functions.
5. Port the remaining quant parsers in small batches:
   - JSON/YAML-like formats first: Sage, AlphaDIA, AlphaPept, Wombat, quantms.
   - XML formats next: MaxQuant / MaxDIA.
   - Text/log/table formats next: DIA-NN, FragPipe, PEAKS, Spectronaut,
     MSAngel, ProlineStudio, i2MassChroQ.
6. Add APB CLI support only after the library API is stable, e.g.:

   ```bash
   anndata-proteomics parse-params <parameter-file> --software Sage
   ```

7. Add a ProteoBench compatibility layer that delegates to APB while preserving
   existing imports during transition.
8. Run APB parser tests and ProteoBench parser tests against the same fixtures.
9. Delete ProteoBench parser implementations only after compatibility tests pass.

## Test Strategy

- Port the existing ProteoBench expected CSV fixtures into APB tests or consume
  them from ProteoBench during the initial migration.
- Add unit tests for the normalized modification model:
  - vendor string to internal model
  - internal model to SDRF key/value encoding
  - internal model to ProForma where unambiguous
  - unsupported/ambiguous cases preserve the original source string
- Add one registry test per parser to verify software detection and parser
  dispatch.
- Add compatibility tests in ProteoBench that prove existing `extract_params`
  calls still return the same serialized values.

## Out of Scope For First Pass

- Rewriting ProteoBench UI.
- Changing submitted JSON schemas.
- Changing benchmark scoring.
- Rich parameter validation beyond current parsed fields.
- Converting every modification to ProForma if localization or target
  information is ambiguous.
- De novo parser migration, unless explicitly added to the approved plan.

## Open Questions

- Should APB expose parsed parameters as plain Python objects, Pydantic models,
  or a DataFrame/Series-compatible object?
- Should APB support auto-detection of parameter-file software, or require
  explicit software selection first?
- Should quant and de novo parameter parsers share one model, or should de novo
  get a narrower model?
- Which modification controlled vocabulary should be canonical internally:
  Unimod, PSI-MOD, or a neutral model that can carry either?
- Should SDRF export be part of the parameter parser API, or a separate exporter
  layered on top of parsed parameters?
