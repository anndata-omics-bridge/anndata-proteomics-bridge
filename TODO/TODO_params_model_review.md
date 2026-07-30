# Parameter Model Review

Date: 2026-07-02

This supersedes the misplaced root-level `TODO_params_model_review.md`.

Scope: review the two uncommitted comments in
`src/anndata_proteomics/params/model.py`:

- `COMMENT: map is for all the tools. Should this not be per tool?`
- `COMMENT : as few possible NONE fields.`

## Online Pydantic Skill Search

Pydantic `SKILL.md` / skill-adjacent URLs checked before the original review:

- <https://github.com/pydantic/skills>
  - Official Pydantic organization skills repository. Useful as the canonical
    first check, but currently focused on Pydantic AI / Logfire rather than core
    `BaseModel` schema review.
- <https://github.com/bobmatnyc/claude-mpm-skills/blob/main/toolchains/python/validation/pydantic/SKILL.md>
  - Useful unofficial Pydantic v2 validation skill covering `ConfigDict`,
    `field_validator`, `model_validator`, `model_dump`, strict validation, and
    model tests.
- <https://github.com/CJHarmath/claude-agents-skills/blob/main/skills/py-pydantic-patterns/SKILL.md>
  - Useful unofficial Pydantic v2 pattern note for validators, field
    constraints, and `extra="forbid"`.
- <https://github.com/microsoft/skills/blob/main/.github/plugins/azure-sdk-python/skills/pydantic-models-py/SKILL.md>
  - API-schema oriented Pydantic model skill. Less directly applicable to APB,
    but useful for thinking about separate required/optional model contracts.
- <https://github.com/pydantic/pydantic/issues/13189>
  - Official discussion about a future core Pydantic `SKILL.md`; useful reminder
    that unofficial skills are background, not normative guidance.

Official Pydantic docs used for semantics:

- <https://docs.pydantic.dev/latest/migration/#required-optional-and-nullable-fields>
- <https://docs.pydantic.dev/latest/concepts/validators/>
- <https://docs.pydantic.dev/latest/concepts/fields/>

## Corrected Findings

### High: required fields must not be represented as optional nullable fields

The `as few possible NONE fields` comment is not just a request to reduce JSON
noise. It is a schema-contract issue: fields that APB requires should be
required and non-empty, not silently accepted as `None` or empty strings.

The current model does not enforce that for many string fields:

- `software_name`, `software_version`, `search_engine`,
  `search_engine_version`, `quantification_method`, `protein_inference`,
  `abundance_normalization_ions`, and `predictors_library` are normalized by
  `_empty_strings_to_none`.
- `enzyme` also treats missing/empty values as `None`.
- `Parameters()` is currently valid, which makes every field effectively
  optional at model construction time.

Pydantic v2 distinction:

- `field: str` means required and non-null, but still needs a constraint or
  validator if empty / whitespace-only strings must be rejected.
- `field: str | None` without a default means required but nullable. That is
  usually not what we want for APB identity fields.
- `field: str | None = None` means optional and nullable.
- For required non-empty strings, use a constrained string type, for example
  `Annotated[str, StringConstraints(strip_whitespace=True, min_length=1)]`.

Recommended action:

- Classify `Parameters` fields into three groups:
  1. Required and non-empty.
  2. Optional, but if present non-empty / validated.
  3. Truly optional and nullable because not every vendor reports them.
- Convert the first group away from `| None = None`.
- Stop converting empty strings to `None` for required fields; empty strings
  should fail validation there.
- Update `from_series()` and parser tests so missing required fields fail
  explicitly instead of being hidden by placeholder normalization.

Likely required/non-empty candidates to review first:

- `software_name`
- `search_engine`
- `search_engine_version`
- `enzyme`

These are the fields most directly tied to dispatch/version selection and
downstream FASTA digestion. Some vendors may force a deliberate policy decision:
either parser-specific defaults are legitimate, or the parser should fail if the
parameter file cannot provide the required value.

`software_version` is a confirmed exception, not a required-field candidate. Two
supported PEAKS parameter files genuinely omit it; see
[`TODO_fragpipe_peaks_version_coverage.md`](Archive/TODO_fragpipe_peaks_version_coverage.md).
Keep it optional-but-observed: parameter JSON always writes the key with `null`,
while `search_parameters_version_status` and rule-selection provenance distinguish
missing data from parse errors.

### Medium: enzyme aliases need vendor provenance

The `_ENZYME_MAP` comment raises a valid boundary question. The corrected
position is:

- Vendor-specific raw syntax belongs in `params/parsers/<vendor>.py`.
- Shared APB canonical display names belong in the model. FASTA digestion does
  not need to know vendor aliases; it only needs the canonical
  `Parameters.enzyme` value and the cleavage rule for that canonical name.
- If a shared alias map remains in the model for now, it must record where each
  alias came from. A flat `raw_alias -> canonical_name` map loses the evidence
  that explains whether an alias is DIA-NN syntax, FragPipe/MSFragger syntax,
  ProteoBench expected-CSV compatibility, or a generic display-name spelling.

The data path is:

1. A vendor parser extracts an enzyme value from the parameter file.
2. `Parameters.enzyme` canonicalizes known aliases through `_ENZYME_MAP`.
3. `annotation/var_fasta.py` reads the stored search parameters and passes
   `params.enzyme` into `resolve_cleavage()`.
4. `fasta/annotation.py` maps canonical names such as `Trypsin`, `Trypsin/P`,
   and `Lys-C` to the actual cleavage rules used for peptide counting.

Evidence in the current repo:

- DIA-NN, Sage, FragPipe, and WOMBAT already do parser-local enzyme
  interpretation where raw syntax depends on vendor structure.
- `K*,R*` / `K*,R*,!P` are DIA-NN command-line syntax and are already handled in
  `params/parsers/diann.py`.
- `stricttrypsin` comes from FragPipe/MSFragger fixtures and is already handled
  in `params/parsers/fragpipe.py`.
- `src/anndata_proteomics/annotation/var_fasta.py` uses `params.enzyme` unless a
  cleavage override is supplied.
- `src/anndata_proteomics/fasta/annotation.py` only needs canonical enzyme
  names; it maps those names to cleavage rules.
- Tests assert ProteoBench parity for aliases such as `KR -> Trypsin/P`.

Recommended action:

- Short term: keep the shared alias map only if each alias carries provenance,
  for example alias, canonical enzyme, vendor(s), parser/test source, and a
  short note.
- Preferred end state: vendor parsers translate their own raw enzyme syntax to
  canonical APB enzyme names before constructing `Parameters`; the model then
  validates that the value is one of the supported canonical names.
- Add tests for the parser-local tokens that have vendor-specific meaning.
- Do not let an empty enzyme become `None` for records where enzyme is required.

### Low: the inline comments should not stay in source

The two `COMMENT` lines are review prompts, not durable source comments.

Recommended action:

- Remove both lines once the model contract change is implemented.
- Replace them with either code-level constraints or short explanatory
  comments/docstrings.

## Implementation Plan

1. Add a local `NonEmptyStr` type alias to `params/model.py`.
2. Mark the agreed required fields as `NonEmptyStr` rather than
   `str | None = None`.
3. Split string normalization so required fields reject empty strings and
   optional fields may still normalize configured missing sentinels to `None`.
4. Update `tests/test_params_model.py`:
   - `Parameters()` should fail once required fields are decided.
   - Empty required strings should fail.
   - Optional present-but-empty fields should follow the chosen policy.
5. Update vendor parser tests and fixtures for any parser that currently relies
   on model-level `None` defaults.
6. Re-run focused parameter and FASTA tests.
