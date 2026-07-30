# TODO: Resolve vendor-to-module sample-name mismatches

Status: superseded on 2026-07-29.

The accepted pipeline contract is now:

```text
apb convert -> apb annotate -> apb proteobench
```

Annotation alone resolves vendor identifiers (including exact
`raw_file_alias` fallback) and writes `obs["sample_name"]` and
`obs["condition"]`. ProteoBench scoring requires those columns and does not
contain a second run-name resolver. The scoring-side resolver proposed below
must not be implemented; the remaining historical analysis is retained only
for context.

## Outcome from the corpus run

Corpus run `6bad57331ccb4aeb9f4b38baf8772a89` completed 527 of 541
Snakemake targets. All 147 `convert` and all 147 `fasta` targets completed.
The 14 failures were:

- 8 `annotate` targets;
- 6 `proteobench` targets;
- 0 conversion targets.

The failures are 14 stage/branch manifestations of three sample-alignment
problems:

| Fixture | APB `obs_names` | Module `[[samples]].raw_file` | Missing mapping |
| --- | --- | --- | --- |
| WOMBAT ion `58dfec05` | `A_1` … `B_3` | `LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01` … | `abundance_A_1` → `Condition_A_Sample_Alpha_01` |
| WOMBAT peptidoform `557e8c64` | `A_1` … `B_3` | `abundance_A_1` … `abundance_B_3` | `abundance_A_1` → `Condition_A_Sample_Alpha_01` |
| PEAKS diaPASEF `5691d485`, `806987cb` | `LFQ_ttSCP_diaPASEF_Condition_A_Sample_Alpha_01` … | `ttSCP_diaPASEF_Condition_A_Sample_Alpha_01_11494` … | `LFQ_ttSCP_…_01` → `Condition_A_Sample_Alpha_01` |

In every case, direct matching against `raw_file` finds zero of six runs. The
upstream mapping resolves all six runs uniquely and completely.

## Root cause

The converted observation names are correct representations of the vendor
exports:

- WOMBAT columns are `abundance_A_1`, … and the wide rule extracts the sample
  token `A_1`, …;
- PEAKS columns are
  `LFQ_ttSCP_diaPASEF_Condition_A_Sample_Alpha_01 Normalized Area`, … and the
  wide rule extracts the `LFQ_…` prefix as the sample token.

The ProteoBench module TOMLs use a different namespace. They describe either
raw acquisition names with module-specific prefixes/numeric suffixes or, for
the peptidoform module, the full WOMBAT abundance column name.

ProteoBench bridges these namespaces with per-module, per-tool `[run_mapper]`
tables. The mappings are already present at APB's pinned settings revision
`2738c47f8d621f0ee1fa4a6d3d358846f2bfa261`, for example:

- `DDA/ion/QExactive/parse_settings_wombat.toml`;
- `DDA/peptidoform/parse_settings_wombat.toml`;
- `DIA/ion/diaPASEF/parse_settings_peaks.toml`.

However, `apb-testdata annotations` downloads only each module's
`module_settings.toml`. It does not acquire the tool-specific parse settings,
so the only authoritative crosswalk is discarded.

Two downstream implementations then try to compensate independently:

1. `annotation.apply.annotate_obs()` joins `obs` to `raw_file` exactly.
2. `proteobench.intermediate.align_runs()` tries cleaned `raw_file`,
   `sample_name`, and a limited wide-rule heuristic.

The ProteoBench heuristic happens to make WOMBAT peptidoform scoring work,
because `abundance_A_1` matches the x-layer regex and yields `A_1`. It cannot
resolve WOMBAT ion or PEAKS diaPASEF, and annotation does not use the heuristic
at all.

The root defect is therefore not the WOMBAT or PEAKS parsing rule. It is the
missing declarative run crosswalk plus duplicated sample-resolution logic.

## Intended design

Keep vendor-derived observation identifiers unchanged. Add one shared,
strict sample resolver used by both annotation and ProteoBench scoring.

Inputs:

- converted observations and the stored effective `ParseRule`;
- the module `[[samples]]` design;
- an optional upstream `[run_mapper]`, whose values name
  `[[samples]].sample_name`.

Resolution order should mirror the pinned ProteoBench contract:

1. explicit per-tool `run_mapper`;
2. cleaned `[[samples]].raw_file`;
3. exact `[[samples]].sample_name`.

For a wide APB rule, normalize each mapper key into the observation namespace:

1. match it against the x-layer source regex and use the `sample` group when
   present;
2. otherwise apply the standard run-name cleanup.

This converts `abundance_A_1` to `A_1`, while leaving the PEAKS `LFQ_…` mapper
keys intact. Resolve the mapper value through the unique module
`sample_name`.

The resolver must reject:

- aliases that normalize to the same observation but point at different
  samples;
- multiple observations resolving to one module sample;
- partial module coverage;
- unknown mapper values;
- ambiguous cleaned raw-file or sample names.

No fuzzy or substring matching.

## Implementation tasks

### 1. Acquire the authoritative mapping

- Extend `apb-testdata annotations` to fetch the matching per-tool parse
  settings from the same pinned ProteoBench revision as the module settings.
- Store them by `(module, software)` rather than as a module-only resource.
- Validate `[run_mapper]` as `dict[str, str]`.
- Do not copy mapper entries into APB's global vendor parsing rules: these
  mappings are module-specific.

### 2. Add one alignment API

- Introduce a typed run-mapping model and a shared resolver in APB.
- Move common cleanup, one-to-one validation, and completeness checks into
  that resolver.
- Replace `proteobench.intermediate._wide_run_mapping()` with the shared
  implementation.
- Make annotation consume the same resolved sample records instead of joining
  directly on `raw_file`.

### 3. Thread the resource explicitly

Add an optional explicit mapping argument to both commands:

```text
apb annotate INPUT MODULE_SETTINGS --run-mapping PARSE_SETTINGS --output OUTPUT
apb proteobench INPUT MODULE_SETTINGS --run-mapping PARSE_SETTINGS --output OUTPUT
```

Thread the path through:

- APB CLI and Python pipeline APIs;
- APB Studio fixture/resource resolution;
- frozen `RunSnapshot`;
- registry commands and Snakemake inputs;
- artifact command/provenance display.

Identity-matching tools must continue to work without `--run-mapping`.

### 4. Preserve provenance

Record:

- mapping source path;
- selected resolution strategy;
- observed identifier → module `sample_name` mapping;
- matched/expected run counts.

Store this once in APB metadata and expose it in `apb summary`.

### 5. Tests

Unit-test the three confirmed crosswalks:

- WOMBAT ion: `A_1` → `Condition_A_Sample_Alpha_01`;
- WOMBAT peptidoform: `A_1` → `Condition_A_Sample_Alpha_01`;
- PEAKS diaPASEF:
  `LFQ_ttSCP_diaPASEF_Condition_A_Sample_Alpha_01`
  → `Condition_A_Sample_Alpha_01`.

Also test:

- identity matching without a mapper;
- AnnData and MuData annotation;
- ProteoBench scoring with the same resolver;
- collision, ambiguity, unknown target, and incomplete-coverage errors;
- CLI command rendering and frozen Studio snapshots.

## Acceptance criteria

- The three mappings above resolve six of six samples, one-to-one.
- The affected 8 annotation and 6 ProteoBench targets succeed.
- The corpus reaches 541 of 541 targets with no sample-name failures.
- Converted `obs_names` remain vendor-derived and unchanged.
- Annotation and scoring report the same canonical `sample_name` and
  `condition`.
- The exact APB command and run-mapping source are visible in artifact
  provenance.

## Non-solutions

- Do not hard-code module names in WOMBAT or PEAKS parsing rules.
- Do not rewrite `obs_names` with ad hoc prefix/suffix stripping.
- Do not edit downloaded `module_settings.toml` files to impersonate vendor
  output.
- Do not add separate annotation-only and ProteoBench-only matching fixes.
- Do not silently accept partial or ambiguous mappings.

## Separate PEAKS data defect

Fixture `5691d485` contains malformed auxiliary headers such as:

```text
LFQ_ttSCP_diaPASEF_Condition_A_Sample_LFQ_ttSCP_diaPASEF_Condition_B_Sample_Alpha_02 m/z
```

The x-layer `Normalized Area` headers are correct, so this does not cause the
14 sample-alignment failures. APB correctly excludes the malformed auxiliary
sample token from the observation axis and logs a warning. Track correction
of that source export separately if complete per-run m/z/RT layers are
required.
