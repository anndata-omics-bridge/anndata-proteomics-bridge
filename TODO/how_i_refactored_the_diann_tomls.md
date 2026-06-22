# How I Refactored the DIA-NN TOMLs

Date: 2026-06-22

## Problem

The DIA-NN peptide and peptidoform TOMLs were creating standalone AnnData modalities by
aggregating precursor-level layers:

- `Precursor.Quantity`
- `Precursor.Normalised`

That was wrong for the parser contract. Those columns are ion/precursor-level measurements in the
DIA-NN report. Reusing them under peptide or peptidoform modalities created new quantitative data
during conversion instead of only representing what the vendor output actually contains.

The corrected principle is:

> Capture all data present in the vendor output, but do not create new quantitative data in the
> parsing conversion.

## What Changed

I removed the active DIA-NN peptide and peptidoform parsing rules:

- `src/anndata_proteomics/parsing_rules/diann/parse_diann_peptide.toml`
- `src/anndata_proteomics/parsing_rules/diann/parse_diann_peptidoform.toml`

This de-registers those levels at the rule-discovery layer. Nothing special is hidden in the GUI.
The GUI and MuData conversion now see only levels with valid report-backed TOMLs.

Current DIA-NN level resolution:

- DIA-NN 1.x: `ion`, `protein`, `fragment`
- DIA-NN 2.x: `ion`, `protein`

`mudata` remains available when at least two report-backed levels resolve.

## What Stayed

Peptide and peptidoform identifiers are still mandatory metadata. They are not optional extras.
They remain in `.var` when the valid rule selects or computes them from vendor output columns.

For DIA-NN ion-level output, the rule still preserves:

- `ProForma_peptidoform`
- `ProForma_peptide`
- `ProForma_ion`
- `Protein_Group`
- `Protein_Ids`
- `Protein_Names`
- `Genes`

The important distinction is that these are metadata/link columns on a valid measurement level.
Their existence does not imply that APB should create a standalone peptide or peptidoform
abundance matrix.

## Why This Is the Root Fix

The previous TOMLs explicitly said the peptide and peptidoform quantities were APB-derived
rollups. That made the parsing rules do two jobs:

1. parse the source format
2. derive new quantitative levels

Those must stay separate. A parsing TOML should describe the vendor file. If APB later needs
derived peptide or peptidoform quantification, that should be a named derivation pipeline with an
explicit algorithm, not a parsing-rule side effect.

## Tests Updated

The tests now encode the corrected behavior:

- `tests/test_rule_resolve.py` expects DIA-NN v1 to resolve `ion`, `protein`, `fragment` and DIA-NN
  v2 to resolve `ion`, `protein`.
- `tests/test_diann_levels.py` converts only report-backed DIA-NN levels and asserts the ion
  `.var` still contains peptide, peptidoform, and protein metadata.
- `tests/test_mudata_levels.py` builds MuData from report-backed levels only and checks
  fragment-to-ion/protein metadata links without expecting peptide or peptidoform modalities.
- `tests/conftest.py` now finds cached DIA-NN fixtures by the corrected report-backed level set.

## Docs Updated

I updated:

- `TODO/HOWTO/test_gui.md`
- `TODO/TODO_ui_test_tool.md`
- `TODO/TODO_to_mu_data.md`

The docs now say that standalone levels are driven by real vendor output layers, while identifiers
and biological links are preserved as `.var` metadata.

## Skill Added

I created a reusable skill documenting this rule. The single source of truth is the
`claude-kaiser-skills` repo; both Claude Code and Codex consume it via symlinks managed by the
skill coordinator — there is no hand-placed copy.

- source of truth:
  `/Users/wolski/projects/wews_skill_coordinator/repos/claude-kaiser-skills/apb-toml-level-design/`
  (`SKILL.md` + `agents/openai.yaml`)
- registered in `/Users/wolski/projects/wews_skill_coordinator/skills.toml` under
  `[skills.claude-kaiser-skills].paths`
- wired into both harnesses by `make install` (`skill_coordinator.py` mirrors every skill to
  `~/.claude/skills/` and `~/.codex/skills/`):
  - `~/.claude/skills/apb-toml-level-design` → repo source (symlink)
  - `~/.codex/skills/apb-toml-level-design` → repo source (symlink)

To pick up edits or re-create the symlinks, run `make install` in the coordinator repo. Do not
edit or copy files under `~/.codex/skills/` or `~/.claude/skills/` directly — they are symlinks.

Core rule from the skill:

> A standalone AnnData level is valid only when the source output has one or more real quantitative
> layer columns for that level. Metadata/link columns such as `ProForma_peptide` are mandatory, but
> they do not prove that a standalone quantitative level exists.

## Verification

Commands run after the refactor:

```bash
uv run pytest tests/test_rule_resolve.py tests/test_diann_levels.py tests/test_mudata_levels.py tests/test_ui_support_converted_runs.py
uv run pytest tests/test_packaged_rules.py tests/test_json_schema_validation.py
uv run ruff check src/anndata_proteomics/scripts/_ui_support.py tests/test_rule_resolve.py tests/test_diann_levels.py tests/test_mudata_levels.py tests/conftest.py
uv run --extra gui python -c "import anndata_proteomics.scripts.ui_test_tool"
uv run --with pyyaml python /Users/wolski/.codex/skills/.system/skill-creator/scripts/quick_validate.py /Users/wolski/projects/wews_skill_coordinator/repos/claude-kaiser-skills/apb-toml-level-design
# then wire it into both harnesses:
make -C /Users/wolski/projects/wews_skill_coordinator install
```
