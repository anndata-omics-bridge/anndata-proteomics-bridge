# Completed: align APB and APB Studio with the FGCZ Python quality reference

> Completed and archived 2026-07-24. APB and APB Studio enforce the aligned
> local and CI quality contract; their full commit, push, security, package, and
> documentation gates pass.

## Goal

Adopt the complete local and CI quality system from
`/Users/wolski/projects/fgcz_python_project_reference/` in both `apb/` and
`apb_studio/`, adapted only for their package names, Python 3.13 environments,
console entry points, documentation, and wheel contents.

The result must be operational, not a set of copied configurations that are
known to fail.

## Reference contract to reproduce

### Pre-commit stage

- Ruff lint with the reference rule set and complexity limits.
- Ruff formatting.
- Pyright strict.
- Deptry dependency-declaration validation.

### Pre-push stage

- Pytest with line and branch coverage.
- 100% changed-line coverage through `diff-cover`.
- Source-distribution/wheel build and package-contract inspection.
- Strict MkDocs build.

### Manual and scheduled stage

- `pip-audit` vulnerability scan over the locked local environment.

### GitHub Actions

- `ci.yml` runs the same pre-commit and pre-push hooks from a frozen uv
  environment on pushes to `main`, pull requests, and manual dispatch.
- `pages.yml` builds documentation strictly and deploys GitHub Pages on
  `main`.
- `security.yml` runs the same manual dependency-audit hook weekly and on
  demand.
- Workflow permissions remain minimal and action versions follow the reference.

## Current baseline

### APB

- Existing tests: 487 passed, 4 skipped.
- Existing strict Pyright and Deptry gates pass.
- Reference Ruff policy currently reports 112 findings.
- Line/branch coverage is 85%; the reference gate requires 100%.
- MkDocs exists, but CI currently contains only a separate `docs.yml`.
- Missing changed-line coverage, package smoke, vulnerability audit, and the
  mirrored CI workflow.

### APB Studio

- Existing tests: 148 passed.
- Reference Ruff policy currently reports 144 findings.
- Strict Pyright, when run with the Studio Python 3.13 environment, reports 59
  errors.
- Deptry is not configured; its raw baseline reports 34 findings, most of which
  are missing first-party/dev-group configuration plus several real dependency
  declarations.
- Coverage, Pyright, Deptry, pre-commit, audit, diff-cover, and MkDocs tooling
  are absent from the declared development environment.
- No documentation site or GitHub Actions workflows exist.

## Implementation plan

### 1. Establish the same dependency and configuration surface

- Add the reference quality tools to each repository's locked development
  environment: `pre-commit`, `pyright`, `deptry`, `pytest-cov`, `diff-cover`,
  and `pip-audit`.
- Add/retain the MkDocs Material documentation group.
- Keep each repository's existing Python 3.13 requirement and uv source
  declarations.
- Configure strict pytest markers/options, Ruff lint/complexity/format rules,
  strict Pyright, line/branch coverage, and Deptry first-party/package maps.
- Refresh and commit each `uv.lock`.
- Add `py.typed` to each shipped package and verify it is included in wheels.

### 2. Make the pre-commit gates pass

- Replace the current hook layouts with the reference staged layout.
- Fix Ruff findings mechanically where safe, then refactor root causes for
  complexity and argument-count findings without changing public APIs.
- Fix Studio's strict typing errors at their source. Use typed protocols or
  narrowly typed adapters for injected subprocess/Dash test objects; do not add
  ignore comments or weaken strict mode.
- Resolve Deptry findings through correct runtime/dev dependency ownership and
  documented narrow ignores only for genuinely dynamic CLI/plugin use.
- Preserve APB's existing prohibition on type-checker ignore comments.

### 3. Reach the coverage contract

- Add focused tests for uncovered behavior and error paths in both packages.
- Set line and branch coverage to the reference's 100% blocking threshold.
- Add `scripts/diff_coverage.py` in both repositories with
  `DIFF_COVER_BASE=main` locally and the pull-request base branch in CI.
- Do not use coverage exclusions, `pragma: no cover`, or broad omit patterns to
  manufacture the result; only structurally unexecutable module-entry guards
  may be excluded if the reference policy requires it consistently.

### 4. Add package-smoke contracts

- Adapt `scripts/package_smoke.py` for each distribution.
- APB smoke checks:
  - the `anndata_proteomics` package and `py.typed`;
  - packaged parsing-rule/schema/license assets;
  - `apb` and `apb-testdata` console entry points.
- APB Studio smoke checks:
  - the `apb_studio` package, registry/workflow assets, and `py.typed`;
  - the two preferred application entry points and compatibility aliases.
- Build into a temporary directory with `uv build`; never leave `dist/`
  artifacts in the repository.

### 5. Make strict documentation builds available

- Retain APB's existing documentation and make its strict build the pre-push
  and Pages source of truth.
- Replace APB's standalone `docs.yml` with the reference-aligned `pages.yml` to
  avoid duplicate Pages deployments.
- Add a minimal maintained APB Studio MkDocs site based on its existing README
  and architecture/status contracts, with development/check instructions.
- Add scoped documentation/workflow agent guidance where the reference relies
  on it.

### 6. Add mirrored GitHub Actions

- Add reference-aligned `ci.yml`, `pages.yml`, and `security.yml` to both
  repositories.
- Use frozen lockfiles and install the exact local hook environments required
  by each project.
- Keep local hooks as the command source of truth; workflows invoke hooks
  rather than duplicate check commands.

### 7. Update developer interfaces and documentation

- Update `AGENTS.md`, `CLAUDE.md`-linked development instructions, README,
  Makefile targets, and `CHANGES.md` so the canonical commands are:
  - fast gate:
    `uv run pre-commit run --hook-stage pre-commit --all-files`;
  - full gate:
    `uv run pre-commit run --hook-stage pre-push --all-files`;
  - manual audit:
    `uv run pre-commit run dependency-audit --hook-stage manual --all-files`.
- Keep APB and APB Studio naming and dependency direction unchanged.

## Validation

For each repository:

1. `uv sync --frozen` with all quality/docs groups required by hooks.
2. `uv run pre-commit run --hook-stage pre-commit --all-files`.
3. `uv run pre-commit run --hook-stage pre-push --all-files`.
4. `uv run pre-commit run dependency-audit --hook-stage manual --all-files`.
5. Parse all workflow YAML and inspect the built wheel contract.
6. Confirm a clean tracked worktree, excluding the existing generated
   `test_data_download/` and `apb_outputs/` caches.

## Done when

- Both packages expose every reference check and all checks pass locally.
- Local pre-commit stages and GitHub Actions invoke the same commands.
- CI, Pages, and scheduled security workflows exist in both repositories.
- Strict typing, dependency validation, full line/branch coverage,
  changed-line coverage, package smoke, strict docs, and vulnerability audit
  are blocking at the same stages as the reference.
- No public APB/APB Studio API or dependency direction is changed merely to
  satisfy tooling.
