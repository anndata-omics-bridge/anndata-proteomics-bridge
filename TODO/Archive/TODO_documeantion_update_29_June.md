# ARCHIVED: documentation update and reconciliation, 29 June

Status: completed by commit `2636a8f`; archived 2026-07-20.

## Current cleanup decision

Moved addressed or historical APB notes out of the active TODO surface:

| Moved to `TODO/Archive/` | Why it is no longer active |
|---|---|
| `TODO_fasta.md` | Implemented in code, CLI, tests, README, and architecture docs. |
| `TODO_ui_test_tool.md` | APB no longer carries this GUI/runtime surface; keep as historical context only. |
| `REVIEW_toml_modifications.md` | Its findings are marked implemented and the behavior is now covered by tests/docs. |
| `REVIEW_qpx.md` | Prior-art review artifact, not an active implementation TODO. |
| `REVIEW_qpx_vs_apb.md` | Prior-art review artifact, not an active implementation TODO. |

`TODO/REVIEW_qpx.html` was ignored generated output, not a source artifact, and was removed from
the working tree.

The active APB TODO surface should stay small. After this cleanup, the only pre-existing active
note is `TODO_modification_homogenization_design.md`, which is intentionally deferred future work.

## Documentation sources to reconcile

- `README.md` -- public quick-start and short user-facing feature summary.
- `docs/index.md` -- MkDocs landing page and routing surface.
- `docs/ARCHITECTURE.md` -- current module map and data-flow diagrams.
- `docs/toml_schema.md` -- parsing-rule authoring contract.
- `docs/parameter_parsers.md` -- vendor parameter parser contract and caveats.
- `docs/parsing_architecture.md` -- lower-level UML/flow diagrams.
- `mkdocs.yml` and the root `make docs` target -- documentation build contract.
- `AGENTS.md` -- maintainer rules; should agree with the public docs but not duplicate them.

## Known reconciliation problems

1. **MkDocs warnings for source/test links.** The build passes, but `ARCHITECTURE.md` and
   `parameter_parsers.md` still contain relative links to files outside `docs/`, which MkDocs warns
   about because they are not documentation pages. Decide one convention:
   - use GitHub source links for source/tests, or
   - keep plain code paths without markdown links.

2. **R-side report package section is stale.** `docs/ARCHITECTURE.md` names
   `~/projects/anndata_bridge/annProtSum/`, but that path is absent. The workspace currently has
   `annProtSum_to_be_dropped/` and `annProtSum_0.0.0.9000.tar.gz`. Treat the R report package as
   non-canonical until a decision is made:
   - drop the section from APB docs, or
   - rewrite it as historical/experimental, or
   - revive it as a maintained package with current APB paths and data model.

3. **R report assumptions are behind current APB.** The stale package/docs still refer to old
   paths such as `anndata_proteomics_bridge`, `tools/generate_h5ads.py`, and `examples/results/`.
   They are also centered on `.h5ad` reports, while APB now commonly writes multi-level `.h5mu`
   and stores FASTA annotation in the protein modality. Do not use that package as a source of
   truth for APB docs until it is either removed or updated.

4. **README and docs overlap.** The README now contains substantial TOML, CLI, modifications, and
   limitations content. Decide what stays in README as a concise entry point and what moves to the
   MkDocs pages.

5. **Architecture docs mix current state and future notes.** `ARCHITECTURE.md` has a "Not yet
   implemented" section plus implementation details. Confirm each item against current code and
   move stale/future material into TODO files or remove it.

6. **Parameter/modification docs need one canonical story.** README, `docs/toml_schema.md`,
   `docs/parameter_parsers.md`, and `docs/parsing_architecture.md` all discuss modifications.
   Reconcile them so:
   - parsing-rule modifications are documented in `toml_schema.md`;
   - parameter-file modification normalization is documented in `parameter_parsers.md`;
   - architecture pages explain only the module flow and data model.

7. **Generated outputs must stay out of git.** Keep `public/`, `docs/*.html`, and `TODO/*.html`
   ignored. Do not restore Pandoc-era generated HTML.

## Update plan

1. **Run a docs inventory.**
   - Build with `uv run --frozen --group docs mkdocs build`.
   - Capture all warnings.
   - Search for stale paths and names:
     `anndata_proteomics_bridge`, `generate_h5ads`, `generate_report.py`, `examples/results`,
     `annProtSum`, old HTML outputs, and obsolete CLI names.

2. **Fix broken or noisy links.**
   - Convert source/test links that MkDocs cannot resolve into GitHub links or plain paths.
   - Keep documentation-page links relative inside `docs/`.
   - Rebuild and require no unexpected MkDocs warnings.

3. **Reframe README as the entry point.**
   - Keep install, common CLI examples, the published docs link, and a compact feature list.
   - Move detailed TOML and modification exposition into MkDocs pages where appropriate.
   - Keep README examples synchronized with the actual `apb` CLI.

4. **Make `docs/index.md` the navigation contract.**
   - Add one-line descriptions for each page.
   - Ensure it points users to README for quick start and to specific docs pages for authoring
     rules, parser behavior, and architecture.

5. **Rewrite the R-side report section.**
   - First decide whether the R package is dead, experimental, or maintained.
   - If dead/experimental, remove it from the current architecture path and mention only that
     report generation is not part of the current APB public contract.
   - If maintained, update it for the current package name, workspace paths, `.h5mu` handling,
     protein FASTA `varm['fasta']`, and current search-parameter storage.

6. **Separate current architecture from future work.**
   - Keep `ARCHITECTURE.md` factual and current.
   - Move unresolved future items into active TODO files, or archive them if already done.
   - Preserve Mermaid diagrams only where they clarify current flow.

7. **Reconcile modification documentation.**
   - Ensure `toml_schema.md` documents only supported parsing-rule modification modes.
   - Ensure `parameter_parsers.md` documents parameter-file modification normalization and the
     deferred registry-backed design as future work only when necessary.
   - Ensure `parsing_architecture.md` diagrams match the actual modules.

8. **Validate and publish.**
   - Run `uv run --frozen --group docs mkdocs build`.
   - Run `git diff --check`.
   - Confirm generated HTML remains ignored.
   - Commit and push the APB documentation-only update.
   - Check the GitHub Pages workflow and published site.

## Acceptance criteria

- Active `TODO/` contains only live work.
- Archived notes remain available under `TODO/Archive/`.
- MkDocs renders Mermaid diagrams and builds without new warnings.
- README links clearly to the published GitHub Pages site.
- Public docs do not present stale R-report-package paths as current architecture.
- No generated HTML or `public/` output is tracked.
