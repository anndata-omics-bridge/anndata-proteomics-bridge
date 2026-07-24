# APB high-level roadmap

> Close the remaining APB TODOs while keeping APB a lossless, JSON-configured bridge into the scverse ecosystem.

Date: 2026-07-20  
Status: active roadmap

## Strategic direction

APB converts proteomics vendor outputs into level-specific AnnData and MuData objects. The next
stage is to make the existing conversion reliable across the corpus, formalize its input contracts,
and represent relationships between proteomics levels without inventing an APB-specific ecosystem.

QPX and alphabase remain useful prior art for schemas, vendor mappings, and modification handling.
For cross-level feature linking, APB will align with **mulink**, not QPX.

## Roadmap

### 1. Complete conversion coverage

Resolve the known FragPipe and PEAKS corpus failures and verify the results at artifact level. A
successful conversion must have correct samples, feature identities, quantitative layers, and
provenance—not merely avoid an exception.

Source: [TODO_failures_260702.md](TODO_failures_260702.md)

### 2. Formalize the vendor-column contract

Add a clear contract for required, optional, and aliased vendor columns. Preserve existing rules,
produce stable missing-column errors, and keep the selected vendor source visible in provenance.

QPX is prior art for explicit required fields; alphabase is prior art for column aliases. Neither
project replaces APB's JSON rule model.

Source:
[TODO_column_extraction_required_columns_alphapep_qpx_2026-07-07.md](TODO_column_extraction_required_columns_alphapep_qpx_2026-07-07.md)

### 3. Align proteomics hierarchies with mulink

Use mulink as the intended scverse representation for relationships among fragment, ion/precursor,
peptidoform, peptide, protein, and gene modalities.

Before committing APB to a serialized linking convention:

- contact mulink maintainer Lucas Diedrich;
- align on edge direction, direct versus transitive links, relation provenance, validation, and
  serialization conventions;
- establish a small shared proteomics example;
- avoid creating a competing APB graph API.

The local mulink checkout is at `related_work/mulink`.

Outreach draft: [DRAFT_email_mulink_developer.md](DRAFT_email_mulink_developer.md)

### 4. Consolidate modification handling

Improve result-table modification correctness using alphabase's mature mapping knowledge while
retaining APB's accession-backed ProForma output and lossless behavior. The strategic concerns are
vendor token coverage, mass-based matching, fixed modifications, and one shared registry-backed
model rather than overlapping mechanisms.

Parameter-side registry enrichment remains deferred until SDRF export or another concrete consumer
needs the richer typed metadata.

Sources:

- [TODO_modification_alphabase_2026-07-07.md](TODO_modification_alphabase_2026-07-07.md)
- [TODO_modification_homogenization_design.md](TODO_modification_homogenization_design.md)

### 5. Clarify the search-parameter contract

Decide which parameter fields APB can require across vendors and which are legitimately absent.
Keep vendor-specific syntax at the parser boundary and canonical values in the shared model. The
contract must preserve APB's deliberate ability to retain quantification data when an attached
parameter file is missing, incomplete, or belongs to the wrong tool.

Source: [TODO_params_model_review.md](TODO_params_model_review.md)

### 6. Validate identifications against FASTA — delivered

APB now performs scalable peptide-presence validation against the searched FASTA database and
preserves complete peptide-protein match evidence for mulink-compatible hierarchy edges.

Optional protein inference remains a separate active milestone.

Sources:

- [Archived FASTA validation specification](Archive/TODO_aho_cor.md)
- [Optional protein inference follow-up](TODO_protein_inference.md)

## Deferred work

- Registry-backed enrichment of searched-parameter modifications remains deferred until consumed.
- APB-derived protein inference remains separate from FASTA validation and must not overwrite
  vendor protein assignments.
- Additional APB Studio workflow features remain demand-driven rather than roadmap commitments.

## Roadmap completion

This roadmap is complete when the linked active TODOs have either been delivered or deliberately
retired, mulink alignment is documented, and the active TODO directory contains no addressed work.
