# Optional protein inference from FASTA evidence

**Status:** open follow-up
**Extracted:** 2026-07-24 from the completed FASTA-validation specification

## Goal

Expose Prozor's deterministic greedy-parsimony inference as an explicit,
opt-in APB operation built on stored peptide-to-protein evidence.

## Decisions required

- Define how unrepresented FASTA proteins remain queryable without inventing a
  misleading quantified modality.
- Define unique/shared/razor, decoy, contaminant, subset, and
  indistinguishable-protein behavior.
- Specify deterministic tie-breaking and record algorithm/version provenance.
- Decide the storage contract without overwriting vendor protein groups.

## Acceptance

- Compare inference against a fixed Prozor reference fixture.
- Preserve vendor protein assignments and all quantified rows.
- Keep inference separate from default FASTA validation.
- Pass APB's full quality gate.
