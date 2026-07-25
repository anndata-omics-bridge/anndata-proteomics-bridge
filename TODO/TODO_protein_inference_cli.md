# Protein inference — a new, opt-in APB CLI command

**Status:** open follow-up (new tool, not yet started)
**Scope clarified:** 2026-07-24

## What this is (and is not)

This is about adding a **new, standalone APB CLI command** (working name
`apb infer-proteins`, this file: `protein_inference_cli`) that exposes Prozor's
deterministic greedy-parsimony protein inference as an explicit, opt-in operation.

It is **not** about the FASTA/mapping work, which is already done:

- Prozor already implements the algorithm (`prozor.greedy.greedy_parsimony`,
  `GreedyResult`) — per AGENTS.md, peptide-to-protein algorithms live in Prozor.
- APB already validates peptides against the FASTA and **already stores the
  peptide-to-protein match matrix** — `varp['feature_mapping']` (MuData
  adjacency) and `varm['fasta_validation']` (per-feature evidence).

So the algorithm and the evidence both exist. What is missing is the APB-side
**command** that runs inference and writes its result.

## Design still open

- **Does the command consume the stored mapping, or run standalone?** It *might*
  build its `PeptideProteinMatrix` from the already-stored `varp['feature_mapping']`,
  or it *might* infer independently from peptides + FASTA. This is deliberately
  left open — decide when the command is designed; do not assume validation ran first.
- Storage contract: where inferred groups land, **without overwriting vendor
  protein assignments** or the stored mapping.
- Which semantics APB surfaces vs. what Prozor already handles internally
  (Prozor's greedy already does indistinguishable grouping + subset `subsume`;
  decoy/contaminant filtering is the caller's job — decide whether the command
  filters or defers).
- Deterministic tie-breaking is already Prozor's contract — just record the
  algorithm/version in provenance.

## Acceptance

- New CLI command with its own tests; compare output against a fixed Prozor
  reference fixture.
- Preserve vendor protein assignments, the stored mapping, and all quantified rows.
- Inference stays separate from (and optional relative to) default FASTA validation.
- Pass APB's full quality gate.
