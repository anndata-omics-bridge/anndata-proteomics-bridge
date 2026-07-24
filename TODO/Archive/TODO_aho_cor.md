# APB FASTA validation with Aho--Corasick

> Implemented and archived 2026-07-24. Optional protein inference remains a
> separate milestone in `../TODO_protein_inference.md`.

## Goal

`apb fasta` annotates proteins and automatically validates every peptide-derived
identification against the exact supplied FASTA database. Validation is
annotation only: unmatched peptides, decoys, contaminants, and their
quantifications remain in the object.

The backend-neutral Aho--Corasick and protein-inference algorithms live in the
standalone Python `prozor` package. APB depends on Prozor, while APB's streaming
FASTA parser and AnnData/MuData storage remain consumer-specific concerns.

## Accepted behavior

- `apb fasta DATA FASTA...` performs peptide validation by default.
- `--no-validate` disables peptide validation while retaining protein annotation.
- The separate `apb validate-fasta` command is removed.
- Protein-only AnnData is annotated; peptide-only AnnData is validated; MuData
  receives both operations wherever the corresponding modalities exist.
- All peptide-derived MuData modalities (`ion`, `fragment`, `peptidoform`, and
  `peptide`) are validated with one union of unique `ProForma_peptide` patterns
  and one FASTA scan.
- No unmatched-fraction gate filters or rejects ordinary output. Validation
  results are evidence, not a quantification filter.

## Per-feature validation

Each peptide-derived modality stores a var-aligned DataFrame in
`varm["fasta_validation"]`:

| Column | Meaning |
| --- | --- |
| `peptide_in_fasta` | The unmodified peptide occurs at least once in the supplied FASTA. |
| `fasta_match_site_count` | Total occurrence sites, including repeated sites in one protein. |
| `fasta_matching_protein_count` | Number of distinct matching FASTA accessions. |
| `fasta_matching_protein_ids` | Sorted, semicolon-separated distinct matching accessions. |
| `leading_protein_in_fasta` | The reported leading protein exists in the FASTA; nullable when no leading-protein field exists. |
| `peptide_in_leading_protein` | The peptide occurs specifically in the reported leading protein; nullable when it cannot be evaluated. |

Site count and distinct-protein count are deliberately separate. Multiple sites
within one protein are useful evidence and are not collapsed.

The in-memory result retains the complete site table (`sequence`, raw FASTA ID,
clean accession, half-open coordinates, length, decoy flag, contaminant flag).
It is not serialized as a nested AnnData or as a ragged table in `uns`.

## Decoy and contaminant configuration

FASTA identifier handling is a typed configuration stored as JSON under
`uns["anndata_proteomics"]["fasta_config"]` on the root AnnData/MuData.

- Omitted patterns are inferred by scanning raw FASTA identifier tokens against
  conservative configured candidates (for example `^REV_`, `^rev_`, and
  `^zz(?:\\||_)`).
- Explicit patterns override inference; an explicit empty pattern disables that
  classification.
- Resolved patterns, whether they were explicit/inferred, and per-candidate hit
  counts are persisted and copied into operation provenance.
- Classification never removes FASTA records or quantified rows.
- APB never attempts to regenerate decoy sequences; validation uses exactly the
  sequences present in the supplied FASTA.
- Repeated raw FASTA identifiers do not suppress later records; every supplied
  sequence is scanned.
- Decoy prefixes preceding `sp|`/`tr|` are preserved when deriving clean
  accessions (`REV_sp|P12345|...` becomes `REV_P12345`, not `P12345`).

Protein annotation keeps targets, decoys, and contaminants and adds
`is_decoy`/`is_contaminant`. Duplicate clean accessions use deterministic
target-before-decoy, `sp|`-before-`tr|`, then input-order precedence.

## MuLink representation

For MuData with an existing protein modality, representable relationships are
stored as a signed-integer CSR adjacency matrix at
`mdata.varp["feature_mapping"]`, following MuLink's `u -> v` convention:

`peptide-derived feature -> existing protein feature`

- Unique-sequence matches are expanded back to every feature row.
- Edges from all peptide-derived modalities are combined.
- Existing `feature_mapping` edges and weights are preserved. APB tracks only
  its additive FASTA contribution in
  `varp["_apb_fasta_feature_mapping_contribution"]`, so revalidation replaces
  APB-owned edges for the selected modalities without leaving stale links or
  deleting relationships owned by another producer.
- Only protein nodes already present on the MuData global var axis can receive
  edges. Matches to other FASTA proteins remain visible in
  `fasta_matching_protein_ids` and are counted in provenance.
- Standalone AnnData and MuData without a protein modality still receive complete
  per-feature validation but no invented protein modality.

APB writes the MuLink-compatible standard representation directly. Consumers
may import `mulink` for its query/plot accessors; APB does not need its heavier
runtime dependency merely to populate `varp`.

## Sequence rules

- Match the unmodified `ProForma_peptide` field; do not strip modification text
  ad hoc from another field.
- Upper-case both sides.
- I/L equivalence is opt-in and recorded in provenance.
- Ambiguous residues remain literal.
- This is sequence-presence validation, not digestion validation; enzyme
  specificity and missed cleavages do not gate a match.

## Future TODO: optional protein inference in APB

Protein inference remains a separate opt-in milestone. It may consume the
peptide--protein evidence and Prozor's greedy-parsimony implementation, but
must not overwrite vendor protein groups or run automatically during FASTA
validation. Before exposing it:

- define how unrepresented FASTA proteins become queryable without inventing a
  misleading quantified modality;
- define unique/shared/razor, decoy, contaminant, subset, and indistinguishable
  protein behavior;
- make tie-breaking deterministic and record algorithm/version provenance;
- compare against a fixed Prozor fixture.

## Acceptance checks

- Automatic/default and `--no-validate` CLI behavior.
- Protein-only, peptide-only, and multi-modality MuData inputs.
- Shared peptides and repeated sites.
- Leading-protein existence and sequence consistency.
- Decoy-only and contaminant-only identifications retained and valid when their
  sequences occur in the supplied FASTA.
- Pattern inference, explicit override, empty override, and JSON round-trip.
- MuLink orientation, multi-modality union, scoped revalidation, preservation of
  existing links and exact weights, large sparse axes, and `.h5mu` round-trip.
- Pure-Python/Rust backend equivalence and invalid backend rejection.
- Generator/path FASTA provenance remains intact.
