# PLAN 2026-05-13 — FASTA protein annotation module

## Goal

Add a pure-Python module that parses a UniProt-style FASTA file and
returns a `pandas.DataFrame` of protein annotations, replicating
prolfquapp's `get_annot_from_fasta()` so APB can offer the same
data to downstream R consumers (prolfquapp) without round-tripping
through R.

Reference: [prolfqua_fml/prolfquapp/R/get_annot_from_FASTA.R](file:///Users/wolski/projects/prolfqua_fml/prolfquapp/R/get_annot_from_FASTA.R)

## Output contract

`fasta_to_dataframe(paths, *, ...)` returns a `pd.DataFrame` with
one row per non-decoy protein record and the columns below.

| Column | Source | Notes |
|---|---|---|
| `fasta.id` | First whitespace-delimited token of the header, with `>` and trailing `;` stripped | Always |
| `fasta.header` | Remainder of the header after the first whitespace | Always |
| `proteinname` | Middle pipe field of `fasta.id` when `is_uniprot=True` (`sp|ACC|NAME` → `ACC`); else equals `fasta.id` | Always |
| `gene_name` | Regex `GN=(.+?) PE=` extracted from `fasta.header` | Added only if >1 records match — matches prolfquapp's gating |
| `protein_length` | `len(sequence)` | Always |
| `nr_tryptic_peptides` | Number of fully-tryptic peptides between K/R cleavage sites (K\|R not followed by P), counted with `min_length ≤ L < max_length` (defaults 7 ≤ L < 30) | Always |
| `sequence` | Raw amino-acid string | Only when `include_sequence=True` |

Decoy records with `fasta.id` matching `decoy_pattern` (default
`r"^REV_\|^rev_"`) are filtered out before any downstream column is
computed. Duplicate `fasta.id` values are deduplicated (first wins).
If `<10%` of records match the decoy pattern, log a warning (matches
prolfquapp's heuristic).

## Module layout

```
src/anndata_proteomics/fasta/
  __init__.py        # empty (per CLAUDE.md rule)
  parser.py          # FastaRecord, iter_fasta(path) → Iterator[FastaRecord]
  annotation.py      # fasta_to_dataframe(paths, **opts) → pd.DataFrame
                     # plus _find_cleavage_sites, _count_tryptic_peptides,
                     # extract_gene_name
```

Public surface (importable from `anndata_proteomics.fasta`):
- `fasta_to_dataframe(paths, *, decoy_pattern, is_uniprot,
  min_length, max_length, include_sequence) -> pd.DataFrame`
- `extract_gene_name(header) -> str` (per-row helper)
- `count_tryptic_peptides(sequence, *, min_length, max_length) -> int`

Keep it minimal — no `ProteinAnnotation` runtime class, no R6
emulation. That class is a wrapper around an LFQData object; the
data-extraction half is the only thing APB needs to own.

## Implementation choices

- **No biopython dependency.** A plain text reader is enough; we
  don't need SeqRecord, alphabets, etc. Adding biopython for one
  parser is overkill.
- **Tryptic peptide counting** uses prolfquapp's exact rule: cleave
  after K or R not followed by P; count peptides with length in
  `[min_length, max_length)` (note: prolfquapp uses `< max_length`,
  not `<=`, despite the docstring saying "maximum length 30" — the
  code is authoritative).
- **Gene-name gating** matches prolfquapp: the `gene_name` column
  is added only when more than one record produces a match. This
  avoids polluting non-UniProt fastas with empty strings.
- Multiple FASTA paths in one call: iterate in order, then dedupe
  on `fasta.id` (first wins).
- Use `pandas.DataFrame` directly; no pydantic validation of the
  bulk frame (would be slow for large fastas; the data here is
  read-only and structural).
- `FastaRecord` is a small `@dataclass(frozen=True, slots=True)`:
  `(header: str, sequence: str)`. No biology semantics there — the
  per-protein interpretation happens in `annotation.py`.

## Tests

`tests/test_fasta_annotation.py`:

1. `test_parser_reads_uniprot_records` — feed a 2-record string,
   assert headers and sequences are correctly joined across wrapped
   lines.
2. `test_extract_gene_name_uniprot` — `extract_gene_name` matches
   the regex on `... GN=ygdT PE=4 SV=1` → `"ygdT"`; returns `""`
   for non-UniProt headers.
3. `test_count_tryptic_peptides_matches_prolfquapp` — use
   prolfquapp's docstring example (`"MKGLPRAKSHGSTGWGKRKRNKPK"`
   with `min_length=5`), assert the count matches what the R code
   would return. Confirm via the cleavage-site rule (K\|R not
   followed by P).
4. `test_fasta_to_dataframe_columns` — pass the prolfquapp fixture
   string (the one in `.getSequences()` of
   `get_annot_from_FASTA.R`) and assert the column set, row count
   after decoy removal (4 REV_ entries excluded → expect N-4 rows),
   and the per-row values for the first record.
5. `test_fasta_to_dataframe_filters_decoys` — `decoy_pattern=""`
   keeps all rows; default keeps only forward.
6. `test_fasta_to_dataframe_non_uniprot` — `is_uniprot=False`,
   `proteinname` equals `fasta.id`.
7. `test_fasta_to_dataframe_gene_name_gating` — a 1-record UniProt
   fasta produces no `gene_name` column (matches prolfquapp's
   "added only if >1 matches" behaviour).
8. `test_fasta_to_dataframe_include_sequence` — when
   `include_sequence=True`, the `sequence` column is present.

Fixture: a small text-string FASTA literal embedded in the test
file. No external file required.

## Out of scope

- Joining the annotation frame into `adata.var` or stuffing it in
  `adata.uns`. That belongs to a separate "annotate AnnData with
  FASTA" wiring step; this plan delivers just the parser + frame
  builder.
- The `ProteinAnnotation` R6 wrapper class (contamination /
  decoy patterns at runtime, peptide-count joining onto LFQData,
  filtering API). Not needed for the data-extraction goal.
- Multi-fasta concordance checks. Just dedupe on `fasta.id`.
