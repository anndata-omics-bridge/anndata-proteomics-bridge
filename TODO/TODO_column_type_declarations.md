# Declared column types in the parsing-rule JSON

**Status:** proposed; not started
**Date:** 2026-07-31
**Scope:** `apb/` — rule schema, both generated JSON schemas, all packaged rule documents, converters
**Origin:** owner's design, prompted by the 2026-07-31 Spectronaut comma-decimal incident (`CHANGES.md`)

## Why

`converters/_axis.build_axis_frame` applies **no typing at all**. Every `var` and `obs` dtype in an
APB object is whatever pandas happened to infer while reading the vendor file. Nothing in a rule
states what a column is supposed to be, so nothing can notice when the inference is wrong.

That is how the comma-decimal Spectronaut exports stayed silent: `PG.Quantity` arrived as a string
column, coerced to NaN, and reached the protein AnnData with 15 of 66 222 values. The layer-occupancy
check added the same day catches the *symptom* at conversion time. The durable fix is for the rule to
declare the intended type, so a column failing to reach it is an attributable contract failure rather
than an invisible inference miss.

## Design

Convert **only** where a type is declared; leave everything else exactly as read.

| Axis | Default | Declaring otherwise |
| --- | --- | --- |
| `layers` | numeric | already expressed by `encoding_mode = "factor"` + `categories` — unchanged |
| `columns.var` | **string** | `type = "double"` / `"int32"` / … ; converted only then |
| `columns.obs` | **string** (or factor) | same — converted only on an explicit type |

The layer mechanism is the precedent to follow: `Layer.encoding_mode` already declares "this column
is not a bare number, and here is how to represent it". This extends the same idea to the two axis
column groups, which currently have no equivalent.

`columns.{obs,var}.select` therefore grows from `Name = "Vendor col"` to a form that can carry a
type. Keep the plain-string form as shorthand for "source only, no conversion", so existing
documents stay readable and only columns that genuinely need a type gain one.

## Migration risks

- **Defaulting `var`/`obs` to string changes current output.** Any column that is numeric today
  *only* because pandas inferred it becomes a string until its rule declares a type. Audit these
  consumers before flipping the default:
  - `converters/assemble._format_charge`
  - the `proforma_ion` and `proforma_fragment` computes
  - ProteoBench role resolution (`proteobench/resolve.py`)
  - the FASTA join keys (`annotation/var_fasta._var_join_keys`)
- **The ProteoBench legacy intermediate hash is pinned** at `9077847f…`
  (`tests/test_proteobench.py:47`). It must not move. If it does, the typing change reached scoring
  input and the cause needs finding before proceeding.
- **All 12 packaged rule documents and both generated JSON schemas** are rewritten by any change to
  `columns.*.select`. Worth landing in one pass.

## Decisions

**Explicit annotation is the mechanism.** There is no alternative that is not inference in disguise
— pandas dtype guessing, content sniffing, or deriving a type from whichever consumer reads the
column later are all the same failure mode. The type is stated in the rule or it is guessed.

Only the **exceptions** are annotated. String is the default and is written nowhere; a column carries
a type only when it must be converted (charge, m/z, retention time, scores, q-values). Sequences,
accessions, protein groups, run names and modification strings stay silent. This mirrors `layers`,
where numeric is the default and only `factor` is declared.

**The output flip is a correction, not a regression.** The affected columns are numeric today only by
accident of inference, and that accident is what produced the all-null ProteoBench scores. Of the
consumers audited: ProteoBench resolves `PG_ProteinAccessions`, the FASTA join keys are accessions,
and the ProForma computes concatenate into strings — none need a numeric var column. Charge is the
real case and gets a declared int; read `converters/assemble._format_charge` first, since it exists
because pandas hands it a float today and should get simpler, not harder. `9077847f…` is the gate: if
it moves, stop and find the cause before continuing.

**One pass, not incremental.** No half-migrated state preserves current behaviour, so incremental
means a compatibility shim for an unreleased layout — ruled out by `CODE_REVIEW_Jul30.md`. Single
commit: schema, both generated JSON schemas, types in all 12 documents where numeric is relied on,
converters applying declared types only. Land this *before* the alias/`required` work so
`columns.*.select` migrates once.

## Related

- [TODO_column_extraction_required_columns_alphapep_qpx_2026-07-07.md](TODO_column_extraction_required_columns_alphapep_qpx_2026-07-07.md)
  also plans to change `columns.*.select` (alias lists, `required` on axis columns). Different
  origin, but the same field — sequence the two so that surface is migrated once, not twice.
- `docs/json_schema.md` — the current `columns.*.select` contract.
