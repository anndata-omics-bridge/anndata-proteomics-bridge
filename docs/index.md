# APB Documentation

APB converts proteomics vendor output into AnnData and MuData through declarative
parsing rules. The package is a library plus the `apb` CLI; GUI and workflow
orchestration live in the sibling `apb_studio` package.

## Where to go next

| If you want to... | Read |
|-------------------|------|
| Understand the package layout and current supported vendors | [Package Architecture](ARCHITECTURE.md) |
| Understand how vendor tables become AnnData or MuData | [Parsing Architecture](parsing_architecture.md) |
| Write or review parsing-rule TOMLs | [TOML Schema](toml_schema.md) |
| Understand search-parameter parsing and rule version selection | [Parameter Parsers](parameter_parsers.md) |

## Common commands

```bash
apb list
apb validate
apb convert report.tsv --params report.log.txt
apb annotate data.h5mu annotation.toml
apb fasta data.h5mu proteome.fasta
```

## Build these docs

```bash
docs/render_docs.sh
```

The static site is written to `public/`; open `public/index.html` to browse it
from disk.
