# APB Documentation

APB converts proteomics vendor output into AnnData and MuData through declarative
parsing rules. The package is a library plus the `apb` CLI; it ships no GUI.

## Where to go next

| If you want to... | Read |
|-------------------|------|
| Understand the package layout and current supported vendors | [Package Architecture](ARCHITECTURE.md) |
| Understand how vendor tables become AnnData or MuData | [Parsing Architecture](parsing_architecture.md) |
| Write or review parsing-rule JSON | [JSON Schema](json_schema.md) |
| Understand search-parameter parsing and rule version selection | [Parameter Parsers](parameter_parsers.md) |

## Common commands

```bash
apb list
apb validate
apb convert report.tsv --params report.log.txt
apb annotate data.h5mu module_settings.toml
apb fasta data.h5mu proteome.fasta
```

## Build these docs

```bash
make docs
```

The static site is written to `public/`; open `public/index.html` to browse it
from disk.
