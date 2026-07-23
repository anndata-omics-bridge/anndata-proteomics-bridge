# ProteoBench protein mapper attribution

`mapper.csv` is redistributed from ProteoBench tag `v0.17.0`, commit
`fc95e712ca0466485814d3895087a048cfc0d2b0`, where it is stored at
`proteobench/io/parsing/io_parse_settings/mapper.csv`.

ProteoBench is licensed under the Apache License 2.0. The mapper is retained to
reproduce ProteoBench's accession-to-description normalization without loading a
FASTA or depending on the ProteoBench Python package. The license text is
included beside this notice as `APACHE-2.0.txt`.

The compatibility-only raw modification renderer in `mapping.py` is adapted to
APB's feature-aligned representation from ProteoBench's
`proteobench/io/parsing/proforma.py` at the same revision. APB's canonical
ProForma identifiers are not changed.
