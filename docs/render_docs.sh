#!/usr/bin/env bash
# Render the parsing docs (Markdown -> standalone HTML) with pandoc.
# Mermaid diagrams render client-side via mermaid.js (needs internet to load it).
#
# Usage:  docs/render_docs.sh
# Requires: pandoc (brew install pandoc).
set -euo pipefail
cd "$(dirname "$0")"

for doc in parsing_architecture parameter_parsers; do
  pandoc -f gfm -s --toc --toc-depth=2 \
    --lua-filter=_md2html_links.lua \
    -H _pandoc_mermaid.html \
    -B _pandoc_nav.html \
    -V pagetitle="$doc" \
    -o "$doc.html" "$doc.md"
  echo "wrote $doc.html"
done
