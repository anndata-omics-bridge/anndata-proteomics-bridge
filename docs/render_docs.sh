#!/usr/bin/env bash
# Build the MkDocs documentation site.
#
# Usage: docs/render_docs.sh
# Output: public/index.html
set -euo pipefail
cd "$(dirname "$0")/.."

uv run --group docs mkdocs build
echo "wrote public/index.html"
