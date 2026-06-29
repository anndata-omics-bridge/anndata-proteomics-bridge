# Documentation build targets. The MkDocs config lives at ../mkdocs.yml.
#
# Requires: uv.

.PHONY: all serve build clean

all: build

serve:
	cd .. && uv run --group docs mkdocs serve

build:
	cd .. && uv run --group docs mkdocs build

clean:
	rm -rf ../public
