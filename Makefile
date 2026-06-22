# APB — top-level developer commands. Run `make` (or `make help`) to list targets.
# GUI targets use the `gui` extra (marimo + plotly); uv installs it on first run.

UI     := src/anndata_proteomics/scripts/ui_test_tool.py
VIEWER := src/anndata_proteomics/scripts/anndataview.py

.DEFAULT_GOAL := help
.PHONY: help ui ui-edit viewer test lint

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) \
		| awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-10s\033[0m %s\n", $$1, $$2}'

ui:  ## Launch the test-data browser GUI (marimo)
	uv run --extra gui marimo run $(UI)

ui-edit:  ## Open the test-data browser GUI in marimo's editor
	uv run --extra gui marimo edit $(UI)

viewer:  ## Open one .h5ad in the AnnData viewer:  make viewer FILE=path.h5ad
	uv run --extra gui marimo run $(VIEWER) -- $(FILE)

test:  ## Run the test suite
	uv run --extra dev pytest -q

lint:  ## Ruff check src/ and tests/
	uv run --extra dev ruff check src/ tests/
