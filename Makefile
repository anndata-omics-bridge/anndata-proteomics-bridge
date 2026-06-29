# APB — top-level developer commands. Run `make` (or `make help`) to list targets.
# apb is a pure library + `apb` CLI; the marimo GUIs live in the sibling apb_studio package.

.DEFAULT_GOAL := help
.PHONY: help clean test lint docs docs-serve

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) \
		| awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-10s\033[0m %s\n", $$1, $$2}'

clean:  ## Remove generated AnnData/MuData outputs
	find . -type f \( -name '*.h5ad' -o -name '*.h5mu' \) -not -path './.venv/*' -delete

test:  ## Run the test suite
	uv run --extra dev pytest -q

lint:  ## Ruff check src/ and tests/
	uv run --extra dev ruff check src/ tests/

docs:  ## Build the documentation site into public/
	uv run --group docs mkdocs build

docs-serve:  ## Preview the documentation site at http://127.0.0.1:8000
	uv run --group docs mkdocs serve
