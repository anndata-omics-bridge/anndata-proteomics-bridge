# APB — top-level developer commands. Run `make` (or `make help`) to list targets.
# apb is a pure library + `apb` CLI; the marimo GUIs live in the sibling apb_studio package.

.DEFAULT_GOAL := help
.PHONY: help clean sync test lint check check-full audit package docs docs-serve

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) \
		| awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-10s\033[0m %s\n", $$1, $$2}'

clean:  ## Remove generated AnnData/MuData outputs
	find . -type f \( -name '*.h5ad' -o -name '*.h5mu' \) -not -path './.venv/*' -delete

sync:  ## Install the frozen development and documentation environment
	uv sync --frozen --extra dev --group docs

test:  ## Run the test suite
	uv run --frozen --extra dev pytest -q

lint:  ## Ruff check src/ and tests/
	uv run --frozen --extra dev ruff check .

check:  ## Run the commit-stage quality gate
	uv run pre-commit run --hook-stage pre-commit --all-files

check-full:  ## Run the push-stage quality gate
	uv run pre-commit run --hook-stage pre-push --all-files

audit:  ## Audit locked dependencies
	uv run pre-commit run dependency-audit --hook-stage manual --all-files

package:  ## Build and inspect the distribution contract
	uv run --frozen --extra dev python scripts/package_smoke.py

docs:  ## Build the documentation site into public/
	uv run --frozen --group docs mkdocs build --strict

docs-serve:  ## Preview the documentation site at http://127.0.0.1:8000
	uv run --frozen --group docs mkdocs serve
