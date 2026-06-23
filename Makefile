# APB — top-level developer commands. Run `make` (or `make help`) to list targets.
# GUI targets use the `gui` extra (marimo + plotly); uv installs it on first run.

UI     := src/anndata_proteomics/scripts/ui_test_tool.py
VIEWER := src/anndata_proteomics/scripts/anndataview.py
CONVERT_ONE := anndata_proteomics.scripts.convert_one
CONVERTED_DIR := logs/ui_converted

DIANN_INPUT := Results_quant_ion_DIA_Astral/4cc229a793e08041348ccff49d7afaf0117a9da0/input_file.parquet
DIANN_PARAMS := test_data_download/json_dir/Results_quant_ion_DIA_Astral/4cc229a793e08041348ccff49d7afaf0117a9da0/param_0..txt

.DEFAULT_GOAL := help
.PHONY: help ui ui-edit viewer spectro_mudata diann_mudata clean test lint

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) \
		| awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-10s\033[0m %s\n", $$1, $$2}'

ui:  ## Launch the test-data browser GUI (marimo)
	uv run --extra gui marimo run $(UI)

ui-edit:  ## Open the test-data browser GUI in marimo's editor
	uv run --extra gui marimo edit $(UI)

viewer:  ## Open one .h5ad in the AnnData viewer:  make viewer FILE=path.h5ad
	uv run --extra gui marimo run $(VIEWER) -- $(FILE)

spectro_mudata:  ## Convert all Spectronaut test files to MuData
	@set -e; \
	rows=$$(uv run python -c 'from anndata_proteomics.scripts import _ui_support as ui; cat = ui.load_catalog(); rows = cat[cat["slug"].eq("spectronaut") & cat["targets_str"].str.contains("mudata", regex=False)]; [print(f"{row.input_file_path}|{row.param_path}") for _, row in rows.iterrows()]'); \
	for row in $$rows; do \
		input="$${row%%|*}"; \
		params="$${row#*|}"; \
		outdir="$(CONVERTED_DIR)/$$(date +%Y%m%dT%H%M%S)_spectronaut_mudata"; \
		mkdir -p "$$outdir"; \
		cmd="uv run python -m $(CONVERT_ONE) --input $$input --slug spectronaut --target mudata --params $$params --outdir $$outdir"; \
		printf '$$ %s\n' "$$cmd" | tee "$$outdir/console.log"; \
		eval "$$cmd" 2>&1 | tee -a "$$outdir/console.log"; \
	done

diann_mudata:  ## Convert the default DIA-NN test file to MuData
	@outdir="$(CONVERTED_DIR)/$$(date +%Y%m%dT%H%M%S)_diann_mudata"; \
	mkdir -p "$$outdir"; \
	cmd="uv run python -m $(CONVERT_ONE) --input $(DIANN_INPUT) --slug diann --target mudata --params $(DIANN_PARAMS) --outdir $$outdir"; \
	printf '$$ %s\n' "$$cmd" | tee "$$outdir/console.log"; \
	eval "$$cmd" 2>&1 | tee -a "$$outdir/console.log"

clean:  ## Remove generated AnnData/MuData outputs and GUI conversion logs
	rm -rf $(CONVERTED_DIR)
	find . -type f \( -name '*.h5ad' -o -name '*.h5mu' \) -not -path './.venv/*' -delete
	mkdir -p $(CONVERTED_DIR)

test:  ## Run the test suite
	uv run --extra dev pytest -q

lint:  ## Ruff check src/ and tests/
	uv run --extra dev ruff check src/ tests/
