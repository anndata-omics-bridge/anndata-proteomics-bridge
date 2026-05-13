.PHONY: help catalog select download fasta clean

.DEFAULT_GOAL := help

help:  ## Show this help message.
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-12s\033[0m %s\n", $$1, $$2}'

catalog: raw_file_db_full.csv  ## Build the full catalog CSV by downloading datapoint JSONs for all 8 quant modules.

raw_file_db_full.csv:
	uv run python extract_raw_file_db.py catalog \
		--json-dir ./json_dir \
		--output ./raw_file_db_full.csv

select: raw_file_db_selected.csv  ## Reduce the catalog to one row per (module, software, version) -> raw_file_db_selected.csv.

raw_file_db_selected.csv: raw_file_db_full.csv
	uv run python extract_raw_file_db.py select \
		--input ./raw_file_db_full.csv \
		--output ./raw_file_db_selected.csv

download: raw_file_db_selected.csv  ## Download raw quantification files for the selected rows -> raw_file_db_downloaded.csv.
	uv run python extract_raw_file_db.py download \
		--input ./raw_file_db_selected.csv \
		--json-dir ./json_dir \
		--output ./raw_file_db_downloaded.csv

fasta: ## Download the ProteoBench FASTA databases (HYE + HY mixed-species) to ./fasta/.
	mkdir -p ./fasta
	cd ./fasta && curl -sLO https://proteobench.cubimed.rub.de/fasta/ProteoBenchFASTA_MixedSpecies_HYE.zip
	cd ./fasta && curl -sLO https://proteobench.cubimed.rub.de/fasta/ProteoBenchFASTA_MixedSpecies_HY.zip
	cd ./fasta && unzip -q -o ProteoBenchFASTA_MixedSpecies_HYE.zip
	cd ./fasta && unzip -q -o ProteoBenchFASTA_MixedSpecies_HY.zip
	rm -rf ./fasta/__MACOSX

clean:  ## Remove downloaded JSON folders, FASTAs, and all CSVs.
	rm -rf ./json_dir ./fasta ./raw_file_db_full.csv ./raw_file_db_selected.csv ./raw_file_db_downloaded.csv
