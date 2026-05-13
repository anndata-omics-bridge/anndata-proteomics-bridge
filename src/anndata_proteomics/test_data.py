"""Locate downloaded benchmark inputs for a packaged ParseRule by software_name.

The cache lives at `<repo_root>/test_data_download/json_dir/...` and is indexed
by `<repo_root>/test_data_download/raw_file_db_downloaded.csv`. The cache is
gitignored — regenerate via `test_data_download/Makefile`. Both the test suite
and `tools/generate_report.py` use this lookup.
"""

from __future__ import annotations

import csv
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
TEST_DATA_DIR = REPO_ROOT / "test_data_download"
DOWNLOADED_DB = TEST_DATA_DIR / "raw_file_db_downloaded.csv"
FASTA_DIR = TEST_DATA_DIR / "fasta"

# Which ProteoBench FASTA pairs with each module. Each module ships a single
# species mix; see ProteoBench's per-module docs under
# docs/available-modules/.
_MODULE_FASTA: dict[str, str] = {
    "dda_qexactive": "ProteoBenchFASTA_MixedSpecies_HYE.fasta",
    "dda_astral": "ProteoBenchFASTA_MixedSpecies_HYE.fasta",
    "dda_peptidoform": "ProteoBenchFASTA_MixedSpecies_HYE.fasta",
    "dia_aif": "ProteoBenchFASTA_MixedSpecies_HYE.fasta",
    "dia_astral": "ProteoBenchFASTA_MixedSpecies_HYE.fasta",
    "dia_diapasef": "ProteoBenchFASTA_MixedSpecies_HYE.fasta",
    "dia_zenotof": "ProteoBenchFASTA_MixedSpecies_HYE.fasta",
    # Single-cell DIA uses the HY (no E. coli) FASTA. The file inside the
    # zip is misleadingly named *_DDAQuantification_noecoli.fasta — that's
    # the HY mix.
    "dia_singlecell": "ProteoBenchFASTA_DDAQuantification_noecoli.fasta",
}

# ProteoBench checks in a representative parameter file per tool. Same workspace
# layout used by the tests in `tests/test_params_*.py`.
PROTEOBENCH_PARAMS_DIR = REPO_ROOT.parent / "ProteoBench" / "test" / "params"

# One canonical sample per packaged tool, keyed by the rule's software_name.
_PROTEOBENCH_PARAM_FIXTURES: dict[str, str] = {
    "DIA-NN": "DIANN_output_20240229_report.log.txt",
    "FragPipe": "fragpipe.workflow",
    "MaxQuant": "mqpar_mq2.6.2.0_1mc_MBR.xml",
    "PEAKS": "PEAKS_parameters_DDA.txt",
    "Spectronaut": "Spectronaut_dynamic.txt",
    "WOMBAT": "wombat_params.yaml",
}


def find_test_data(software_name: str) -> Path | None:
    """Return the first cached input for `software_name`, or None if absent.

    Matches rows where `status == "ok"`. Returns None when the cache index
    file does not exist (cache not regenerated yet).
    """
    if not DOWNLOADED_DB.exists():
        return None
    with open(DOWNLOADED_DB) as f:
        for row in csv.DictReader(f):
            if row["software_name"] == software_name and row.get("status") == "ok":
                return TEST_DATA_DIR / "json_dir" / row["input_file_path"]
    return None


def find_fasta(
    *, dataset_dir: Path | None = None, module: str | None = None
) -> Path | None:
    """Return the cached ProteoBench FASTA for a module, or None if missing.

    Pass either ``module`` (the ProteoBench module key as it appears in
    the ``module`` column of ``raw_file_db_downloaded.csv`` — e.g.
    ``"dia_singlecell"``) or ``dataset_dir`` (a path inside the cached
    ``test_data_download/json_dir/...`` tree; the module is read from
    the index). When both are given, ``module`` wins.

    Returns the absolute path to the unzipped FASTA, or ``None`` when
    the FASTA cache has not been downloaded yet
    (``make -C test_data_download fasta``).
    """
    if module is None and dataset_dir is not None:
        module = _module_for_dataset(dataset_dir)
    if module is None:
        return None
    fasta_name = _MODULE_FASTA.get(module)
    if fasta_name is None:
        return None
    path = FASTA_DIR / fasta_name
    return path if path.exists() else None


def _module_for_dataset(dataset_dir: Path) -> str | None:
    """Look up the ``module`` column for a cached dataset path."""
    if not DOWNLOADED_DB.exists():
        return None
    target = str(dataset_dir.resolve())
    with open(DOWNLOADED_DB) as f:
        for row in csv.DictReader(f):
            cached = TEST_DATA_DIR / "json_dir" / row["input_file_path"]
            if (
                str(cached.resolve()) == target
                or str(cached.parent.resolve()) == target
            ):
                return row.get("module")
    return None


def find_param_file(software_name: str) -> Path | None:
    """Return a sample parameter file for ``software_name``, or None.

    Resolves against the workspace-local ProteoBench fixtures directory.
    Returns None when no fixture is registered for the tool or the file
    is missing on disk.
    """
    fixture = _PROTEOBENCH_PARAM_FIXTURES.get(software_name)
    if fixture is None:
        return None
    path = PROTEOBENCH_PARAMS_DIR / fixture
    return path if path.exists() else None
