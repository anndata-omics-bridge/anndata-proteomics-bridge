"""Locate downloaded benchmark inputs for a packaged ParseRule by software_name.

The canonical cache lives at `<repo_root>/test_data_download/json_dir/...` and
is indexed by `<repo_root>/test_data_download/raw_file_db_downloaded.csv`. It is
gitignored, regenerated with ``apb-testdata``, and consumed by the test suite.
"""

from __future__ import annotations

import csv
from pathlib import Path

from anndata_proteomics.rules.loader import software_version_matches

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
TEST_DATA_DIR = REPO_ROOT / "test_data_download"
DOWNLOADED_DB = TEST_DATA_DIR / "raw_file_db_downloaded.csv"
FASTA_DIR = TEST_DATA_DIR / "fasta"
ANNOTATION_DIR = TEST_DATA_DIR / "annotations"

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

_MODULE_ANNOTATION: dict[str, str] = {
    "dda_qexactive": "dda_qexactive.toml",
    "dda_astral": "dda_astral.toml",
    "dda_peptidoform": "dda_peptidoform.toml",
    "dia_aif": "dia_aif.toml",
    "dia_astral": "dia_astral.toml",
    "dia_diapasef": "dia_diapasef.toml",
    "dia_zenotof": "dia_zenotof.toml",
    "dia_singlecell": "dia_singlecell.toml",
}

# A representative parameter file per tool, committed in-repo under tests/params/
# (the same fixtures the tests/test_params_*.py suite reads). Kept in-repo so the
# CLI integration test and the report generator need no external ProteoBench
# checkout.
PARAM_FIXTURE_DIR = REPO_ROOT / "tests" / "params"

# One canonical sample per packaged tool, keyed by the rule's software_name.
_PROTEOBENCH_PARAM_FIXTURES: dict[str, str] = {
    "AlphaPept": "alphapept_0.4.9.yaml",
    "DIA-NN": "DIANN_output_20240229_report.log.txt",
    "FragPipe": "fragpipe.workflow",
    "MaxQuant": "mqpar_mq2.6.2.0_1mc_MBR.xml",
    "PEAKS": "PEAKS_parameters_DDA.txt",
    "Sage": "sage_parameterfile.json",
    "Spectronaut": "Spectronaut_dynamic.txt",
    "WOMBAT": "wombat_params.yaml",
}


def find_test_data(software_name: str, version_pattern: str | None = None) -> Path | None:
    """Return the first cached input for `software_name`, or None if absent.

    Matches rows where `status == "ok"`. Returns None when the cache index
    file does not exist (cache not regenerated yet).

    Pass a rule's ``software_version`` as ``version_pattern`` to restrict the result to a
    submission that rule actually covers. Without it, a vendor whose cached submissions
    span several export schemas can return a file no packaged rule was written for — e.g.
    Sage's charge-collapsed 0.14.6 ``lfq.tsv``, which the charge-resolved ion rule cannot
    convert.
    """
    if not DOWNLOADED_DB.exists():
        return None
    with open(DOWNLOADED_DB) as f:
        for row in csv.DictReader(f):
            if row["software_name"] != software_name or row.get("status") != "ok":
                continue
            if version_pattern is not None and not software_version_matches(
                version_pattern,
                row["software_version"],
            ):
                continue
            return TEST_DATA_DIR / "json_dir" / row["input_file_path"]
    return None


def find_fasta(
    *,
    dataset_dir: Path | None = None,
    module: str | None = None,
    test_data_dir: Path = TEST_DATA_DIR,
) -> Path | None:
    """Return the cached ProteoBench FASTA for a module, or None if missing.

    Pass either ``module`` (the ProteoBench module key as it appears in
    the ``module`` column of ``raw_file_db_downloaded.csv`` — e.g.
    ``"dia_singlecell"``) or ``dataset_dir`` (a path inside the cached
    ``test_data_download/json_dir/...`` tree; the module is read from
    the index). ``test_data_dir`` selects an explicit cache root. When both
    ``module`` and ``dataset_dir`` are given, ``module`` wins.

    Returns the absolute path to the unzipped FASTA, or ``None`` when
    the FASTA cache has not been downloaded yet
    (``apb-testdata fasta``).
    """
    test_data_dir = test_data_dir.expanduser().resolve()
    if module is None and dataset_dir is not None:
        module = _module_for_dataset(dataset_dir, test_data_dir=test_data_dir)
    if module is None:
        return None
    fasta_name = _MODULE_FASTA.get(module)
    if fasta_name is None:
        return None
    path = test_data_dir / "fasta" / fasta_name
    return path if path.exists() else None


def find_annotation(
    *,
    module: str,
    test_data_dir: Path = TEST_DATA_DIR,
) -> Path | None:
    """Return the downloaded ProteoBench module annotation TOML, if available."""
    filename = _MODULE_ANNOTATION.get(module)
    if filename is None:
        return None
    path = test_data_dir.expanduser().resolve() / "annotations" / filename
    return path if path.is_file() else None


def _module_for_dataset(dataset_dir: Path, *, test_data_dir: Path = TEST_DATA_DIR) -> str | None:
    """Look up the ``module`` column for a cached dataset path."""
    test_data_dir = test_data_dir.expanduser().resolve()
    downloaded_db = test_data_dir / "raw_file_db_downloaded.csv"
    if not downloaded_db.exists():
        return None
    target = str(dataset_dir.resolve())
    with open(downloaded_db) as f:
        for row in csv.DictReader(f):
            cached = test_data_dir / "json_dir" / row["input_file_path"]
            if str(cached.resolve()) == target or str(cached.parent.resolve()) == target:
                return row.get("module")
    return None


def find_param_file(software_name: str) -> Path | None:
    """Return a sample parameter file for ``software_name``, or None.

    Resolves against the in-repo ``tests/params/`` fixtures directory. Returns
    None when no fixture is registered for the tool or the file is missing on
    disk.
    """
    fixture = _PROTEOBENCH_PARAM_FIXTURES.get(software_name)
    if fixture is None:
        return None
    path = PARAM_FIXTURE_DIR / fixture
    return path if path.exists() else None
