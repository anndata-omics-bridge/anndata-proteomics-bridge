"""Locate downloaded benchmark inputs for a packaged ParseRule by software_name.

The canonical cache lives at `<repo_root>/test_data_download/json_dir/...` and
is indexed by `<repo_root>/test_data_download/raw_file_db_downloaded.csv`. It is
gitignored, regenerated with ``apb-testdata``, and consumed by the test suite.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
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
    "AlphaDIA": "alphadia_1.10.3.log.txt",
    "AlphaPept": "alphapept_0.4.9.yaml",
    "DIA-NN": "DIANN_output_20240229_report.log.txt",
    "FragPipe": "fragpipe.workflow",
    "MaxQuant": "mqpar_mq2.6.2.0_1mc_MBR.xml",
    "PEAKS": "PEAKS_parameters_DDA.txt",
    "Sage": "sage_parameterfile.json",
    "Spectronaut": "Spectronaut_dynamic.txt",
    "WOMBAT": "wombat_params.yaml",
}


@dataclass(frozen=True, slots=True)
class VendorDataUnavailable:
    """A cached vendor input could not be located."""

    software_name: str


@dataclass(frozen=True, slots=True)
class FastaUnavailable:
    """A cached FASTA could not be located for the requested module."""

    module: str


@dataclass(frozen=True, slots=True)
class DatasetFastaUnavailable:
    """A cached FASTA could not be located for a requested dataset."""

    dataset_dir: Path


@dataclass(frozen=True, slots=True)
class AnnotationUnavailable:
    """A cached sample annotation could not be located for a module."""

    module: str


@dataclass(frozen=True, slots=True)
class DatasetModuleUnavailable:
    """A cached dataset path could not be associated with a module."""

    dataset_dir: Path


@dataclass(frozen=True, slots=True)
class ParameterFileUnavailable:
    """No registered parameter fixture exists for a software name."""

    software_name: str


def find_test_data(software_name: str) -> Path | VendorDataUnavailable:
    """Return the first cached input for ``software_name``.

    Matches rows where ``status == "ok"``. A typed unavailable result records
    that the cache index or matching artifact is absent.
    """
    if not DOWNLOADED_DB.exists():
        return VendorDataUnavailable(software_name)
    with open(DOWNLOADED_DB) as f:
        for row in csv.DictReader(f):
            if row["software_name"] == software_name and row.get("status") == "ok":
                return TEST_DATA_DIR / "json_dir" / row["input_file_path"]
    return VendorDataUnavailable(software_name)


def find_test_data_for_version(
    software_name: str,
    version_pattern: str,
) -> Path | VendorDataUnavailable:
    """Return the first cached input covered by one software-version pattern.

    This is deliberately distinct from :func:`find_test_data`: version-aware
    selection is a different lookup, not optional behavior on the same function.
    """
    if not DOWNLOADED_DB.exists():
        return VendorDataUnavailable(software_name)
    with open(DOWNLOADED_DB) as f:
        for row in csv.DictReader(f):
            if row["software_name"] != software_name or row.get("status") != "ok":
                continue
            if not software_version_matches(
                version_pattern,
                row["software_version"],
            ):
                continue
            return TEST_DATA_DIR / "json_dir" / row["input_file_path"]
    return VendorDataUnavailable(software_name)


def find_fasta_for_module(
    module: str,
    *,
    test_data_dir: Path = TEST_DATA_DIR,
) -> Path | FastaUnavailable:
    """Return the cached ProteoBench FASTA for one explicit module.

    ``test_data_dir`` selects an explicit cache root. A typed unavailable
    result records an unknown module or an absent downloaded FASTA.
    """
    test_data_dir = test_data_dir.expanduser().resolve()
    fasta_name = _MODULE_FASTA.get(module)
    if fasta_name is None:
        return FastaUnavailable(module)
    path = test_data_dir / "fasta" / fasta_name
    return path if path.exists() else FastaUnavailable(module)


def find_fasta_for_dataset(
    dataset_dir: Path,
    *,
    test_data_dir: Path = TEST_DATA_DIR,
) -> Path | DatasetFastaUnavailable:
    """Resolve a cached dataset to its module and then locate that module's FASTA."""
    module = _module_for_dataset(dataset_dir, test_data_dir=test_data_dir)
    if isinstance(module, DatasetModuleUnavailable):
        return DatasetFastaUnavailable(dataset_dir)
    fasta = find_fasta_for_module(module, test_data_dir=test_data_dir)
    return DatasetFastaUnavailable(dataset_dir) if isinstance(fasta, FastaUnavailable) else fasta


def find_annotation(
    *,
    module: str,
    test_data_dir: Path = TEST_DATA_DIR,
) -> Path | AnnotationUnavailable:
    """Return the downloaded ProteoBench module annotation TOML."""
    filename = _MODULE_ANNOTATION.get(module)
    if filename is None:
        return AnnotationUnavailable(module)
    path = test_data_dir.expanduser().resolve() / "annotations" / filename
    return path if path.is_file() else AnnotationUnavailable(module)


def _module_for_dataset(
    dataset_dir: Path,
    *,
    test_data_dir: Path = TEST_DATA_DIR,
) -> str | DatasetModuleUnavailable:
    """Look up the ``module`` column for a cached dataset path."""
    test_data_dir = test_data_dir.expanduser().resolve()
    downloaded_db = test_data_dir / "raw_file_db_downloaded.csv"
    if not downloaded_db.exists():
        return DatasetModuleUnavailable(dataset_dir)
    target = str(dataset_dir.resolve())
    with open(downloaded_db) as f:
        for row in csv.DictReader(f):
            cached = test_data_dir / "json_dir" / row["input_file_path"]
            if str(cached.resolve()) == target or str(cached.parent.resolve()) == target:
                module = row.get("module")
                if module:
                    return module
                return DatasetModuleUnavailable(dataset_dir)
    return DatasetModuleUnavailable(dataset_dir)


def find_param_file(software_name: str) -> Path | ParameterFileUnavailable:
    """Return a sample parameter file for ``software_name``.

    Resolves against the in-repo ``tests/params/`` fixtures directory. A typed
    unavailable result records an unregistered or missing fixture.
    """
    fixture = _PROTEOBENCH_PARAM_FIXTURES.get(software_name)
    if fixture is None:
        return ParameterFileUnavailable(software_name)
    path = PARAM_FIXTURE_DIR / fixture
    return path if path.exists() else ParameterFileUnavailable(software_name)
