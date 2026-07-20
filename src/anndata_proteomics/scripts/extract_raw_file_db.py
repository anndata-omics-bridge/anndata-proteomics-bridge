"""Build a local test-data cache from ProteoBench submissions.

Run the commands in order: ``catalog`` collects submission metadata, ``select``
chooses representative submissions without downloading them, and ``download``
fetches the selected vendor files. All generated files default to APB's
``test_data_download`` directory.
"""

import io
import json
import os
import shutil
import zipfile
from contextlib import chdir
from pathlib import Path
from typing import Literal

import pandas as pd
import requests
from bs4 import BeautifulSoup
from cyclopts import App

from anndata_proteomics.test_data import TEST_DATA_DIR


CONFIGS = {
    "dda_qexactive": {
        "repo_url": "https://github.com/Proteobench/Results_quant_ion_DDA/archive/refs/heads/main.zip"
    },
    "dda_astral": {
        "repo_url": "https://github.com/Proteobench/Results_quant_ion_DDA_Astral/archive/refs/heads/main.zip"
    },
    "dda_peptidoform": {
        "repo_url": "https://github.com/Proteobench/Results_quant_peptidoform_DDA/archive/refs/heads/main.zip"
    },
    "dia_astral": {
        "repo_url": "https://github.com/Proteobench/Results_quant_ion_DIA_Astral/archive/refs/heads/main.zip"
    },
    "dia_diapasef": {
        "repo_url": "https://github.com/Proteobench/Results_quant_ion_DIA_diaPASEF/archive/refs/heads/main.zip"
    },
    "dia_aif": {
        "repo_url": "https://github.com/Proteobench/Results_quant_ion_DIA_AIF/archive/refs/heads/main.zip"
    },
    "dia_zenotof": {
        "repo_url": "https://github.com/Proteobench/Results_quant_ion_DIA_ZenoTOF/archive/refs/heads/main.zip"
    },
    "dia_singlecell": {
        "repo_url": "https://github.com/Proteobench/Results_quant_ion_DIA_singlecell/archive/refs/heads/main.zip"
    },
}

ModuleKey = Literal[
    "dda_qexactive",
    "dda_astral",
    "dda_peptidoform",
    "dia_astral",
    "dia_diapasef",
    "dia_aif",
    "dia_zenotof",
    "dia_singlecell",
]
SelectionStrategy = Literal[
    "smallest-per-software-version",
    "smallest-per-software",
    "smallest-per-module",
    "all",
]

app = App(name="apb-testdata", help=__doc__, help_on_error=True)

GENERATED_CSV_NAMES = (
    "raw_file_db_full.csv",
    "raw_file_db_selected.csv",
    "raw_file_db_downloaded.csv",
)


def _feature_count(data: dict) -> int | float | None:
    """Return ProteoBench's canonical feature count from one submission.

    ProteoBench renamed the serialized ``nr_prec`` field to ``nr_feature`` so
    the metric also applies to non-precursor quantification modules. Older
    result repositories can still contain the legacy field.
    """
    value = data.get("nr_feature")
    return value if value is not None else data.get("nr_prec")


# --- ProteoBench download primitives -----------------------------------------
# Ported from `proteobench.utils.server_io` so this catalog tool needs no
# proteobench install. Both helpers depend only on requests + beautifulsoup4,
# which APB already requires (see pyproject.toml). Behaviour matches upstream
# `get_merged_json` (GitHub results-repo ZIP -> extract) and `get_raw_data`
# (scrape + download raw files from the ProteoBench datasets server) as of
# 2026-07; the unused DataFrame return of the original get_merged_json and its
# unused content-length read were dropped.

DATASETS_BASE_URL = "https://proteobench.cubimed.rub.de/datasets/"


def get_merged_json(repo_url: str) -> Path:
    """Download a results-repo ZIP from GitHub and extract it into the cwd.

    Returns the archive's actual extracted root. GitHub can redirect a historical
    repository URL to a renamed repository, so the ZIP root is discovered from its
    members instead of being derived from the request URL.
    """
    response = requests.get(repo_url)
    response.raise_for_status()
    output_directory = Path(repo_url.split("/")[-5])
    with zipfile.ZipFile(io.BytesIO(response.content)) as zip_ref:
        archive_roots = {
            Path(member.filename).parts[0]
            for member in zip_ref.infolist()
            if member.filename and Path(member.filename).parts
        }
        if len(archive_roots) != 1:
            raise RuntimeError(
                f"Expected one root folder in {repo_url}, found {sorted(archive_roots)}"
            )
        zip_ref.extractall(output_directory)
    return output_directory / archive_roots.pop()


def get_raw_data(
    df: pd.DataFrame,
    base_url: str = DATASETS_BASE_URL,
    output_directory: str = "extracted_files",
) -> dict[str, str]:
    """Download raw quantification files for the submissions listed in `df`.

    Scrapes the datasets-server directory listing, matches folder names to
    `df["intermediate_hash"]`, downloads each matching folder's ZIP(s), and
    extracts them to `{output_directory}/{hash}/`. Returns
    `{intermediate_hash: extract_dir}` for the folders found. Folders already
    present and non-empty are skipped (idempotent re-runs).
    """
    hash_vis_dir: dict[str, str] = {}
    hash_list = df["intermediate_hash"].tolist()

    response = requests.get(base_url)
    response.raise_for_status()
    soup = BeautifulSoup(response.text, "html.parser")
    folder_links = [
        link["href"].strip("/") for link in soup.find_all("a") if link["href"].endswith("/")
    ]
    matching_folders = [folder for folder in folder_links if folder in hash_list]

    for folder in matching_folders:
        extract_dir = f"{output_directory}/{folder}"
        if os.path.exists(extract_dir) and os.listdir(extract_dir):
            print(f"Folder already exists and is not empty, skipping download: {extract_dir}")
            hash_vis_dir[folder] = extract_dir
            continue

        folder_url = f"{base_url}{folder}/"
        print(f"Processing folder: {folder_url}")

        folder_response = requests.get(folder_url)
        folder_response.raise_for_status()
        folder_soup = BeautifulSoup(folder_response.text, "html.parser")
        zip_files = [
            link["href"] for link in folder_soup.find_all("a") if link["href"].endswith(".zip")
        ]

        for zip_file in zip_files:
            zip_url = f"{folder_url}{zip_file}"
            print(f"Downloading: {zip_url}")

            zip_response = requests.get(zip_url, stream=True)
            zip_response.raise_for_status()

            zip_filename = os.path.basename(zip_file)
            with open(zip_filename, "wb") as f:
                for data in zip_response.iter_content(1024):
                    f.write(data)

            os.makedirs(extract_dir, exist_ok=True)
            with zipfile.ZipFile(zip_filename, "r") as zip_ref:
                zip_ref.extractall(extract_dir)
                print(f"Extracted contents to: {extract_dir}")

            os.remove(zip_filename)
            hash_vis_dir[folder] = extract_dir

    return hash_vis_dir


def _download_module_jsons(repo_url: str, json_dir_root: Path) -> Path:
    """Download JSONs for one repo into `{json_dir_root}/{repo_name}/{repo_name}-main/`.

    Wipes any pre-existing `{json_dir_root}/{repo_name}` first so the download is
    authoritative. Returns the folder containing the *.json files.
    """
    repo_name = repo_url.split("/")[-5]
    target_parent = json_dir_root / repo_name
    target_json_dir = target_parent / f"{repo_name}-main"

    if target_parent.exists():
        shutil.rmtree(target_parent)

    json_dir_root.mkdir(parents=True, exist_ok=True)
    with chdir(json_dir_root):
        extracted_dir = get_merged_json(repo_url=repo_url)

    extracted_dir = json_dir_root / extracted_dir
    if extracted_dir != target_json_dir:
        extracted_dir.rename(target_json_dir)

    if not target_json_dir.is_dir():
        raise RuntimeError(f"Expected extracted folder not found: {target_json_dir}")
    return target_json_dir


@app.command
def catalog(
    *,
    catalog_csv: Path = TEST_DATA_DIR / "raw_file_db_full.csv",
    cache_dir: Path = TEST_DATA_DIR / "json_dir",
) -> None:
    """Refresh the complete ProteoBench submission catalog.

    Args:
        catalog_csv: CSV file to write with one row per ProteoBench submission.
        cache_dir: Directory for repository metadata and downloaded submission files.
    """
    cache_dir = cache_dir.resolve()

    rows: list[dict] = []
    for module_key in CONFIGS:
        repo_url = CONFIGS[module_key]["repo_url"]
        repo_name = repo_url.split("/")[-5]
        print(f"[{module_key}] downloading {repo_url}")
        metadata_dir = _download_module_jsons(repo_url, cache_dir)

        json_files = sorted(metadata_dir.glob("*.json"))
        print(f"[{module_key}] read {len(json_files)} JSON file(s) from {metadata_dir}")

        for jf in json_files:
            try:
                with open(jf, "r", encoding="utf-8") as f:
                    data = json.load(f)
            except json.JSONDecodeError as e:
                print(f"  skip {jf.name}: {e}")
                continue

            rows.append(
                {
                    "module": module_key,
                    "repo_name": repo_name,
                    "intermediate_hash": data.get("intermediate_hash", jf.stem),
                    "software_name": data.get("software_name", ""),
                    "software_version": data.get("software_version", ""),
                    "nr_feature": _feature_count(data),
                    "is_temporary": data.get("is_temporary"),
                    "old_new": data.get("old_new", ""),
                }
            )

    df = pd.DataFrame(rows)
    df.to_csv(catalog_csv, index=False)

    print(f"\nTotal rows: {len(df)}")
    print("\nRows per module:")
    print(df["module"].value_counts().to_string())
    print("\nUnique (software_name, software_version) per module:")
    unique_counts = df.groupby("module").apply(
        lambda g: g[["software_name", "software_version"]].drop_duplicates().shape[0],
        include_groups=False,
    )
    print(unique_counts.to_string())
    print(f"\nWritten to {catalog_csv}")


@app.command
def select(
    *,
    catalog_csv: Path = TEST_DATA_DIR / "raw_file_db_full.csv",
    selection_csv: Path = TEST_DATA_DIR / "raw_file_db_selected.csv",
    strategy: SelectionStrategy = "smallest-per-software-version",
    module: ModuleKey | None = None,
) -> None:
    """Select ProteoBench submissions for download.

    Use --strategy to choose how rows are reduced and --module to restrict the
    selection to one ProteoBench module. The default keeps the smallest submission
    for every module, software, and software version.

    Args:
        catalog_csv: CSV containing all cataloged ProteoBench submissions.
        selection_csv: CSV to write with the submissions selected for download.
        strategy: Row-selection policy; smallest means the lowest nr_feature.
        module: Restrict selection to this ProteoBench module.
    """
    df = pd.read_csv(catalog_csv)
    n_input = len(df)

    print("Available catalog rows per module:")
    print(df["module"].value_counts().to_string())

    if module is not None:
        df = df[df["module"] == module]

    n_considered = len(df)
    selected, n_missing, rule = _select_rows(df, strategy)
    print(f"\nSelection rule: {rule}")
    if module is not None:
        print(f"Module filter:   {module}")

    selected.to_csv(selection_csv, index=False)

    print(f"Catalog rows:         {n_input}")
    print(f"Rows considered:      {n_considered}")
    print(f"Dropped (nr_feature NA): {n_missing}")
    print(f"Selected rows:        {len(selected)}")
    print("\nSelected rows per module:")
    print(selected["module"].value_counts().to_string())
    print(f"\nWritten to {selection_csv}")


def _select_rows(
    df: pd.DataFrame,
    strategy: SelectionStrategy,
) -> tuple[pd.DataFrame, int, str]:
    """Apply a download-selection strategy to catalog rows."""
    if strategy == "all":
        return df.reset_index(drop=True), 0, "all catalog rows"

    group_columns = {
        "smallest-per-software-version": ["module", "software_name", "software_version"],
        "smallest-per-software": ["module", "software_name"],
        "smallest-per-module": ["module"],
    }[strategy]
    group_label = " + ".join(column.removesuffix("_name") for column in group_columns)
    n_missing = int(df["nr_feature"].isna().sum())
    ranked = df.dropna(subset=["nr_feature"]).sort_values(
        ["nr_feature", "intermediate_hash"], kind="stable"
    )
    selected = (
        ranked.groupby(group_columns, dropna=False, as_index=False).head(1).reset_index(drop=True)
    )
    rule = (
        f"smallest nr_feature per {group_label}; "
        "ties use the lexicographically smallest intermediate_hash"
    )
    return selected, n_missing, rule


def get_datasets_to_download(df: pd.DataFrame, output_directory: Path) -> tuple[pd.DataFrame, dict]:
    """Check which datasets are already present in `output_directory`.

    Copied from `benchmark_analysis.py` — same idempotency logic.
    Returns (rows_still_to_download, {hash: extracted_dir} for rows already present).
    """
    hash_list = df["intermediate_hash"].tolist()
    existing_hashes = []
    hash_vis_dir: dict[str, str] = {}

    if output_directory.exists():
        for hash_dir in os.listdir(output_directory):
            if hash_dir in hash_list:
                existing_hashes.append(hash_dir)
                hash_vis_dir[hash_dir] = str(output_directory / hash_dir)

    if existing_hashes:
        df_to_download = df[~df["intermediate_hash"].isin(existing_hashes)]
        return df_to_download, hash_vis_dir
    return df, hash_vis_dir


@app.command
def download(
    *,
    selection_csv: Path = TEST_DATA_DIR / "raw_file_db_selected.csv",
    cache_dir: Path = TEST_DATA_DIR / "json_dir",
    manifest_csv: Path = TEST_DATA_DIR / "raw_file_db_downloaded.csv",
) -> None:
    """Download raw quantification files listed in a selection CSV.

    Run ``catalog`` and ``select`` first to create the selection CSV.

    Args:
        selection_csv: CSV listing the submissions to download.
        cache_dir: Directory for repository metadata and downloaded submission files.
        manifest_csv: CSV to write with download paths, sizes, and statuses.
    """
    df = pd.read_csv(selection_csv)
    required = {"module", "repo_name", "intermediate_hash"}
    missing_cols = required - set(df.columns)
    if missing_cols:
        raise SystemExit(f"Input CSV missing required columns: {sorted(missing_cols)}")

    cache_dir = cache_dir.resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)

    hash_to_dir: dict[str, str] = {}
    for repo_name, group in df.groupby("repo_name"):
        module_output_dir = cache_dir / repo_name
        module_output_dir.mkdir(parents=True, exist_ok=True)

        df_to_download, already_present = get_datasets_to_download(group, module_output_dir)
        hash_to_dir.update(already_present)

        if len(df_to_download) > 0:
            print(
                f"[{repo_name}] downloading {len(df_to_download)} dataset(s) (already present: {len(already_present)})"
            )
            new_dirs = get_raw_data(df_to_download, output_directory=str(module_output_dir))
            hash_to_dir.update(new_dirs)
        else:
            print(f"[{repo_name}] all {len(group)} dataset(s) already present — skipping")

    out_rows = []
    for row in df.to_dict(orient="records"):
        h = row["intermediate_hash"]
        extract_dir = hash_to_dir.get(h)
        if extract_dir is None:
            row["input_file_path"] = ""
            row["input_file_size_bytes"] = None
            row["status"] = "not_on_server"
        else:
            inputs = sorted(Path(extract_dir).glob("input_file.*"))
            if inputs:
                inp = inputs[0]
                row["input_file_path"] = str(inp.relative_to(cache_dir))
                row["input_file_size_bytes"] = inp.stat().st_size
                row["status"] = "ok"
            else:
                row["input_file_path"] = ""
                row["input_file_size_bytes"] = None
                row["status"] = "input_file_missing"
        out_rows.append(row)

    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(manifest_csv, index=False)

    print(f"\nTotal rows:         {len(out_df)}")
    print("\nStatus breakdown:")
    print(out_df["status"].value_counts().to_string())
    print(f"\nWritten to {manifest_csv}")


@app.command
def clean(*, data_dir: Path = TEST_DATA_DIR) -> None:
    """Remove generated test-data artifacts from one explicit data directory.

    Args:
        data_dir: Root containing the generated cache, FASTAs, catalogs, and manifests.
    """
    removed = clean_generated_data(data_dir.resolve())
    print(f"Removed {len(removed)} generated path(s).")
    for path in removed:
        print(f"  {path}")


def clean_generated_data(test_data_dir: Path) -> list[Path]:
    """Remove only known generated test-data artifacts below test_data_dir."""
    test_data_dir = test_data_dir.expanduser().resolve()
    if test_data_dir in {Path(test_data_dir.anchor), Path.home().resolve()}:
        raise ValueError("Refusing to clean the filesystem or home root.")
    targets = [test_data_dir / "json_dir", test_data_dir / "fasta"]
    targets.extend(test_data_dir / name for name in GENERATED_CSV_NAMES)
    removed = []
    for path in targets:
        if path.is_dir():
            shutil.rmtree(path)
            removed.append(path)
        elif path.exists():
            path.unlink()
            removed.append(path)
    return removed


def build_database(results_dir: Path) -> pd.DataFrame:
    """Legacy: walk `results_dir` for locally-present module/hash folders and pair
    each intermediate_hash with its input_file.txt. Kept for the future `download`
    subcommand — not used by `catalog`.
    """
    rows = []

    for mod_dir in sorted(results_dir.iterdir()):
        if not mod_dir.is_dir():
            continue
        module = mod_dir.name

        main_dirs = list(mod_dir.glob("*-main"))
        if not main_dirs:
            continue
        json_dir = main_dirs[0]

        for hash_dir in sorted(mod_dir.iterdir()):
            if not hash_dir.is_dir() or hash_dir.name.endswith("-main"):
                continue

            intermediate_hash = hash_dir.name
            inp = hash_dir / "input_file.txt"
            jp = json_dir / f"{intermediate_hash}.json"

            if not jp.exists():
                rows.append(
                    {
                        "module": module,
                        "intermediate_hash": intermediate_hash,
                        "software_name": "",
                        "software_version": "",
                        "input_file_path": "",
                        "input_file_size_bytes": None,
                        "status": "json_missing",
                    }
                )
                continue

            meta = json.load(open(jp))
            software = meta.get("software_name", "")
            version = meta.get("software_version", "")

            if inp.exists():
                rel_path = inp.relative_to(results_dir)
                size = inp.stat().st_size
                status = "ok"
            else:
                rel_path = ""
                size = None
                status = "input_file_missing"

            rows.append(
                {
                    "module": module,
                    "intermediate_hash": intermediate_hash,
                    "software_name": software,
                    "software_version": version,
                    "input_file_path": str(rel_path),
                    "input_file_size_bytes": size,
                    "status": status,
                }
            )

    return pd.DataFrame(rows)


def main() -> None:
    """Run the test-data command-line application."""
    app()


if __name__ == "__main__":
    main()
