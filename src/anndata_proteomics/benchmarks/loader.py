"""
Benchmark data loader for proteomics datasets.

Downloads test data from ProteoBench repositories for testing converters.
Based on proteobench.utils.server_io.
"""

import io
import json
import zipfile
from pathlib import Path
from typing import Optional

import pandas as pd
import requests
from bs4 import BeautifulSoup


# Available benchmark repositories
BENCHMARK_REPOS = {
    "dda_qexactive": {
        "name": "DDA Q-Exactive Ion",
        "repo_url": "https://github.com/Proteobench/Results_quant_ion_DDA/archive/refs/heads/main.zip",
        "data_url": "https://proteobench.cubimed.rub.de/datasets/",
    },
    "dda_astral": {
        "name": "DDA Astral Ion",
        "repo_url": "https://github.com/Proteobench/Results_quant_ion_DDA_Astral/archive/refs/heads/main.zip",
        "data_url": "https://proteobench.cubimed.rub.de/datasets/",
    },
    "dia_aif": {
        "name": "DIA AIF Ion",
        "repo_url": "https://github.com/Proteobench/Results_quant_ion_DIA_AIF/archive/refs/heads/main.zip",
        "data_url": "https://proteobench.cubimed.rub.de/datasets/",
    },
    "dia_astral": {
        "name": "DIA Astral Ion",
        "repo_url": "https://github.com/Proteobench/Results_quant_ion_DIA_Astral/archive/refs/heads/main.zip",
        "data_url": "https://proteobench.cubimed.rub.de/datasets/",
    },
    "dia_diapasef": {
        "name": "DIA diaPASEF Ion",
        "repo_url": "https://github.com/Proteobench/Results_quant_ion_DIA_diaPASEF/archive/refs/heads/main.zip",
        "data_url": "https://proteobench.cubimed.rub.de/datasets/",
    },
}


class BenchmarkDataLoader:
    """
    Downloads and manages ProteoBench test datasets.

    Example:
        loader = BenchmarkDataLoader("dda_qexactive", output_dir="test_data")
        loader.load_metadata()
        loader.list_software()
        df = loader.get_datasets(software_filter=["MaxQuant"], max_per_software=1)
    """

    def __init__(
        self,
        benchmark: str = "dda_qexactive",
        output_dir: str = "test_data",
    ):
        if benchmark not in BENCHMARK_REPOS:
            raise ValueError(f"Unknown benchmark: {benchmark}. Available: {list(BENCHMARK_REPOS.keys())}")

        self.benchmark = benchmark
        self.config = BENCHMARK_REPOS[benchmark]
        self.output_dir = Path(output_dir)
        self.benchmark_dir = self.output_dir / benchmark
        self.benchmark_dir.mkdir(parents=True, exist_ok=True)

        self._metadata: Optional[pd.DataFrame] = None
        self._hash_to_dir: dict[str, Path] = {}

    @property
    def repo_url(self) -> str:
        return self.config["repo_url"]

    @property
    def data_url(self) -> str:
        return self.config["data_url"]

    @property
    def repo_name(self) -> str:
        return self.repo_url.split("/")[-5]

    @property
    def metadata(self) -> pd.DataFrame:
        """Lazy-load metadata on first access."""
        if self._metadata is None:
            self.load_metadata()
        return self._metadata

    def load_metadata(self) -> pd.DataFrame:
        """Download and parse benchmark metadata from GitHub."""
        print(f"Downloading metadata from {self.repo_url}...")
        response = requests.get(self.repo_url)
        response.raise_for_status()

        with zipfile.ZipFile(io.BytesIO(response.content)) as zip_ref:
            zip_ref.extractall(self.benchmark_dir)

        json_data = self._parse_json_files()
        self._metadata = pd.json_normalize(json_data)
        print(f"Found {len(self._metadata)} benchmark submissions")
        return self._metadata

    def _parse_json_files(self) -> list[dict]:
        """Parse all JSON files from extracted repo."""
        base_path = self.benchmark_dir / f"{self.repo_name}-main"
        results = []

        for json_file in base_path.rglob("*.json"):
            try:
                with open(json_file, "r", encoding="utf-8") as f:
                    results.append(json.load(f))
            except json.JSONDecodeError as e:
                print(f"Error reading {json_file}: {e}")

        return results

    def list_software(self) -> dict[str, list[str]]:
        """List available software and versions."""
        result = {}
        for _, row in self.metadata.iterrows():
            sw, ver = row["software_name"], row["software_version"]
            result.setdefault(sw, [])
            if ver not in result[sw]:
                result[sw].append(ver)
        return result

    def get_datasets(
        self,
        software_filter: Optional[list[str]] = None,
        max_per_software: Optional[int] = 1,
    ) -> pd.DataFrame:
        """
        Get test datasets, downloading if needed.

        Args:
            software_filter: Only include these software tools
            max_per_software: Max datasets per software (None for all)

        Returns:
            DataFrame with local_dir and input_file paths added
        """
        df = self.metadata.copy()

        if software_filter:
            df = df[df["software_name"].isin(software_filter)]

        if max_per_software:
            df = df.groupby("software_name").head(max_per_software)

        # Build list of datasets with software info for organized folder structure
        datasets = [
            {
                "hash": row["intermediate_hash"],
                "software": row["software_name"],
                "version": str(row["software_version"]),
            }
            for _, row in df.iterrows()
        ]
        self._download_datasets(datasets)

        df["local_dir"] = df["intermediate_hash"].map(
            lambda h: str(self._hash_to_dir.get(h, ""))
        )
        df["input_file"] = df["local_dir"].apply(self._get_input_file_path)

        return df

    def _download_datasets(self, datasets: list[dict]) -> None:
        """Download datasets organized by software/version."""
        available = self._fetch_available_folders()
        to_download = [d for d in datasets if d["hash"] in available]

        print(f"Found {len(to_download)} datasets to process")

        for i, dataset in enumerate(to_download, 1):
            self._download_single(dataset, i, len(to_download))

    def _fetch_available_folders(self) -> set[str]:
        """Get list of available dataset folders from server."""
        response = requests.get(self.data_url)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, "html.parser")
        return {
            link["href"].strip("/")
            for link in soup.find_all("a")
            if link.get("href", "").endswith("/")
        }

    def _download_single(self, dataset: dict, current: int, total: int) -> None:
        """Download a single dataset to software/version folder."""
        hash_id = dataset["hash"]
        software = dataset["software"]
        version = dataset["version"]

        # Organize by benchmark/software/version
        extract_dir = self.benchmark_dir / software / version

        if extract_dir.exists() and any(extract_dir.iterdir()):
            print(f"[{current}/{total}] Already exists: {software}/{version}")
            self._hash_to_dir[hash_id] = extract_dir
            return

        print(f"[{current}/{total}] Downloading: {software}/{version}")
        folder_url = f"{self.data_url}{hash_id}/"

        for zip_url in self._find_zip_files(folder_url):
            self._download_and_extract_zip(zip_url, extract_dir)

        self._hash_to_dir[hash_id] = extract_dir

    def _find_zip_files(self, folder_url: str) -> list[str]:
        """Find all zip files in a folder."""
        response = requests.get(folder_url)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, "html.parser")
        return [
            f"{folder_url}{link['href']}"
            for link in soup.find_all("a")
            if link.get("href", "").endswith(".zip")
        ]

    def _download_and_extract_zip(self, zip_url: str, extract_dir: Path) -> None:
        """Download and extract a single zip file."""
        response = requests.get(zip_url, stream=True)
        response.raise_for_status()

        zip_data = io.BytesIO(response.content)
        extract_dir.mkdir(parents=True, exist_ok=True)

        with zipfile.ZipFile(zip_data) as zip_ref:
            zip_ref.extractall(extract_dir)

    @staticmethod
    def _get_input_file_path(local_dir: str) -> Optional[str]:
        """Get input file path if it exists."""
        if not local_dir:
            return None
        input_file = Path(local_dir) / "input_file.txt"
        return str(input_file) if input_file.exists() else None


def list_available_software(benchmark: str = "dda_qexactive") -> dict[str, list[str]]:
    """List available software and versions in a benchmark."""
    return BenchmarkDataLoader(benchmark).list_software()


def get_test_datasets(
    benchmark: str = "dda_qexactive",
    software_filter: Optional[list[str]] = None,
    output_dir: str = "test_data",
    max_per_software: int = 1,
) -> pd.DataFrame:
    """Get test datasets for specific software tools."""
    loader = BenchmarkDataLoader(benchmark, output_dir)
    return loader.get_datasets(software_filter, max_per_software)
