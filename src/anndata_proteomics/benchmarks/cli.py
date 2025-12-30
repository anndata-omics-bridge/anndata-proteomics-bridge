"""CLI for downloading ProteoBench test data."""

from pathlib import Path
from typing import Annotated, Optional

import cyclopts

from .loader import BENCHMARK_REPOS, BenchmarkDataLoader

app = cyclopts.App(
    name="download-benchmark",
    help="Download ProteoBench benchmark datasets for testing.",
)


@app.default
def download(
    output: Annotated[Path, cyclopts.Parameter(help="Output directory for downloaded data")] = Path("benchmark_data"),
    benchmark: Annotated[str, cyclopts.Parameter(help="Benchmark type")] = "dda_qexactive",
    software: Annotated[Optional[list[str]], cyclopts.Parameter(help="Filter by software name(s)")] = None,
    max_per_software: Annotated[Optional[int], cyclopts.Parameter(help="Max datasets per software")] = None,
    list_only: Annotated[bool, cyclopts.Parameter(help="Only list available software, don't download")] = False,
):
    """Download benchmark datasets from ProteoBench."""
    loader = BenchmarkDataLoader(benchmark, output_dir=str(output))

    if list_only:
        print(f"\nAvailable software in '{benchmark}':\n")
        for sw, versions in sorted(loader.list_software().items()):
            print(f"  {sw}: {', '.join(versions)}")
        return

    df = loader.get_datasets(software_filter=software, max_per_software=max_per_software)

    print(f"\nDownloaded {len(df)} datasets to {output}/")
    print("\nDatasets:")
    for _, row in df.iterrows():
        print(f"  {row['software_name']} {row['software_version']}: {row['local_dir']}")


@app.command
def list_benchmarks():
    """List available benchmark types."""
    print("\nAvailable benchmarks:\n")
    for key, info in BENCHMARK_REPOS.items():
        print(f"  {key}: {info['name']}")


def main():
    app()


if __name__ == "__main__":
    main()
