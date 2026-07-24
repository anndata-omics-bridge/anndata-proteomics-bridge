"""Build and inspect the APB distribution without installing from PyPI."""

from __future__ import annotations

import configparser
import subprocess
import tempfile
import zipfile
from pathlib import Path

_ENTRY_POINTS = {
    "apb": "anndata_proteomics.scripts.cli:main",
    "apb-testdata": "anndata_proteomics.scripts.extract_raw_file_db:main",
}
_REQUIRED_FILES = {
    "anndata_proteomics/__init__.py",
    "anndata_proteomics/py.typed",
    "anndata_proteomics/modifications/unimod_registry.toml",
    "anndata_proteomics/parsing_rules/_schema/parse_rule.schema.json",
    "anndata_proteomics/parsing_rules/maxquant/rules.json",
    "anndata_proteomics/proteobench/APACHE-2.0.txt",
    "anndata_proteomics/proteobench/mapper.csv",
}


def _verify_wheel(wheel: Path) -> None:
    with zipfile.ZipFile(wheel) as archive:
        names = set(archive.namelist())
        if missing := _REQUIRED_FILES - names:
            raise RuntimeError(f"Wheel is missing package files: {sorted(missing)}")

        entry_point_files = [name for name in names if name.endswith(".dist-info/entry_points.txt")]
        if len(entry_point_files) != 1:
            raise RuntimeError(f"Expected one entry_points.txt, found {entry_point_files}")

        parser = configparser.ConfigParser()
        parser.read_string(archive.read(entry_point_files[0]).decode())
        console_scripts = parser["console_scripts"]
        for name, expected in _ENTRY_POINTS.items():
            actual = console_scripts.get(name)
            if actual != expected:
                raise RuntimeError(f"Unexpected {name} entry point: {actual!r}")


def main() -> None:
    """Build an sdist and wheel, then inspect the wheel's public contract."""
    with tempfile.TemporaryDirectory(prefix="apb-package-") as temp_dir:
        output_dir = Path(temp_dir)
        subprocess.run(["uv", "build", "--out-dir", str(output_dir)], check=True)
        wheels = list(output_dir.glob("*.whl"))
        if len(wheels) != 1:
            raise RuntimeError(f"Expected one wheel, found {wheels}")
        _verify_wheel(wheels[0])
        print(f"Package smoke passed: {wheels[0].name}")


if __name__ == "__main__":
    main()
