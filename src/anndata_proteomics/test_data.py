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
