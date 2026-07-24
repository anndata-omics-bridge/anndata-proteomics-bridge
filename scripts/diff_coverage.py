"""Run changed-line coverage against a configurable Git comparison branch."""

from __future__ import annotations

import os
import subprocess


def main() -> None:
    """Run diff-cover with the local or CI comparison branch."""
    compare_branch = os.environ.get("DIFF_COVER_BASE", "main")
    subprocess.run(
        [
            "diff-cover",
            "coverage.xml",
            f"--compare-branch={compare_branch}",
            "--fail-under=100",
        ],
        check=True,
    )


if __name__ == "__main__":
    main()
