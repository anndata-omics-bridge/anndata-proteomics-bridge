"""Convert one cached ProteoBench dataset to AnnData/MuData — run as a subprocess.

The GUI (ui_test_tool.py) launches this with ``python -m anndata_proteomics.scripts.convert_one``
as a background job (see scripts/jobrunner.py), streaming this module's stdout+stderr to a log
file so the UI can tail it live. Heavy targets (fragment, mudata) can take seconds and GBs, which
is exactly why they run out-of-process rather than blocking the marimo event loop.

    python -m anndata_proteomics.scripts.convert_one \\
        --input <input_file_path> --slug <slug> --target <ion|...|mudata> \\
        --params <param_file> --outdir <dir>

The param file is mandatory: it provides the software version that selects the rule variant.
"""

from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path

from anndata_proteomics.scripts import _ui_support as ui


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="path under test_data_download/json_dir/")
    parser.add_argument("--slug", required=True, help="vendor slug, e.g. diann")
    parser.add_argument("--target", required=True, help="a level or 'mudata'")
    parser.add_argument("--params", required=True, help="param file (gives the software version)")
    parser.add_argument("--outdir", required=True, help="directory for result + log")
    args = parser.parse_args(argv)

    outdir = Path(args.outdir).expanduser()
    outdir.mkdir(parents=True, exist_ok=True)

    def log(message: str) -> None:
        print(message, flush=True)

    log(f"convert slug={args.slug} target={args.target}")
    try:
        obj = ui.convert_target(
            args.input, args.slug, args.target, param_path=args.params, log=log
        )
        if args.target == ui.MUDATA:
            out = outdir / "result.h5mu"
            obj.write(out)
        else:
            out = outdir / "result.h5ad"
            obj.write_h5ad(out)
        log(f"wrote {out}")
        log("summary:")
        log(json.dumps(ui.summarize(obj), indent=2, default=str))
        log("DONE")
        return 0
    except Exception:  # full traceback → stderr → merged into the log file
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
