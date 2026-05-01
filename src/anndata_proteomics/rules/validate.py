"""Validate parsing-rule TOMLs and produce CLI-friendly results."""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

from anndata_proteomics.rules.loader import load_rule
from anndata_proteomics.rules.registry import iter_packaged_rules, packaged_rules_root
from anndata_proteomics.rules.schema import ParseRule


@dataclass(frozen=True)
class ValidationResult:
    path: Path
    ok: bool
    error: str | None = None
    rule: ParseRule | None = None


def validate_file(path: Path | str) -> ValidationResult:
    """Validate one TOML file. Never raises; failures come back in the result."""
    p = Path(path)
    try:
        rule = load_rule(p)
    except Exception as e:
        notes = getattr(e, "__notes__", [])
        msg = "; ".join([str(e)] + list(notes)) if notes else str(e)
        return ValidationResult(path=p, ok=False, error=msg)
    return ValidationResult(path=p, ok=True, rule=rule)


def validate_all_packaged() -> list[ValidationResult]:
    """Validate every packaged rule. Sorted by path."""
    return [validate_file(p) for p in iter_packaged_rules()]


def main(argv: list[str] | None = None) -> int:
    """Console entrypoint: PASS/FAIL per rule, summary line, exit 1 if any failed."""
    pkg_parent = packaged_rules_root().parent
    results = validate_all_packaged()
    for r in results:
        rel = r.path.relative_to(pkg_parent) if r.path.is_relative_to(pkg_parent) else r.path
        if r.ok:
            print(f"PASS  {rel}")
        else:
            err_summary = (r.error or "(no error message)").splitlines()[0]
            print(f"FAIL  {rel}: {err_summary}")
    failed = sum(1 for r in results if not r.ok)
    print(f"{len(results)} packaged rule(s) checked, {failed} failed.")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
