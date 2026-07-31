"""Validate software-version parsing-rule documents and produce CLI-friendly results."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from loguru import logger

from anndata_proteomics.rules.loader import load_rule_document
from anndata_proteomics.rules.registry import iter_packaged_documents, packaged_rules_root
from anndata_proteomics.rules.schema import ParseRuleDocument


@dataclass(frozen=True, slots=True)
class ValidationResult:
    """Operational validation outcome for one source document."""

    path: Path
    ok: bool
    error: str | None = None
    document: ParseRuleDocument | None = None


def validate_file(path: Path | str) -> ValidationResult:
    """Validate one JSON document and every merged level without raising."""
    source = Path(path)
    try:
        document = load_rule_document(source)
    except (OSError, UnicodeDecodeError, ValueError) as exc:
        notes = getattr(exc, "__notes__", [])
        message = "; ".join([str(exc), *notes]) if notes else str(exc)
        return ValidationResult(path=source, ok=False, error=message)
    return ValidationResult(path=source, ok=True, document=document)


def validate_all_packaged() -> list[ValidationResult]:
    """Validate every packaged software-version document in path order."""
    return [validate_file(path) for path in iter_packaged_documents()]


def log_and_exit_code(results: list[ValidationResult]) -> int:
    """Log PASS/FAIL per document and return zero only when all are valid."""
    package_parent = packaged_rules_root().parent
    for result in results:
        path = result.path
        relative = path.relative_to(package_parent) if path.is_relative_to(package_parent) else path
        if result.ok:
            levels = ", ".join(result.document.levels) if result.document is not None else ""
            logger.info(f"PASS  {relative} [{levels}]")
        else:
            summary = (result.error or "(no error message)").splitlines()[0]
            logger.error(f"FAIL  {relative}: {summary}")
    failed = sum(not result.ok for result in results)
    logger.info(f"{len(results)} document(s) checked, {failed} failed.")
    return 0 if failed == 0 else 1
