"""Generate source-document and effective-rule JSON Schemas from Pydantic."""

from __future__ import annotations

import json
from pathlib import Path

from loguru import logger

from anndata_proteomics.rules.schema import ParseRule, ParseRuleDocument


def main() -> None:
    output_directory = Path(__file__).resolve().parent.parent / "parsing_rules" / "_schema"
    output_directory.mkdir(parents=True, exist_ok=True)
    schemas = {
        "parse_rule.schema.json": ParseRule.model_json_schema(),
        "parse_rule_document.schema.json": ParseRuleDocument.model_json_schema(),
    }
    for filename, schema in schemas.items():
        output = output_directory / filename
        output.write_text(json.dumps(schema, indent=2) + "\n")
        logger.info(f"wrote {output}")


if __name__ == "__main__":
    main()
