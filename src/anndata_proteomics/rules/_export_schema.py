"""Generate parse_rule.schema.json from the pydantic ParseRule model."""

from __future__ import annotations

import json
from pathlib import Path

from loguru import logger

from anndata_proteomics.rules.schema import ParseRule


def main() -> None:
    schema = ParseRule.model_json_schema()
    out = (
        Path(__file__).resolve().parent.parent
        / "parsing_rules"
        / "_schema"
        / "parse_rule.schema.json"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(schema, indent=2) + "\n")
    logger.info(f"wrote {out}")


if __name__ == "__main__":
    main()
