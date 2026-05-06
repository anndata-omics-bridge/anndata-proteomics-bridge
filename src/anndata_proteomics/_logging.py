"""Default loguru configuration for the package's CLIs and tools."""

from __future__ import annotations

import sys

from loguru import logger

DEFAULT_FORMAT = "<level>{level: <7}</level> | {message}"


def configure_default_sink(level: str = "INFO") -> None:
    """Reset loguru and install one stderr sink with a plain format.

    Called once at the top of each console-script's `main()`. Library modules
    do not call this — they just `from loguru import logger` and emit messages.
    Idempotent.
    """
    logger.remove()
    logger.add(sys.stderr, format=DEFAULT_FORMAT, level=level)
