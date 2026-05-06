"""Bridge loguru output into pytest's capsys/capfd.

Loguru's default sink captures `sys.stderr` at handler-registration time, which
runs once at module import. Pytest's `capsys` monkeypatches `sys.stderr` per
test, so without a bridge, loguru output bypasses capture and tests can't
assert on it. We replace the default sink with one whose writer callable
looks up `sys.stderr` at *write* time, picking up whatever pytest has patched
in for the current test.
"""

from __future__ import annotations

import sys

import pytest
from loguru import logger


@pytest.fixture(autouse=True)
def _loguru_to_pytest_capsys():
    logger.remove()
    logger.add(
        lambda msg: sys.stderr.write(msg),
        format="{level: <7} | {message}",
        level="DEBUG",
    )
    yield
    # Drop every sink — including any added by code under test (e.g.
    # configure_default_sink() or per-run file sinks). Avoids leaking sinks
    # across tests and tolerates main()-style code that calls logger.remove().
    logger.remove()
