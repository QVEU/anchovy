"""
conftest.py -- shared pytest fixtures.

pytest automatically discovers this file and makes the fixtures below available
to every test without importing them. Fixtures are pytest's dependency-injection
mechanism: a test function that takes an argument named `data_dir` receives
whatever the `data_dir` fixture returns.

We centralize path resolution here so no test hardcodes a path, and so moving the
test data means editing one place.
"""

from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture
def data_dir() -> Path:
    """Directory holding test inputs and golden outputs.

    Resolved relative to THIS file, not the current working directory, so tests
    pass no matter where pytest is invoked from.
    """
    return Path(__file__).parent / "data"


@pytest.fixture
def golden_dir(data_dir: Path) -> Path:
    """Directory holding the frozen 'golden' reference outputs."""
    return data_dir / "golden"


@pytest.fixture
def tmp_run_dir(tmp_path: Path) -> Path:
    """A fresh temporary directory for a test to write outputs into.

    `tmp_path` is a built-in pytest fixture giving a unique temp dir per test.
    We wrap it so tests express intent ('a place to run the pipeline') rather
    than plumbing. Nothing here persists between tests -- that isolation is the
    point: a test can never accidentally read another test's leftovers.
    """
    return tmp_path
