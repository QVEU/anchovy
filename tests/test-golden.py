"""
test_golden.py -- characterization test guarding the extract stage.

This asserts that the migrated package produces the SAME extract output as the
original anchovy.py did on the fixture input. It's the safety net for the whole
migration: run it after every change to extract.py / barcodes.py / io.py, and if
it stays green, behavior was preserved.

Right now `run_extract` is written against the FUTURE package API
(anchovy.extract.run) that we haven't built yet, so this test will fail to import
or be skipped until extract.py exists. That's intentional and correct for
tests-first: the test is written before the code it tests, describing the
behavior the code must satisfy. We watch it go from red to green as we migrate.
"""

from __future__ import annotations

import pandas as pd
import pytest

# The canonicalizer lives with the golden generator; reuse it so the test and the
# frozen output are canonicalized identically. If you prefer, move canonicalize()
# into anchovy.io and import it from there instead.
from make_golden import canonicalize, SIGNATURE


def _load_golden(golden_dir):
    path = golden_dir / "extract_anchovy.csv"
    if not path.exists():
        pytest.skip(
            "Golden output missing. Generate it first:\n"
            "  python tests/make_fixtures.py\n"
            "  python tests/make_golden.py /path/to/original/anchovy.py"
        )
    return canonicalize(pd.read_csv(path))


def test_extract_matches_golden(data_dir, golden_dir, tmp_run_dir):
    """Migrated extract stage reproduces the original's output exactly."""
    # Imported lazily so the whole test module doesn't error out before the
    # package function exists -- lets the rest of the suite run in the meantime.
    extract = pytest.importorskip("anchovy.extract")

    golden = _load_golden(golden_dir)

    sam = str(data_dir / "test.sam")
    whitelist = str(data_dir / "whitelist.txt")

    # FUTURE API: extract.run returns the anchovy DataFrame (or writes a CSV we
    # read back). We'll define this signature when we build extract.py; the test
    # documents the contract we want.
    result = extract.run(sam=sam, whitelist=whitelist, signature=SIGNATURE)
    result = canonicalize(result)

    # Compare column sets first for a clearer failure message than a full-frame
    # diff when the schema itself drifted.
    assert list(result.columns) == list(golden.columns), (
        f"Column mismatch.\n got: {list(result.columns)}\n exp: {list(golden.columns)}"
    )

    pd.testing.assert_frame_equal(
        result, golden,
        check_dtype=False,   # CSV round-trip can shift int/float dtypes harmlessly
        check_like=False,    # order already canonicalized
    )
