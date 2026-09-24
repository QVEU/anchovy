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

ONE DELIBERATE DIVERGENCE, and the golden CSV is deliberately NOT re-frozen to
absorb it -- see test_extract_diverges_from_the_original_only_where_intended.
The frozen file stays the original's real output, because it is the evidence
that says what the original did; a re-freeze would quietly delete that.
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

    # The four columns below diverge on purpose; they get their own test, which
    # pins the exact divergence rather than tolerating any difference.
    unchanged = [c for c in golden.columns if c not in DIVERGENT_COLUMNS]
    pd.testing.assert_frame_equal(
        result[unchanged], golden[unchanged],
        check_dtype=False,   # CSV round-trip can shift int/float dtypes harmlessly
        check_like=False,    # order already canonicalized
    )


# --------------------------------------------------------------------------- #
# The one place the port is deliberately not the original
# --------------------------------------------------------------------------- #
# barcodes.make_read_blocks used to reproduce blockDist's windowing exactly,
# including its off-by-one: the upper bound was len(seq) - query_len, so the
# LAST block -- the one flush with the end of the sequence -- was never
# generated. The signature lands exactly there whenever the barcode construct
# runs up to the cDNA that aligns, because the search window ends at the read's
# alignment offset. With the right block excluded, the best remaining candidate
# starts one base early, and the consequences are all silent:
#
#   readPos   one too low: the signature is reported at the wrong position
#   matchseq  the block shifted one base left
#   UMI       READ ONE BASE OFF, so reads are still grouped -- just into the
#             wrong molecules. Nothing downstream can tell.
#   minD      +2 on every read: a leading insertion and a trailing deletion, so
#             a PERFECT barcode is scored as two errors. That is not cosmetic
#             either: max_barcode_errors switches assignment to exact-or-within-k
#             substitutions, which drops indels by design, so setting it
#             discarded these reads outright.
#
# CBC and BC_ID do NOT diverge: the block-wise lookup recovered the right cell
# from the shifted block anyway, which is exactly why this went unnoticed.
DIVERGENT_COLUMNS = ["minD", "readPos", "matchseq", "UMI"]


def test_extract_diverges_from_the_original_only_where_intended(
        data_dir, golden_dir, tmp_run_dir):
    """The fix is a one-base shift, uniform across every read. Pin it."""
    extract = pytest.importorskip("anchovy.extract")
    golden = _load_golden(golden_dir)
    result = canonicalize(extract.run(
        sam=str(data_dir / "test.sam"),
        whitelist=str(data_dir / "whitelist.txt"),
        signature=SIGNATURE))

    assert (result.readPos - golden.readPos).eq(1).all(), (
        "the signature should now be found one base later -- at its true "
        "position, the last block the original never generated")

    # A perfect barcode scores len(signature) - 48 + 10 ... which for this
    # fixture's 26-N v2 signature is 10. Every fixture read carries an exact
    # whitelist barcode, so every one should score exactly that.
    assert result.minD.eq(10).all(), (
        f"fixture reads carry exact barcodes but scored {sorted(set(result.minD))}; "
        f"the original scored them all 12")
    assert golden.minD.eq(12).all(), (
        "the frozen golden was re-generated. Keep the original's output: it is "
        "the record of what the off-by-one did.")

    for col in ("matchseq", "UMI"):
        assert all(o[1:] == n[:-1] for o, n in zip(golden[col], result[col])), (
            f"{col} is no longer the golden's shifted one base right; the "
            f"divergence is not the windowing fix any more")
