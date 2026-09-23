"""
test_io.py -- unit tests for the pure functions in anchovy.io.

These test only functions with no file I/O, no multiprocessing, and no global
state: inputs in, outputs out. That's what makes them fast (milliseconds),
deterministic (no flakiness), and precise (a failure names the exact function).

The heavier machinery -- reading real SAM/BAM, the multiprocessing pools, the
full extract pipeline -- is covered by the golden/integration test instead,
because it can't be exercised as a pure function.

WHY A FAKE READ CLASS
---------------------
parse_cigar_lengths takes a pysam read object, but it only ever touches read.cigar
(a list of (op, length) tuples). So we don't need real pysam or a real file -- a
tiny stand-in with a .cigar attribute is enough. This is a general testing move:
depend on the smallest surface you can fake, not the whole heavy object.
"""

from __future__ import annotations

import pytest

from anchovy.io import parse_cigar_lengths, parse_description


class FakeRead:
    """Minimal stand-in for a pysam read: only exposes .cigar."""
    def __init__(self, cigar):
        self.cigar = cigar


# --------------------------------------------------------------------------- #
# parse_cigar_lengths
# CIGAR ops: 0=match(M) 1=insertion(I) 2=deletion(D) 3=skip(N)
#            4=soft-clip(S) 5=hard-clip(H) 6=padding(P)
# Contract (preserved from the original loadSAM loop):
#   M/I/N -> add to BOTH readLen and clipReadLen
#   D     -> flips switch only (no length added)
#   S     -> add to readLen ONLY if leading (switch still 0)
#   H/P   -> ignored
#   offset = readLen - clipReadLen
# --------------------------------------------------------------------------- #
def test_cigar_all_match():
    # 100M: pure match. Both lengths 100, offset 0.
    assert parse_cigar_lengths(FakeRead([(0, 100)])) == [100, 100, 0]


def test_cigar_leading_soft_clip_counts_toward_readlen_only():
    # 10S90M: leading soft clip adds to readLen but not clipReadLen.
    # readLen = 10 + 90 = 100; clipReadLen = 90; offset = 10.
    assert parse_cigar_lengths(FakeRead([(4, 10), (0, 90)])) == [100, 90, 10]


def test_cigar_trailing_soft_clip_ignored():
    # 90M10S: the soft clip is NOT leading (switch already 1 after the M),
    # so it's ignored entirely. readLen = clipReadLen = 90; offset = 0.
    assert parse_cigar_lengths(FakeRead([(0, 90), (4, 10)])) == [90, 90, 0]


def test_cigar_insertion_and_skip_count_both():
    # 50M5I45N: I and N both add to both lengths like M does.
    # both = 50 + 5 + 45 = 100; offset = 0.
    assert parse_cigar_lengths(FakeRead([(0, 50), (1, 5), (3, 45)])) == [100, 100, 0]


def test_cigar_deletion_adds_nothing():
    # 50M5D50M: deletion contributes no length, only flips switch.
    # both = 50 + 50 = 100; offset = 0.
    assert parse_cigar_lengths(FakeRead([(0, 50), (2, 5), (0, 50)])) == [100, 100, 0]


def test_cigar_hard_clip_and_padding_ignored():
    # 5H100M5P: hard clip and padding are ignored. both = 100; offset = 0.
    assert parse_cigar_lengths(FakeRead([(5, 5), (0, 100), (6, 5)])) == [100, 100, 0]


def test_cigar_empty():
    # No ops -> all zeros. Documents the edge case explicitly.
    assert parse_cigar_lengths(FakeRead([])) == [0, 0, 0]


# --------------------------------------------------------------------------- #
# parse_description  (the shared Python/R contract for FASTA metadata)
# Format: whitespace-separated tokens including 'coverage:<n>' and 'length:<n>'
# --------------------------------------------------------------------------- #
def test_description_basic():
    result = parse_description("AAACCC ref coverage:37 length:9800")
    assert result == {"coverage": 37.0, "length": 9800.0}


def test_description_order_independent():
    # Parsing is by prefix, not position -- order and extra tokens don't matter.
    # This is the whole point of the rewrite vs the original's positional split.
    result = parse_description("length:120 junk coverage:5 more_junk")
    assert result == {"coverage": 5.0, "length": 120.0}


def test_description_missing_fields_return_none():
    assert parse_description("no metadata here") == {"coverage": None, "length": None}


def test_description_partial():
    assert parse_description("x coverage:12") == {"coverage": 12.0, "length": None}


def test_description_non_numeric_value_is_none():
    # A malformed value is tolerated (returns None) rather than raising,
    # matching the tolerant spirit of the original parsing.
    assert parse_description("coverage:abc length:50") == {
        "coverage": None, "length": 50.0,
    }


def test_description_float_values():
    assert parse_description("coverage:12.5 length:9.0") == {
        "coverage": 12.5, "length": 9.0,
    }


# --------------------------------------------------------------------------- #
# Whitelist validation
# --------------------------------------------------------------------------- #
# A whitelist anchovy cannot search with does not fail on its own.
# build_barcode_query_blocks builds one template per barcode as
# signature[0:22] + barcode + N*umi + signature[-10:], so a wrong-width barcode
# shifts everything after it -- and without max_barcode_errors each read is
# still assigned to its NEAREST entry rather than rejected. The run finishes,
# every cell is wrong, and nothing says so.
import pytest

from anchovy.io import read_whitelist, validate_whitelist
from anchovy.schema import SIGNATURE_BARCODE_LEN


def test_the_barcode_width_is_the_same_in_both_chemistries():
    """The check is chemistry-independent, and that is the point.

    v2 and v3 differ only in UMI width (10 vs 12 nt); the cell barcode is 16 in
    both. So a length that is not 16 is wrong whatever `chemistry` says, and a
    correct length tells you nothing about which chemistry you have.
    """
    assert SIGNATURE_BARCODE_LEN == 16


def test_a_good_whitelist_is_accepted():
    validate_whitelist(["AAACCTGAGAAACCAT", "AAACCTGAGAAACCGC"])


def test_the_cellranger_gem_suffix_is_named_with_its_fix(tmp_path):
    """The likeliest thing to be wrong, for anyone exporting from Seurat."""
    p = tmp_path / "barcodes.tsv"
    p.write_text("AAACCTGAGAAACCAT-1\nAAACCTGAGAAACCGC-1\n")
    with pytest.raises(ValueError, match="GEM suffix") as exc:
        read_whitelist(str(p))
    # The message has to carry the remedy, not just the diagnosis.
    assert "cut -d- -f1" in str(exc.value)


def test_a_header_line_is_caught(tmp_path):
    p = tmp_path / "barcodes.tsv"
    p.write_text("barcode\nAAACCTGAGAAACCAT\n")
    with pytest.raises(ValueError, match="not all the same length"):
        read_whitelist(str(p))


def test_a_right_width_non_barcode_is_caught():
    # A header would have to be a plausible 16-mer to reach this, but a quoted
    # CSV export manages it.
    # 14 bases in quotes is exactly 16 characters, so it passes the width
    # check and only the alphabet catches it.
    quoted = '"AAACCTGAGAAACC"'
    assert len(quoted) == SIGNATURE_BARCODE_LEN
    with pytest.raises(ValueError, match="other than A/C/G/T/N"):
        validate_whitelist([quoted, "AAACCTGAGAAACCAT"])


def test_a_uniform_but_wrong_width_says_it_is_not_a_chemistry_thing():
    with pytest.raises(ValueError, match="NOT a v2/v3 difference"):
        validate_whitelist(["AAACCTGAGAAACCA", "AAACCTGAGAAACCG"])


def test_an_empty_whitelist_is_caught(tmp_path):
    p = tmp_path / "barcodes.tsv"
    p.write_text("\n\n")
    with pytest.raises(ValueError, match="no barcodes"):
        read_whitelist(str(p))


def test_blank_lines_do_not_become_empty_barcodes(tmp_path):
    """An empty barcode is an entry every poor read can match."""
    p = tmp_path / "barcodes.tsv"
    p.write_text("AAACCTGAGAAACCAT\n\nAAACCTGAGAAACCGC\n")
    assert list(read_whitelist(str(p)).CBC) == ["AAACCTGAGAAACCAT",
                                                "AAACCTGAGAAACCGC"]
