"""
The block window the signature search runs over.

barcodes.make_read_blocks slices the search window into overlapping candidates
and best_query_match picks the closest. The bug these guard against is not a
wrong answer that looks wrong -- it is a candidate that was never offered, so
the search returned the best of what was left and every field derived from it
was one base off, silently.
"""

from __future__ import annotations

import pandas as pd

from anchovy.barcodes import best_query_match, make_read_blocks
from anchovy.config import ExtractConfig
from anchovy.extract import find_signature_positions

# v3: 22-base handle + 16-base barcode + 12-base UMI + 10-base handle.
SIGNATURE = "CTACACGACGCTCTTCCGATCT" + "N" * 28 + "TTTCTTATAT"


def test_the_final_flush_block_is_generated():
    """The last candidate is the one the signature actually occupies."""
    seq, query_len = "ACGT" * 25, 60          # 100 nt, 41 legal start positions
    blocks = make_read_blocks(seq, query_len)
    assert len(blocks) == len(seq) - query_len + 1
    assert blocks[-1] == seq[-query_len:], (
        "the flush block is missing; extract searches a window that ENDS at the "
        "read's alignment offset, which is exactly where the signature sits")


def test_a_signature_flush_with_the_window_end_is_found_exactly():
    """The case the off-by-one broke, and the reason it was invisible.

    Nothing failed: the search returned the block starting one base early, which
    still contains the barcode, so the right cell came out. What did not come
    out right was the UMI (read one base off, grouping the wrong molecules) and
    the distance (+2 on a perfect barcode -- a leading insertion and a trailing
    deletion -- which max_barcode_errors then drops as an indel).
    """
    window = "G" * 30 + SIGNATURE
    distance, position, block = best_query_match(window, SIGNATURE)
    assert (distance, position) == (0, 30)
    assert block == SIGNATURE


def test_a_signature_short_of_the_window_end_was_never_affected():
    """Which is why this survived real runs: it only bites when flush."""
    window = "G" * 30 + SIGNATURE + "ACGTACGTAC"
    assert best_query_match(window, SIGNATURE)[:2] == (0, 30)


def test_a_window_shorter_than_the_query_still_yields_one_block():
    """The max(1, ...) guard, kept from the original."""
    blocks = make_read_blocks("ACGTACGTAC", 60)
    assert len(blocks) == 1
    assert blocks[0] == "ACGTACGTAC"


# --------------------------------------------------------------------------- #
# Where the window ENDS, not how many blocks it yields
# --------------------------------------------------------------------------- #
# The tests above fix the window's left-to-right block range. This section is
# about its right edge, which extract pins at the read's alignment offset --
# i.e. it assumes minimap2 starts the alignment exactly where the construct
# ends. It does not always. The 3' handle is T-rich (TTTCTTATAT) and the cDNA
# that follows often starts with polyT, so the aligner can extend a base or
# two back into the handle. The window then truncates the construct, every
# block is shifted, and the barcode slice comes out of the wrong 16 bases.
#
# It does not look like a bug in the output. A one-base overrun reports as a
# clean +2 on every affected read -- indistinguishable from two sequencing
# errors in the barcode, which is how it survived a full run.
PREFIX, SUFFIX = SIGNATURE[:22], SIGNATURE[-10:]


def _read(barcode, umi, overrun, left=150, cdna=300):
    """A read whose alignment starts `overrun` bases inside the 3' handle."""
    construct = PREFIX + barcode + umi + SUFFIX
    seq = "A" * left + construct + "C" * cdna
    return seq, construct, left + len(construct) - overrun


def _matched_block(overrun, **cfg):
    barcode, umi = "ACGTACGTACGTACGT", "TTGGCCAATTGG"
    seq, construct, offset = _read(barcode, umi, overrun)
    frame = pd.DataFrame({"seq": [seq], "offset": [offset], "template": ["t"]})
    out = find_signature_positions(
        frame, SIGNATURE, ExtractConfig(signature=SIGNATURE, nthreads=1, **cfg))
    return out.matchseq.iloc[0], construct, barcode


def test_an_alignment_overrunning_the_handle_still_finds_the_construct():
    """One base of overrun used to shift the whole block, and the barcode with it."""
    for overrun in (0, 1, 2, 3):
        block, construct, barcode = _matched_block(overrun)
        assert block == construct, (
            f"alignment overrunning the 3' handle by {overrun} base(s) shifted "
            f"the matched block; the barcode slice becomes "
            f"{block[22:38]!r} instead of {barcode!r}")
        assert block[22:38] == barcode


def test_the_window_reaches_past_the_offset_by_the_configured_pad():
    """The pad is what makes the construct interior rather than flush."""
    assert ExtractConfig().downstream_window > 0

    # With no pad the window stops dead at the offset and the overrun bites.
    block, construct, _ = _matched_block(2, downstream_window=0)
    assert block != construct, (
        "with downstream_window=0 this is the pre-fix behaviour and the block "
        "should still be shifted -- if it is not, the fixture stopped "
        "reproducing the bug")


def test_padding_cannot_disturb_a_read_whose_alignment_starts_cleanly():
    """The common case must be untouched by the pad."""
    block, construct, barcode = _matched_block(0)
    assert block == construct and block[22:38] == barcode


# --------------------------------------------------------------------------- #
# Telling a misplaced window apart from a misread barcode
# --------------------------------------------------------------------------- #
def _run_pass_two(matchseqs, capsys):
    """Drive pass 2 over blocks we control, and return what it printed."""
    from anchovy.extract import assign_barcodes

    barcode = "ACGTACGTACGTACGT"
    whitelist = pd.DataFrame({"CBC": [barcode, "TTTTTTTTTTTTTTTT"]})
    frame = pd.DataFrame({
        "read": [f"r{i}" for i in range(len(matchseqs))],
        "seq": ["A" * 80] * len(matchseqs),
        "minPos": [0] * len(matchseqs),
        "matchseq": matchseqs,
        "readLen": [80] * len(matchseqs),
        "offset": [0] * len(matchseqs),
        "minD": [0] * len(matchseqs),
    })
    assign_barcodes(frame, SIGNATURE, whitelist,
                    ExtractConfig(signature=SIGNATURE, nthreads=1))
    return capsys.readouterr().out


def test_a_shifted_block_is_reported_as_a_window_problem(capsys):
    """The line that would have caught this from the run's own output."""
    barcode, umi = "ACGTACGTACGTACGT", "TTGGCCAATTGG"
    construct = PREFIX + barcode + umi + SUFFIX
    shifted = ("GG" + construct)[:len(construct)]    # window truncated by two

    out = _run_pass_two([shifted] * 4, capsys)
    assert "do not begin with the signature's 5' handle" in out
    assert "100.0%" in out
    assert "misplaced search window, not misread barcodes" in out


def test_clean_blocks_say_nothing_about_the_window(capsys):
    """The warning must not cry wolf on an ordinary run."""
    barcode, umi = "ACGTACGTACGTACGT", "TTGGCCAATTGG"
    construct = PREFIX + barcode + umi + SUFFIX

    out = _run_pass_two([construct] * 4, capsys)
    assert "5' handle" not in out, (
        "a clean run must not be told its window is misplaced")
