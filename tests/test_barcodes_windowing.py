"""
The block window the signature search runs over.

barcodes.make_read_blocks slices the search window into overlapping candidates
and best_query_match picks the closest. The bug these guard against is not a
wrong answer that looks wrong -- it is a candidate that was never offered, so
the search returned the best of what was left and every field derived from it
was one base off, silently.
"""

from __future__ import annotations

from anchovy.barcodes import best_query_match, make_read_blocks

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
