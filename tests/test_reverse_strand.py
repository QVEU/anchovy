"""
test_reverse_strand.py -- reads stored on the reverse strand are extracted
exactly as they would be on the forward strand.

A SAM stores a read that maps to the reverse strand (FLAG 0x10) reverse-
complemented, with its CIGAR in reference order. The 10X construct then sits at
the RIGHT end of SEQ, beside the trailing soft clip, while the signature search
looks only LEFT of the alignment start. Before the fix every such read was
searched in the wrong place. These tests take the forward-strand fixture, store
each mapped read the way an aligner would store its reverse-strand twin, and
check that extract gives the same answer for both.
"""

from __future__ import annotations

import re

import pandas as pd
import pytest

from anchovy import extract
from anchovy.io import parse_cigar_lengths, read_sam, reverse_complement
from anchovy.schema import AnchovyColumns as C, SamColumns as S

from make_golden import SIGNATURE


class FakeRead:
    def __init__(self, cigar, query_name="test_read"):
        self.cigar = cigar
        self.query_name = query_name


def _reverse_cigar(cigar: str) -> str:
    return "".join(reversed(re.findall(r"\d+[MIDNSHP=X]", cigar)))


def _as_reverse_strand(src, dst) -> int:
    """Rewrite every mapped forward read as its reverse-strand SAM record."""
    flipped = 0
    out = []
    for line in src.read_text().splitlines():
        fields = line.split("\t")
        if line.startswith("@") or fields[2] == "*" or fields[5] == "*":
            out.append(line)
            continue
        fields[1] = str(int(fields[1]) | 16)
        fields[5] = _reverse_cigar(fields[5])
        fields[9] = reverse_complement(fields[9])
        fields[10] = fields[10][::-1]
        out.append("\t".join(fields))
        flipped += 1
    dst.write_text("\n".join(out) + "\n")
    return flipped


def test_reverse_complement():
    assert reverse_complement("AACGTN") == "NACGTT"
    assert reverse_complement("acgT") == "Acgt"


def test_reversed_cigar_puts_the_trailing_clip_in_front():
    # Stored as 100M70S: in sequencing order that is 70S100M.
    assert parse_cigar_lengths(FakeRead([(0, 100), (4, 70)]), reverse=True) \
        == [170, 100, 70]
    # The forward reading is unchanged by the new argument's default.
    assert parse_cigar_lengths(FakeRead([(0, 100), (4, 70)])) == [100, 100, 0]


def test_reverse_strand_reads_are_read_back_in_sequencing_order(
        data_dir, tmp_path):
    rev_sam = tmp_path / "rev.sam"
    assert _as_reverse_strand(data_dir / "test.sam", rev_sam) > 0

    fwd = read_sam(str(data_dir / "test.sam"), len(SIGNATURE))
    rev = read_sam(str(rev_sam), len(SIGNATURE))

    cols = [S.READ, S.SEQ, S.READ_LEN, S.CLIP_READ_LEN, S.OFFSET, S.LENGTH]
    pd.testing.assert_frame_equal(fwd[cols], rev[cols])
    # FLAG still says which reads were turned.
    assert set(rev[S.FLAG]) == {"16"}


def test_extract_assigns_reverse_strand_reads_like_forward_ones(
        data_dir, tmp_path, capsys):
    rev_sam = tmp_path / "rev.sam"
    _as_reverse_strand(data_dir / "test.sam", rev_sam)
    # Reads below the length filter are dropped before they are turned.
    n = len(read_sam(str(rev_sam), len(SIGNATURE)))
    capsys.readouterr()
    whitelist = str(data_dir / "whitelist.txt")

    fwd = extract.run(sam=str(data_dir / "test.sam"), whitelist=whitelist,
                      signature=SIGNATURE)
    capsys.readouterr()
    rev = extract.run(sam=str(rev_sam), whitelist=whitelist,
                      signature=SIGNATURE)
    out = capsys.readouterr().out

    assert len(fwd) > 0
    assert "Reverse-strand reads re-oriented: {}".format(n) in out
    pd.testing.assert_frame_equal(
        fwd.sort_values(C.READ).reset_index(drop=True),
        rev.sort_values(C.READ).reset_index(drop=True))


def test_a_mixed_strand_file_loses_no_reads(data_dir, tmp_path):
    """Half the reads flipped, as in an unoriented long-read run."""
    lines = (data_dir / "test.sam").read_text().splitlines()
    rev_sam = tmp_path / "rev.sam"
    _as_reverse_strand(data_dir / "test.sam", rev_sam)
    rev_lines = rev_sam.read_text().splitlines()
    mixed = [r if i % 2 else f for i, (f, r) in
             enumerate(zip(lines, rev_lines))]
    mixed_sam = tmp_path / "mixed.sam"
    mixed_sam.write_text("\n".join(mixed) + "\n")

    whitelist = str(data_dir / "whitelist.txt")
    fwd = extract.run(sam=str(data_dir / "test.sam"), whitelist=whitelist,
                      signature=SIGNATURE)
    got = extract.run(sam=str(mixed_sam), whitelist=whitelist,
                      signature=SIGNATURE)
    assert len(got) == len(fwd)
    pd.testing.assert_frame_equal(
        fwd.sort_values(C.READ).reset_index(drop=True),
        got.sort_values(C.READ).reset_index(drop=True))
