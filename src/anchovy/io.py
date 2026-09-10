"""
io.py -- shared readers and writers for every pipeline stage.

WHY THIS FILE EXISTS
--------------------
The original coupling wasn't just positional indices; it was that each stage
re-implemented how to read the previous stage's output. anchovy.py parsed SAM by
hand AND via pysam; ConsensusTool.py parsed FASTA description strings by
whitespace offset; the R script parsed those same descriptions AGAIN, differently.
Every one of those parsers was a private copy of an implicit format.

This module makes each read/write operation exist exactly once, built on the
column names in schema.py. When a stage needs the SAM table, it calls read_sam().
When it needs coverage/length out of a description, it calls parse_description() --
the SAME function the R stage should call (via a documented format), so the two
languages can't drift apart.

BEHAVIOR NOTE
-------------
These functions reproduce the original parsing behavior exactly, including its
quirks (e.g. the manual text parse skips only lines beginning with '@'). Where
the original had a latent bug, this module preserves the behavior and flags it in
a comment rather than silently "fixing" it -- behavior changes belong in their own
clearly-labeled commits, validated against the golden tests.
"""

from __future__ import annotations

import sys

import pandas as pd
import pysam

from anchovy.schema import (
    SamColumns,
    AnchovyColumns,
    DESCRIPTION_COVERAGE_KEY,
    DESCRIPTION_LENGTH_KEY,
)


# --------------------------------------------------------------------------- #
# SAM reading
# --------------------------------------------------------------------------- #
def parse_cigar_lengths(read) -> list[int]:
    """Compute [readLen, clipReadLen, offset] from one pysam read's CIGAR.

    Extracted from the original loadSAM inline loop, unchanged in semantics:
      - op 0 (match), 1 (insertion), 3 (skip): add to both readLen and CreadLen
      - op 2 (deletion): advances the 'switch' flag only
      - op 4 (soft clip): add to readLen only if leading (switch == 0)
      - op 5 (hard clip), 6 (padding): ignored
    'offset' is returned as readLen - clipReadLen.

    The bare try/except that swallowed errors as "Problem" is preserved to match
    original behavior; a stricter version belongs in a separate error-handling
    commit (Phase 6), not here.
    """
    read_len = 0
    clip_read_len = 0
    switch = 0
    for cigar_type, cigar_length in read.cigar:
        try:
            if cigar_type == 0:            # match
                read_len += cigar_length
                clip_read_len += cigar_length
                switch = 1
            elif cigar_type == 1:          # insertion
                read_len += cigar_length
                clip_read_len += cigar_length
                switch = 1
            elif cigar_type == 2:          # deletion
                switch = 1
            elif cigar_type == 3:          # skip
                read_len += cigar_length
                clip_read_len += cigar_length
                switch = 1
            elif cigar_type == 4:          # soft clipping
                if switch == 0:            # only count leading soft clipping
                    read_len += cigar_length
                switch = 1
            elif cigar_type == 5:          # hard clipping
                pass
            elif cigar_type == 6:          # padding
                pass
            else:
                print("Wrong CIGAR number")
                sys.exit(1)
        except Exception:
            print("Problem")
    return [read_len, clip_read_len, read_len - clip_read_len]


def read_sam(path: str, min_read_length: int) -> pd.DataFrame:
    """Read a mapped SAM into a DataFrame with CIGAR-derived length columns.

    Reproduces the original loadSAM exactly:
      1. Manual text parse of the first 11 tab-separated fields, skipping only
         lines starting with '@' (pandas struggled with SAM headers, per the
         original's note).
      2. pysam pass to compute CIGAR lengths for mapped reads.
      3. Drop unmapped (template == '*'), attach length columns, filter by length.

    The two passes over the file are the original's design. Collapsing them into
    one is a performance change for a later, separately-validated commit, not a
    behavior-preserving migration.

    Args:
        path: path to the input SAM file.
        min_read_length: keep reads with length strictly greater than this.
            (Callers pass config.extract.effective_min_read_length(), which
            defaults to the signature length -- the original's `minL = quL`.)
    """
    # 1. Manual text parse of the 11 core SAM fields.
    with open(path, "r") as handle:
        rows = [line.split("\t")[0:11] for line in handle if not line.startswith("@")]

    # 2. pysam pass for CIGAR accounting on mapped reads only.
    print("Parsing Cigars...")
    sam_fp = pysam.Samfile(path, "rb")
    cigars = [parse_cigar_lengths(read) for read in sam_fp if not read.is_unmapped]
    print("Done.")

    print("Total Candidate Reads: {}".format(len(rows)))

    # 3. Assemble the frame, attach derived columns, filter.
    df = pd.DataFrame(rows, columns=SamColumns.ORDER)
    df = df[df[SamColumns.TEMPLATE] != "*"]
    df[[SamColumns.READ_LEN, SamColumns.CLIP_READ_LEN, SamColumns.OFFSET]] = cigars
    df[SamColumns.LENGTH] = df[SamColumns.SEQ].apply(len)
    df = df[df[SamColumns.LENGTH] > min_read_length]
    return df


def read_whitelist(path: str) -> pd.DataFrame:
    """Read a 10X barcode whitelist into a single-column DataFrame.

    Reproduces loadBC: take the first tab-separated field of each line, stripped.
    """
    with open(path, "r") as handle:
        barcodes = [line.split("\t")[0].strip() for line in handle]
    print("Total Cell Barcodes: {}".format(len(barcodes)))
    return pd.DataFrame(barcodes, columns=["CBC"])


# --------------------------------------------------------------------------- #
# anchovy CSV writing / reading (the extract -> fasta handoff)
# --------------------------------------------------------------------------- #
def write_anchovy_csv(df: pd.DataFrame, path: str, chunksize: int = 50000) -> None:
    """Write the anchovy extract table, enforcing the canonical column order.

    Writing through schema.AnchovyColumns.ORDER is what makes the downstream
    fasta stage safe: even code that still used positional access would stay
    correct because the order is now guaranteed here. chunksize mirrors the
    original's OOM-avoidance setting.
    """
    df = df[AnchovyColumns.ORDER]
    df.to_csv(path, chunksize=chunksize)


def read_anchovy_csv(path: str) -> pd.DataFrame:
    """Read an anchovy CSV produced by write_anchovy_csv."""
    return pd.read_csv(path)


# --------------------------------------------------------------------------- #
# FASTA description parsing (the shared Python/R contract)
# --------------------------------------------------------------------------- #
def parse_description(description: str) -> dict[str, float | None]:
    """Parse coverage and length out of a consensus FASTA description line.

    THE SHARED CONTRACT. The original had two divergent parsers for this same
    string: ConsensusTool.readSeqs did description.split(" ")[2].split("coverage:")
    and the R extractCoverageStats did tstrsplit + a split on "coverage:" -- and
    the R version had a copy-paste bug parsing 'length' out of the coverage token.
    Both languages should now agree on THIS definition of the format.

    Format assumed (unchanged from the original convention):
        "<barcode> <ref> coverage:<NN> length:<MM> ..."
    i.e. whitespace-separated tokens, with 'coverage:' and 'length:' prefixes
    somewhere among them.

    Returns a dict with 'coverage' and 'length' as floats, or None if a key is
    absent. Returning None rather than raising preserves the tolerant spirit of
    the original while being explicit about a missing field.
    """
    coverage: float | None = None
    length: float | None = None
    for token in description.split():
        if token.startswith(DESCRIPTION_COVERAGE_KEY):
            try:
                coverage = float(token[len(DESCRIPTION_COVERAGE_KEY):])
            except ValueError:
                coverage = None
        elif token.startswith(DESCRIPTION_LENGTH_KEY):
            try:
                length = float(token[len(DESCRIPTION_LENGTH_KEY):])
            except ValueError:
                length = None
    return {"coverage": coverage, "length": length}
