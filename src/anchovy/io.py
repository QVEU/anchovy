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

import glob
import os
import sys
from pathlib import Path

import pandas as pd
import pysam

from anchovy.schema import (
    SamColumns,
    AnchovyColumns,
    DESCRIPTION_COVERAGE_KEY,
    DESCRIPTION_LENGTH_KEY,
    SIGNATURE_BARCODE_LEN,
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

    ONE STREAMING PASS. The original loadSAM read the file twice -- a manual
    text parse of the 11 core fields, then a pysam pass for CIGAR accounting --
    and built a list of every non-header line before filtering anything. That
    is bounded by the WHOLE FILE rather than by the reads that survive, which on
    a real run is the difference that matters: of 3,457,184 records in a PacBio
    run, 402,303 reached the next stage. The other 88% were parsed, stored, and
    copied into a DataFrame before being discarded.

    pysam already parses the record and exposes the CIGAR, so a single
    `for read in fp:` gets everything both passes got. Filtering inside the loop
    means a discarded read is never retained at all.

    WHICH READS SURVIVE IS UNCHANGED. The old code filtered rows on
    `RNAME != "*"` (the text pass) while computing CIGARs for
    `not read.is_unmapped` (the pysam pass), then zipped the two together. Those
    are different predicates: a read that is flagged unmapped but still carries
    a reference name -- legal SAM, an unmapped mate placed at its partner's
    locus -- is kept by the first and dropped by the second. When they disagree
    the old code raised `ValueError: Length of values (N) does not match length
    of index (M)` from pandas, so it could not silently misalign, but it also
    could not proceed. This keeps the row predicate, `reference_name is not
    None`, and computes the CIGAR for exactly those rows, so the two can no
    longer disagree.

    Args:
        path: path to the input SAM file.
        min_read_length: keep reads with length strictly greater than this.
            (Callers pass config.extract.effective_min_read_length(), which
            defaults to the signature length -- the original's `minL = quL`.)
    """
    rows: list[list[str]] = []
    read_lens: list[int] = []
    clip_lens: list[int] = []
    offsets: list[int] = []
    lengths: list[int] = []
    n_records = 0

    print("Parsing Cigars...")
    # check_sq=False so a SAM with no @SQ header still opens; the old text pass
    # never looked at the header at all.
    with pysam.AlignmentFile(path, "rb", check_sq=False) as handle:
        for read in handle:
            n_records += 1

            # The old text filter, expressed on the parsed record: pysam reports
            # RNAME "*" as reference_name None. Applied before anything is kept,
            # so a rejected read never occupies memory.
            if read.reference_name is None:
                continue

            # len() of the SEQ field, matching the old df[SEQ].apply(len); a
            # read with SEQ "*" measured 1 there, and query_sequence is None.
            seq = read.query_sequence or "*"
            if len(seq) <= min_read_length:
                continue

            # to_string() re-renders the record as its SAM text line, so the 11
            # fields are the same STRINGS the text parse produced -- same values
            # and same object dtype. Taking pysam's native attributes instead
            # would silently turn FLAG and POS into ints and change the frame.
            rows.append(read.to_string().split("\t")[:len(SamColumns.ORDER)])

            read_len, clip_len, offset = parse_cigar_lengths(read)
            read_lens.append(read_len)
            clip_lens.append(clip_len)
            offsets.append(offset)
            lengths.append(len(seq))
    print("Done.")

    print("Total Candidate Reads: {}".format(n_records))

    df = pd.DataFrame(rows, columns=SamColumns.ORDER)
    df[SamColumns.READ_LEN] = read_lens
    df[SamColumns.CLIP_READ_LEN] = clip_lens
    df[SamColumns.OFFSET] = offsets
    df[SamColumns.LENGTH] = lengths
    return df


def validate_whitelist(barcodes: list[str], path: str = "whitelist") -> None:
    """Check a barcode list is one anchovy can search with (raises, pure).

    WHY THIS IS FATAL RATHER THAN A WARNING. build_barcode_query_blocks builds
    one search template per barcode as

        signature[0:22] + barcode + "N"*umi + signature[-10:]

    so a barcode of the wrong width shifts everything after it. The reads then
    match nothing well -- and unless max_barcode_errors is set, each one is
    still assigned to its NEAREST entry rather than rejected. The run finishes,
    every cell is wrong, and nothing says so. That is the failure this refuses
    to let through.

    THE BARCODE IS 16 nt IN BOTH v2 AND v3 (schema.SIGNATURE_BARCODE_LEN); the
    chemistries differ only in UMI width. So this check is chemistry-
    independent, and a length that is not 16 is wrong whatever the config says.
    """
    if not barcodes:
        raise ValueError(f"{path} has no barcodes in it.")

    lengths = {len(b) for b in barcodes}
    if lengths == {SIGNATURE_BARCODE_LEN}:
        # Right width. A header line would have to be a plausible 16-mer to get
        # here, so only the alphabet is left to check, and one bad entry is
        # enough to look at -- scanning 6.8M barcodes to find the rest is waste.
        bad = next((b for b in barcodes if set(b.upper()) - set("ACGTN")), None)
        if bad is not None:
            raise ValueError(
                f"{path}: {bad!r} is not a barcode -- it is the right length "
                f"but contains something other than A/C/G/T/N.\n"
                f"  A header line or a quoted CSV export will do this.")
        return

    offender = next(b for b in barcodes if len(b) != SIGNATURE_BARCODE_LEN)

    # The GEM suffix first: it is much the likeliest thing to be wrong here and
    # the only one with a one-line fix. Cell Ranger's barcodes.tsv, and
    # anything taken from Seurat's column names, carry it.
    if "-" in offender:
        raise ValueError(
            f"{path}: barcodes carry a GEM suffix, e.g. {offender!r}.\n"
            f"  That suffix is not part of the barcode. Left on, every search "
            f"template is built\n"
            f"  the wrong width and reads are assigned to an arbitrary nearest "
            f"entry rather than\n"
            f"  failing -- the run finishes and every cell is wrong. Strip it:\n"
            f"    cut -d- -f1 {path} > barcodes_stripped.txt")

    if len(lengths) > 1:
        raise ValueError(
            f"{path}: barcodes are not all the same length "
            f"(saw {sorted(lengths)}), e.g. {offender!r}.\n"
            f"  A mixed-width list usually means a header line, a blank line "
            f"mid-file, or two\n"
            f"  files concatenated.")

    raise ValueError(
        f"{path}: barcodes are {offender and len(offender)} nt, but a 10X cell "
        f"barcode is {SIGNATURE_BARCODE_LEN} nt.\n"
        f"  This is NOT a v2/v3 difference -- both use a "
        f"{SIGNATURE_BARCODE_LEN} nt barcode and differ only in UMI width. "
        f"Check\n  you have a barcode list rather than, say, a feature or "
        f"UMI list.")


def read_whitelist(path: str) -> pd.DataFrame:
    """Read a 10X barcode whitelist into a single-column DataFrame.

    Reproduces loadBC: take the first tab-separated field of each line,
    stripped -- then validate, because an unusable whitelist does not fail on
    its own. See validate_whitelist.

    Blank lines are dropped rather than becoming empty barcodes: a trailing
    blank line is the commonest way a hand-edited list acquires one, and an
    empty barcode silently becomes an entry every poor read can match.
    """
    with open(path, "r") as handle:
        barcodes = [line.split("\t")[0].strip() for line in handle]
    barcodes = [b for b in barcodes if b]
    validate_whitelist(barcodes, path)
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

    # Write beside the target, then rename. to_csv streams in chunks straight
    # to the destination, so a run that dies partway -- an OOM kill, a wall
    # clock limit -- leaves a truncated but entirely plausible-looking CSV at
    # exactly the name the next stage reads. Snakemake usually deletes a failed
    # job's outputs, but a SIGKILL gives it no chance to. os.replace within the
    # same directory is atomic on POSIX, so the file is either complete or
    # absent, never half.
    target = Path(path)
    tmp = target.with_name(target.name + ".partial")
    try:
        df.to_csv(tmp, chunksize=chunksize)
        os.replace(tmp, target)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


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


# --------------------------------------------------------------------------- #
# Reference sequence reading
# --------------------------------------------------------------------------- #
def read_reference_sequence(path: str | Path) -> str:
    """Read a reference sequence from a FASTA file or a raw sequence file.

    WHY THIS EXISTS -- it fixes a silent data-corruption hazard. Both callers
    used to do the obvious thing, `open(path).read().strip()`. Handed a FASTA --
    the format every reference genome actually ships in -- that put the header
    line INTO the sequence:

        ">testref\\nAGGGTATCTAAA..."

    Nothing raised. Every genome coordinate simply shifted by the length of the
    header, and reference[i] could land on a newline, so the whole run produced
    confidently wrong calls. A reference that is merely too SHORT at least fails
    loudly with an IndexError; this one failed silently, which is worse.

    Accepts either form:
      - FASTA: a '>' header followed by sequence lines, which may be wrapped.
      - Raw: the bare sequence, optionally wrapped across lines.

    Internal whitespace and line breaks are removed in both cases, so a wrapped
    sequence reads the same as a single-line one. Case is preserved; callers
    apply their own normalization.

    Raises ValueError on an empty file, or on a multi-record FASTA -- the
    pipeline assumes a single reference (see the Snakefile's cell_consensus
    note), so a segmented genome should fail loudly here rather than silently
    annotate everything against whichever record happened to come first.
    """
    text = Path(path).read_text()

    if ">" in text:
        records: list[list[str]] = []
        for line in text.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                records.append([])
            elif records:
                records[-1].append(line)
            else:
                raise ValueError(
                    f"{path}: sequence data appears before the first '>' header.")
        if len(records) > 1:
            raise ValueError(
                f"{path}: expected one sequence, found {len(records)} FASTA "
                f"records. The pipeline assumes a single reference; split the "
                f"file or pass the one record you mean to use.")
        sequence = "".join(records[0])
    else:
        sequence = "".join(text.split())

    if not sequence:
        raise ValueError(f"{path}: no reference sequence found.")
    return sequence


# --------------------------------------------------------------------------- #
# Sample discovery (the workflow's input contract)
# --------------------------------------------------------------------------- #
FASTQ_PATTERNS = ("*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz")


def sample_name_from_fastq(path: str | Path) -> str:
    """A FASTQ path's sample name: the basename with the read suffixes stripped.

    Lives here rather than in the Snakefile so it can be tested, and so there is
    ONE definition. It used to have two -- a shell copy for the fetch/run
    scripts and the workflow's own idea -- which disagreed the moment a run was
    pointed at different reads: the reads were mapped under one name and the
    workflow then looked for the other, found the previous run's output still in
    place, and reported it up to date without failing.
    """
    base = os.path.basename(str(path))
    for ext in (".gz", ".fastq", ".fq"):
        if base.endswith(ext):
            base = base[: -len(ext)]
    return base


def discover_fastqs(input_dir: str | Path) -> dict[str, str]:
    """{sample name: path} for every FASTQ directly in `input_dir`.

    Sorted, so the DAG is stable between runs. Raises ValueError when the
    directory holds no reads, or when two of them reduce to the same sample
    name -- reads.fastq and reads.fq.gz do, and their results would overwrite
    each other silently rather than collide.

    Does not recurse: a nested layout is far more likely to be a mistake about
    what `input_dir` means than an intent to run every FASTQ under a tree.
    """
    found: dict[str, str] = {}
    for pattern in FASTQ_PATTERNS:
        for path in sorted(glob.glob(os.path.join(str(input_dir), pattern))):
            name = sample_name_from_fastq(path)
            if name in found:
                raise ValueError(
                    f"two FASTQs in {input_dir} reduce to the sample name "
                    f"{name!r}:\n  {found[name]}\n  {path}\n"
                    f"Their results would overwrite each other. Rename one.")
            found[name] = path
    if not found:
        raise ValueError(
            f"no FASTQs in input_dir {str(input_dir)!r}.\n"
            f"  Looked for {', '.join(FASTQ_PATTERNS)}. Check the path, and "
            f"that the reads are not\n  still compressed in another format "
            f"(.bz2, .zst) or nested in a subdirectory -- this does not "
            f"recurse.")
    return found
