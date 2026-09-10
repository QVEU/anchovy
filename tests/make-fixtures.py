"""
make_fixtures.py -- generate a tiny, valid SAM + whitelist for golden testing.

Run once (or whenever you intentionally change the test scenario):
    python tests/make_fixtures.py

WHY GENERATE RATHER THAN HAND-WRITE
-----------------------------------
A SAM with a correct header, valid CIGARs, and reads that actually contain the
10X signature + real barcodes is tedious and easy to get subtly wrong by hand.
Generating it makes the scenario explicit and reproducible: you can read this
file and know exactly what the pipeline is being tested on.

WHAT THE SCENARIO IS
--------------------
- A whitelist of a few known 16-nt barcodes.
- A handful of reads, each built as: [random padding] + signature-with-a-known
  barcode-and-UMI spliced in + [random padding], then mapped to a single
  reference contig. Each barcode appears on enough reads that it survives the
  downstream min-reads-per-cell filter.
- One read deliberately unmapped, and one deliberately too short, to exercise
  the filtering paths.

This is NOT trying to be biologically realistic. It's the smallest input that
drives every branch of anchovy.py's logic so the golden output is meaningful.
"""

from __future__ import annotations

import random
from pathlib import Path

# Deterministic: a fixed seed means the generated fixture is identical every run,
# which is essential -- the golden output must correspond to a fixed input.
random.seed(1234)

DATA = Path(__file__).parent / "data"

# The default 10X signature (matches config.DEFAULT_TENX_SIGNATURE).
SIGNATURE = ("CTACACGACGCTCTTCCGATCT"
             "NNNNNNNNNNNNNNNNNNNNNNNNNN"
             "TTTCTTATAT")

# A few known 16-nt cell barcodes.
BARCODES = [
    "AAACCCAAGAAACACT",
    "AAACCCAAGAAACCAT",
    "AAACCCAAGAAACCCA",
    "AAACCCAAGAAACCTG",
]

REFERENCE_NAME = "testref"
REFERENCE_LENGTH = 2000
NT = "ACGT"


def random_seq(n: int) -> str:
    return "".join(random.choice(NT) for _ in range(n))


def build_signature_instance(barcode: str) -> str:
    """Splice a real barcode + random UMI into the signature's N-region.

    The signature's run of N's is the CBC+UMI slot. We fill the first 16 with the
    barcode and the rest with a random UMI, so anchovy can recover the barcode by
    distance matching.
    """
    n_start = SIGNATURE.index("N")
    n_end = SIGNATURE.rindex("N") + 1
    n_len = n_end - n_start
    umi_len = n_len - len(barcode)
    filled = barcode + random_seq(max(0, umi_len))
    return SIGNATURE[:n_start] + filled[:n_len] + SIGNATURE[n_end:]


def build_read(barcode: str):
    """Build a read and the CIGAR that makes anchovy able to find the signature.

    KEY INVARIANT (discovered via golden testing):
    anchovy's poolBlocks searches a window UPSTREAM of the read's 'offset', where
    offset = readLen - clipReadLen comes from the CIGAR. For a pure-match CIGAR
    (e.g. "100M") offset is 0, so the search window seq[max(0,offset-200):offset]
    is empty and NO signature is ever found. The real pipeline's reads have
    LEADING SOFT CLIPS (the cDNA/adapter portion that doesn't map to the viral
    reference), which is what produces a nonzero offset pointing just past the
    signature. We reproduce that here.

    Layout of the read:
        [ left padding + signature-with-barcode ]  <- soft-clipped (S)
        [ right region ]                           <- matched (M)
    The soft-clip length = everything up to and including the signature, so
    offset lands at the end of the signature and the upstream window covers it.

    Returns (sequence, cigar_string).
    """
    left = random_seq(random.randint(40, 80))
    sig = build_signature_instance(barcode)
    right = random_seq(random.randint(120, 200))

    seq = left + sig + right
    soft_clip = len(left) + len(sig)   # clip through the end of the signature
    matched = len(right)
    cigar = f"{soft_clip}S{matched}M"
    return seq, cigar


def write_whitelist() -> None:
    path = DATA / "whitelist.txt"
    path.write_text("\n".join(BARCODES) + "\n")
    print(f"wrote {path} ({len(BARCODES)} barcodes)")


def write_sam() -> None:
    """Write a minimal valid SAM: one @SQ header line + read records."""
    lines = [f"@SQ\tSN:{REFERENCE_NAME}\tLN:{REFERENCE_LENGTH}"]

    read_id = 0
    # Several reads per barcode so each survives min-reads-per-cell filtering.
    for barcode in BARCODES:
        for _ in range(6):
            seq, cigar = build_read(barcode)
            # SAM fields: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL
            lines.append("\t".join([
                f"read{read_id}", "0", REFERENCE_NAME, "1", "60",
                cigar, "*", "0", "0", seq, "I" * len(seq),
            ]))
            read_id += 1

    # One unmapped read (FLAG 4, RNAME '*') to exercise the unmapped filter.
    unmapped, _ = build_read(BARCODES[0])
    lines.append("\t".join([
        f"read{read_id}", "4", "*", "0", "0",
        "*", "*", "0", "0", unmapped, "I" * len(unmapped),
    ]))
    read_id += 1

    # One very short read to exercise the length filter.
    short = random_seq(20)
    lines.append("\t".join([
        f"read{read_id}", "0", REFERENCE_NAME, "1", "60",
        f"{len(short)}M", "*", "0", "0", short, "I" * len(short),
    ]))

    path = DATA / "test.sam"
    path.write_text("\n".join(lines) + "\n")
    print(f"wrote {path} ({len(lines)} lines)")


if __name__ == "__main__":
    DATA.mkdir(parents=True, exist_ok=True)
    write_whitelist()
    write_sam()
    print("\nFixtures generated. Next: generate golden output with")
    print("    python tests/make_golden.py")
