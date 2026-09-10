"""
make_golden.py -- run the CURRENT anchovy.py and freeze its output as golden.

Run once, BEFORE migrating anchovy.py into the package:
    python tests/make_golden.py /path/to/current/anchovy.py

WHAT "GOLDEN" MEANS
-------------------
A characterization (a.k.a. golden-master) test captures what the code CURRENTLY
does -- not what it "should" do -- and asserts future versions still do the same.
It's the safety net that lets us refactor aggressively: if migrated code produces
byte-for-byte the same canonicalized output, behavior was preserved. If it
differs, we look immediately, before the change compounds.

This deliberately runs the ORIGINAL script as a black box (via subprocess) rather
than importing it, because at this stage the original is loose top-level script
code, not an importable module. We capture its real output on the fixture input.

INPUT FORM NOTE
---------------
anchovy.py's loadSAM opens the file with pysam in "rb" (BAM) mode while also
text-parsing it. Depending on your pysam/htslib build, the pysam pass may require
an actual BAM. If a plain .sam fails, convert first:
    samtools view -bS tests/data/test.sam > tests/data/test.bam
and point this script at the .bam. This mirrors the real upstream flow
(minimap2 -> SAM -> samtools) and previews why Phase 4 pins samtools via conda.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd

DATA = Path(__file__).parent / "data"
GOLDEN = DATA / "golden"

SIGNATURE = ("CTACACGACGCTCTTCCGATCT"
             "NNNNNNNNNNNNNNNNNNNNNNNNNN"
             "TTTCTTATAT")


def canonicalize(df: pd.DataFrame) -> pd.DataFrame:
    """Put an anchovy output frame into a canonical, comparable form.

    Golden comparison must be robust to nondeterminism that doesn't reflect a
    real behavior change:
      - multiprocessing can affect row order -> sort rows deterministically
      - a leading unnamed index column from to_csv/read_csv -> drop it
    We sort by a stable set of columns and reset the index so two runs that
    computed the SAME assignments compare equal regardless of row order.
    """
    # Drop the unnamed index column pandas writes/reads, if present.
    df = df.loc[:, [c for c in df.columns if not str(c).startswith("Unnamed")]]
    # Sort by read then CBC (stable, present in every anchovy output).
    sort_cols = [c for c in ["read", "CBC", "UMI"] if c in df.columns]
    df = df.sort_values(sort_cols).reset_index(drop=True)
    return df


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: python tests/make_golden.py /path/to/anchovy.py [input.bam]")
        return 2

    script = sys.argv[1]
    # Allow overriding the input (e.g. a converted .bam); default to the .sam.
    sam_input = sys.argv[2] if len(sys.argv) > 2 else str(DATA / "test.sam")
    whitelist = str(DATA / "whitelist.txt")

    GOLDEN.mkdir(parents=True, exist_ok=True)

    # Run the original script exactly as a user would.
    cmd = [sys.executable, script, sam_input, whitelist, SIGNATURE]
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print("STDERR:\n", result.stderr)
        print("\nThe original script failed. If the error mentions the SAM/BAM,")
        print("convert to BAM and re-run (see the module docstring).")
        return 1

    # anchovy writes <input>_anchovy.csv next to the input.
    produced = Path(sam_input).with_suffix("").as_posix()
    produced = Path(sam_input.replace(".sam", "_anchovy.csv")
                    .replace(".bam", "_anchovy.csv"))
    if not produced.exists():
        print(f"Expected output not found: {produced}")
        return 1

    df = canonicalize(pd.read_csv(produced))
    golden_path = GOLDEN / "extract_anchovy.csv"
    df.to_csv(golden_path, index=False)
    print(f"\nWrote golden output: {golden_path} ({len(df)} rows)")
    print("Commit tests/data/ so the fixture and golden output are version-controlled.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
