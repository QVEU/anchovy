"""
make_consensus_golden.py -- freeze the original ConsensusTool.py output as golden.

Run once, BEFORE migrating ConsensusTool.py:
    python tests/make_consensus_golden.py legacy/ConsensusTool_original.py

(Adjust the path to wherever the original ConsensusTool.py lives -- e.g. if you
moved it to legacy/ alongside the original anchovy.py.)

WHAT IT CAPTURES
----------------
The original writes three siblings next to the input, replacing
'_allConsensus.fasta' with:
    _consensus_reference.txt   (the computed consensus sequence)
    _filtConsensus.csv         (CBC_ID, genotype, sequence, description table)
    _filtConsensus.fasta       (filtered sequences)
We freeze the reference text and the CSV -- those are what the R stage consumes
and what the migrated genotype_summary() must reproduce on the no-reference path.

VERIFICATION
------------
Unlike the extract golden (which we trusted because we had no independent
prediction), here we DID predict the exact output by hand. So this script checks
the frozen values against those predictions and refuses to write a golden that
disagrees. A golden you've verified against an independent calculation is far
stronger than one you merely captured.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd

DATA = Path(__file__).parent / "data"
GOLDEN = DATA / "golden"

ORF_START = "3"
ORF_END = "27"

# Predicted values (from the design calculation). The freeze must match these.
EXPECTED_CONSENSUS = "TACGTACGTACGTACGTACGTACG"
EXPECTED_GENOTYPES = {
    "AAACCCAAGAAACACT": "",
    "AAACCCAAGAAACCAT": "3T",
    "AAACCCAAGAAACCCA": "8A",
    "AAACCCAAGAAACCTG": "3T_18G",
    "AAACCCAAGAAACGGG": "",
}


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: python tests/make_consensus_golden.py /path/to/ConsensusTool.py")
        return 2

    script = sys.argv[1]
    infile = DATA / "consensus_test_allConsensus.fasta"
    if not infile.exists():
        print(f"Fixture missing: {infile}\nRun: python tests/make_consensus_fixtures.py")
        return 1

    GOLDEN.mkdir(parents=True, exist_ok=True)

    cmd = [sys.executable, script, str(infile), ORF_START, ORF_END]
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print("STDERR:\n", result.stderr)
        return 1

    ref_path = DATA / "consensus_test_consensus_reference.txt"
    csv_path = DATA / "consensus_test_filtConsensus.csv"
    if not (ref_path.exists() and csv_path.exists()):
        print("Expected tool outputs not found next to the fixture.")
        return 1

    consensus = ref_path.read_text().strip()
    df = pd.read_csv(csv_path)

    # --- Verify against hand-predicted values before freezing --------------- #
    problems = []
    if consensus != EXPECTED_CONSENSUS:
        problems.append(f"consensus mismatch:\n  got: {consensus}\n  exp: {EXPECTED_CONSENSUS}")

    got_genos = dict(zip(df["CBC_ID"], df["genotype"].fillna("")))
    for cid, exp in EXPECTED_GENOTYPES.items():
        got = got_genos.get(cid, "<missing>")
        if got != exp:
            problems.append(f"genotype[{cid}]: got '{got}', expected '{exp}'")

    if problems:
        print("\nFREEZE REFUSED -- output disagrees with hand prediction:")
        for p in problems:
            print("  -", p)
        print("\nEither the fixture/prediction is wrong, or the tool behaves\n"
              "differently than we modeled. Resolve before freezing.")
        return 1

    # --- Freeze --------------------------------------------------------------#
    (GOLDEN / "consensus_reference.txt").write_text(consensus + "\n")
    # Keep only the columns the R stage relies on, canonicalized by CBC_ID.
    keep = df[["CBC_ID", "genotype"]].copy()
    keep["genotype"] = keep["genotype"].fillna("")
    keep = keep.sort_values("CBC_ID").reset_index(drop=True)
    keep.to_csv(GOLDEN / "consensus_genotypes.csv", index=False)

    print("\nVerified against hand prediction. Wrote:")
    print(f"  {GOLDEN / 'consensus_reference.txt'}")
    print(f"  {GOLDEN / 'consensus_genotypes.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
