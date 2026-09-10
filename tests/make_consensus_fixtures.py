"""
make_consensus_fixtures.py -- generate the consensus-stage fixture.

Run once:
    python tests/make_consensus_fixtures.py

WHAT ConsensusTool.py EXPECTS (learned by reading the code)
-----------------------------------------------------------
- An *_allConsensus.fasta with headers formatted so readSeqs can parse them:
      >ID <token> coverage:<NN> length:<MM>
  because readSeqs does description.split(" ")[2].split("coverage:")[1] for depth
  and [3].split("length:")[1] for length. So coverage must be the 3rd whitespace
  token and length the 4th.
- All sequences the SAME LENGTH. genotypeSummary builds a rectangular column
  matrix via np.transpose; ragged sequences would fail or misbehave.
- selectSeqs keeps depth > depthMin (default 10) AND < 3 gaps in the target
  region, then trims each sequence to [start:end].
- Variant sites are computed on the TRIMMED sequences, and reported positions
  are 1-based within the trimmed region.

THE SCENARIO (predicted exactly; see the design notes in the PR discussion)
---------------------------------------------------------------------------
Base sequence length 30, ORF region start=3 end=27 (trimmed length 24).
Five cells, all coverage > 10 so none are filtered:
    AAACCCAAGAAACACT  cov=37  -> genotype ""        (matches consensus)
    AAACCCAAGAAACCAT  cov=42  -> genotype "3T"
    AAACCCAAGAAACCCA  cov=25  -> genotype "8A"
    AAACCCAAGAAACCTG  cov=31  -> genotype "3T_18G"
    AAACCCAAGAAACGGG  cov=55  -> genotype ""        (reinforces majority)
Computed consensus (trimmed) = "TACGTACGTACGTACGTACGTACG".

This scenario deliberately does NOT trip the coverage/gap filters (all cells
pass), so the golden captures the clean happy path. Filtering and gap edge cases
are covered by unit tests instead, where inputs are built to trip them.
"""

from __future__ import annotations

from pathlib import Path

DATA = Path(__file__).parent / "data"

# Region of interest passed to the tool as start/end.
ORF_START = 3
ORF_END = 27

# Deterministic base sequence of length 30.
BASE = ("ACGT" * 8)[:30]

# (barcode_id, coverage, [(pos0based, newbase), ...])
CELLS = [
    ("AAACCCAAGAAACACT", 37, []),
    ("AAACCCAAGAAACCAT", 42, [(5, "T")]),
    ("AAACCCAAGAAACCCA", 25, [(10, "A")]),
    ("AAACCCAAGAAACCTG", 31, [(5, "T"), (20, "G")]),
    ("AAACCCAAGAAACGGG", 55, []),
]


def mutate(seq: str, changes) -> str:
    chars = list(seq)
    for pos, base in changes:
        chars[pos] = base
    return "".join(chars)


def write_all_consensus() -> None:
    """Write NAME_allConsensus.fasta with tool-parseable headers.

    Header layout: ">{id} ref coverage:{cov} length:{len}"
    so that split(" ")[2] == 'coverage:NN' and [3] == 'length:MM'.
    """
    records = []
    for cid, cov, changes in CELLS:
        seq = mutate(BASE, changes)
        header = f">{cid} ref coverage:{cov} length:{len(seq)}"
        records.append(f"{header}\n{seq}")

    # The tool derives sibling output names by replacing '_allConsensus.fasta',
    # so the input MUST end with that suffix.
    path = DATA / "consensus_test_allConsensus.fasta"
    path.write_text("\n".join(records) + "\n")
    print(f"wrote {path} ({len(CELLS)} sequences, length {len(BASE)})")
    print(f"run the tool with start={ORF_START} end={ORF_END}")


if __name__ == "__main__":
    DATA.mkdir(parents=True, exist_ok=True)
    write_all_consensus()
