"""
make_annot_fixtures.py -- fixture for the annotation stage (rich enough to
actually exercise the codon translation, syn/non-syn calls, and network edges).

Run once:
    python tests/make_annot_fixtures.py

WHY THIS FIXTURE IS DESIGNED THE WAY IT IS
------------------------------------------
The consensus fixture produces only one trivial variant genotype -- too thin to
test annotation. This fixture is built to stress every branch of the annotation
logic, with EVERY expected output hand-computed (see the analysis-tool work in
the PR discussion) so the R golden can be verified rather than trusted.

Reference (30 nt, 10 codons):  ATGAAAGATTTTCGTGGGCATACTTGGGAA
Protein:                       M  K  D  F  R  G  H  T  W  E

Cells / genotypes:
    cell01  ""        reference
    cell02  ""        reference
    cell03  13A       R5S   (non-syn)
    cell04  13A       R5S   (shares 13A)
    cell05  13A_8T    R5S + D3V  (multi-mutation; shares 13A -> single-step edge)
    cell06  6G        K2K   (SYNONYMOUS)
    cell07  19T       H7Y   (non-syn)

This gives: reference rows, single + multi mutation genotypes, a mutation shared
across three cells (13A -> real network overlap), and both syn and non-syn calls.

Predicted annotations (verified by hand against the standard genetic code):
    13A: CGT->AGT  R5S  Non-Syn
    8T : GAT->GTT  D3V  Non-Syn
    6G : AAA->AAG  K2K  Syn
    19T: CAT->TAT  H7Y  Non-Syn
"""

from __future__ import annotations

import json
from pathlib import Path

DATA = Path(__file__).parent / "data" / "annot"

REFERENCE = "ATGAAAGATTTTCGTGGGCATACTTGGGAA"   # MKDFRGHTWE

CELLS = [
    ("cell01", ""),
    ("cell02", ""),
    ("cell03", "13A"),
    ("cell04", "13A"),
    ("cell05", "13A_8T"),
    ("cell06", "6G"),
    ("cell07", "19T"),
]

# Hand-verified expected annotations (see docstring). subName = refAA+resPos+mutAA.
EXPECTED_ANNOTATIONS = {
    "13A": {"refCodon": "CGT", "mutCodon": "AGT", "subName": "R5S", "subClass": "Non-Syn"},
    "8T":  {"refCodon": "GAT", "mutCodon": "GTT", "subName": "D3V", "subClass": "Non-Syn"},
    "6G":  {"refCodon": "AAA", "mutCodon": "AAG", "subName": "K2K", "subClass": "Syn"},
    "19T": {"refCodon": "CAT", "mutCodon": "TAT", "subName": "H7Y", "subClass": "Non-Syn"},
}

# Predicted per-mutation counts/frequencies (DEPTH = 7 unique cells).
DEPTH = len(CELLS)
EXPECTED_MUTATION_FREQ = {
    "13A": {"count": 3, "freq": 3 / DEPTH},
    "8T":  {"count": 1, "freq": 1 / DEPTH},
    "6G":  {"count": 1, "freq": 1 / DEPTH},
    "19T": {"count": 1, "freq": 1 / DEPTH},
}

# Predicted per-genotype names and haplotype frequencies.
EXPECTED_GENOTYPES = {
    "":        {"name": "reference", "haploFreq": 2 / DEPTH},
    "13A":     {"name": "R5S",       "haploFreq": 2 / DEPTH},
    "13A_8T":  {"name": "R5S_D3V",   "haploFreq": 1 / DEPTH},
    "6G":      {"name": "K2K",       "haploFreq": 1 / DEPTH},
    "19T":     {"name": "H7Y",       "haploFreq": 1 / DEPTH},
}

# Key network relationship the fixture is designed to produce.
EXPECTED_NETWORK_NOTES = {
    "single_step_edge": "13A -> 13A_8T (overlap 1, differ by adding 8T)",
    "reference_edges": ["13A -> reference", "6G -> reference", "19T -> reference"],
}


def build():
    DATA.mkdir(parents=True, exist_ok=True)

    # Reference sequence file (plaintext, single sequence -- what readReference reads).
    (DATA / "reference.txt").write_text(REFERENCE)

    # filtConsensus.csv: the columns the R analysis reads are CBC_ID and genotype.
    # (sequence/description columns aren't used by haploanalysis, but we include
    # minimal placeholders so the CSV shape matches a real filtConsensus.csv.)
    lines = ["CBC_ID,genotype,sequence,description"]
    for cbc, geno in CELLS:
        lines.append(f"{cbc},{geno},{REFERENCE},{cbc} ref coverage:20 length:30")
    (DATA / "filtConsensus.csv").write_text("\n".join(lines) + "\n")

    # Freeze all hand predictions for the golden-verification step.
    expected = {
        "reference": REFERENCE,
        "protein": "MKDFRGHTWE",
        "depth": DEPTH,
        "annotations": EXPECTED_ANNOTATIONS,
        "mutation_freq": EXPECTED_MUTATION_FREQ,
        "genotypes": EXPECTED_GENOTYPES,
        "network": EXPECTED_NETWORK_NOTES,
    }
    (DATA / "annot_expected.json").write_text(json.dumps(expected, indent=2))

    print(f"wrote {DATA}/reference.txt ({len(REFERENCE)} nt, protein MKDFRGHTWE)")
    print(f"wrote {DATA}/filtConsensus.csv ({len(CELLS)} cells)")
    print(f"wrote {DATA}/annot_expected.json (hand-verified predictions)")
    print("\nNext: run the R analysis to freeze a golden, then verify it against")
    print("these predictions:")
    print(f"  Rscript tests/annotate_analysis.R {DATA}/reference.txt "
          f"{DATA}/filtConsensus.csv {DATA}/annot")


if __name__ == "__main__":
    build()
