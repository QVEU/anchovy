"""
test_annotate.py -- unit + golden tests for the annotation stage.

Two guarantees, as everywhere:
  - unit tests: hand-verified amino-acid calls on known codons (fast, precise)
  - golden test: run() output vs the frozen R reference (proves equivalence to
    the original R analysis)
"""

from __future__ import annotations

import pandas as pd
import pytest

from anchovy.annotate import codon, annotate_mutation, haplo_analysis, run

REFERENCE = "ATGAAAGATTTTCGTGGGCATACTTGGGAA"   # MKDFRGHTWE


# --------------------------------------------------------------------------- #
# codon / annotate_mutation -- hand-verified against the standard genetic code
# --------------------------------------------------------------------------- #
def test_codon_first_residue():
    c = codon(REFERENCE, 1)          # ATG -> M, residue 1
    assert c["codon"] == "ATG" and c["resPos"] == 1 and c["AA"] == "M"


def test_codon_residue_five():
    c = codon(REFERENCE, 13)         # CGT -> R, residue 5
    assert c["codon"] == "CGT" and c["resPos"] == 5 and c["AA"] == "R"


def test_annotate_synonymous():
    a = annotate_mutation("6G", REFERENCE)   # AAA->AAG, K2K
    assert a["subName"] == "K2K" and a["subClass"] == "Syn"


def test_annotate_nonsynonymous_R5S():
    a = annotate_mutation("13A", REFERENCE)  # CGT->AGT, R5S
    assert a["ref"]["codon"] == "CGT" and a["mut"]["codon"] == "AGT"
    assert a["subName"] == "R5S" and a["subClass"] == "Non-Syn"


def test_annotate_nonsynonymous_D3V():
    a = annotate_mutation("8T", REFERENCE)   # GAT->GTT, D3V
    assert a["subName"] == "D3V" and a["subClass"] == "Non-Syn"


def test_annotate_nonsynonymous_H7Y():
    a = annotate_mutation("19T", REFERENCE)  # CAT->TAT, H7Y
    assert a["subName"] == "H7Y" and a["subClass"] == "Non-Syn"


def test_annotate_reference_sentinel():
    a = annotate_mutation("reference", REFERENCE)
    assert a["subName"] == "ref" and a["subClass"] == "WT"


# --------------------------------------------------------------------------- #
# haplo_analysis -- frequencies
# --------------------------------------------------------------------------- #
def test_haplo_analysis_frequencies():
    cons = pd.DataFrame({
        "CBC_ID": ["c1", "c2", "c3", "c4"],
        "genotype": ["", "13A", "13A", "6G"],
    })
    t = haplo_analysis(cons)
    # DEPTH = 4; 13A appears in 2 cells -> count 2, freq 0.5
    freq_13A = t.loc[t["mutants"] == "13A", "freq"].iloc[0]
    assert freq_13A == 0.5
    # one reference cell present
    assert (t["genotype"] == "reference").sum() == 1


# --------------------------------------------------------------------------- #
# GOLDEN: run() vs the frozen R output on the annotation fixture
# --------------------------------------------------------------------------- #
def test_annotate_matches_r_golden(tmp_path):
    """Python run() reproduces the R _annot_v3.csv annotation calls."""
    import pathlib
    data = pathlib.Path(__file__).parent / "data" / "annot"
    fixture_csv = data / "filtConsensus.csv"
    ref = data / "reference.txt"
    golden = data / "annot_annot_v3.csv"
    for p in (fixture_csv, ref, golden):
        if not p.exists():
            pytest.skip(f"missing {p.name}; run make_annot_fixtures.py + R analysis first.")

    out_prefix = str(tmp_path / "py")
    result = run(str(fixture_csv), str(ref), out_prefix, network=True)

    got = result["annot"]
    exp = pd.read_csv(golden)

    # Compare the annotation calls per mutation token (the analytical payload).
    # Build {mutants -> (subName, subClass)} from each and compare.
    def calls(df):
        d = df[df["mutants"].notna() & (df["mutants"] != "")]
        return {row["mutants"]: (row["subName"], row["subClass"])
                for _, row in d.iterrows()}

    assert calls(got) == calls(exp), (
        f"annotation calls differ:\n got: {calls(got)}\n exp: {calls(exp)}"
    )


def test_network_matches_r_golden(tmp_path):
    """Python network CSVs reproduce the R epistatic/genotype network output."""
    import pathlib
    data = pathlib.Path(__file__).parent / "data" / "annot"
    fixture_csv = data / "filtConsensus.csv"
    ref = data / "reference.txt"
    epi_golden = data / "annot_epistaticNetwork.csv"
    for p in (fixture_csv, ref, epi_golden):
        if not p.exists():
            pytest.skip(f"missing {p.name}; run make_annot_fixtures.py + R analysis first.")

    out_prefix = str(tmp_path / "py")
    run(str(fixture_csv), str(ref), out_prefix, network=True)

    # Compare numerically (not as text) so int/float rendering can't cause a
    # spurious mismatch; sort rows so ordering differences don't either -- the
    # SET of edges and their values is what matters. Note the source-node column
    # is "source" (renamed from the R's "genotype" for Cytoscape auto-detection).
    got = pd.read_csv(f"{out_prefix}_epistaticNetwork.csv")
    exp = pd.read_csv(epi_golden)

    sort_cols = ["source", "target", "overlap", "mutNumSource", "mutNumTarget"]
    got_s = got.sort_values(sort_cols).reset_index(drop=True)
    exp_s = exp.sort_values(sort_cols).reset_index(drop=True)

    pd.testing.assert_frame_equal(got_s, exp_s, check_dtype=False)
