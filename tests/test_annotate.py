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


# --------------------------------------------------------------------------- #
# Region-aware path: single-CDS GFF reproduces the legacy annotation calls,
# and additionally writes the long-format region table.
# --------------------------------------------------------------------------- #
def test_region_aware_single_cds_matches_legacy_calls(tmp_path):
    """With a single CDS at position 1, the region-aware path fills the legacy
    annotation columns with the same subName/subClass the frame-1 path produced."""
    import pathlib
    data = pathlib.Path(__file__).parent / "data" / "annot"
    fixture_csv = data / "filtConsensus.csv"
    ref = data / "reference.txt"
    for p in (fixture_csv, ref):
        if not p.exists():
            pytest.skip(f"missing {p.name}")

    # Single-CDS GFF spanning the whole 30 nt reference, frame from position 1.
    gff = tmp_path / "single_cds.gff3"
    gff.write_text("ref\ttest\tCDS\t1\t30\t.\t+\t0\tID=cds;Name=cds\n")

    out_prefix = str(tmp_path / "py")
    result = run(str(fixture_csv), str(ref), out_prefix, network=False, gff=str(gff))

    # The legacy annotation columns must match the hand-verified calls.
    got = result["annot"]
    def calls(df):
        d = df[df["mutants"].notna() & (df["mutants"] != "")]
        return {r["mutants"]: (r["subName"], r["subClass"]) for _, r in d.iterrows()}
    c = calls(got)
    assert c["6G"] == ("K2K", "Syn")
    assert c["13A"] == ("R5S", "Non-Syn")
    assert c["8T"] == ("D3V", "Non-Syn")
    assert c["19T"] == ("H7Y", "Non-Syn")

    # The long-format region table was written and has a row per variant.
    region_csv = pathlib.Path(result["written"]["regions"])
    assert region_csv.exists()
    rt = pd.read_csv(region_csv)
    assert set(rt["mutants"]) == {"6G", "8T", "13A", "19T"}
    # every row is annotated against the one CDS region
    assert set(rt["region"]) == {"cds"}


def test_region_aware_noncoding_annotation(tmp_path):
    """A mutation outside the CDS gets a non-coding region row (no amino acid)."""
    import pathlib
    data = pathlib.Path(__file__).parent / "data" / "annot"
    fixture_csv = data / "filtConsensus.csv"
    ref = data / "reference.txt"
    for p in (fixture_csv, ref):
        if not p.exists():
            pytest.skip(f"missing {p.name}")

    # CDS covers 1-12; declare a 3' UTR over 13-30 so mutation 19T lands in it.
    gff = tmp_path / "two_region.gff3"
    gff.write_text(
        "ref\tt\tCDS\t1\t12\t.\t+\t0\tID=cds;Name=cds\n"
        "ref\tt\tthree_prime_UTR\t13\t30\t.\t+\t.\tID=3utr;Name=3UTR\n"
    )
    out_prefix = str(tmp_path / "py")
    result = run(str(fixture_csv), str(ref), out_prefix, network=False, gff=str(gff))

    rt = pd.read_csv(pathlib.Path(result["written"]["regions"]))
    # 19T (genome pos 19) is in the 3' UTR -> non-coding, no amino acid.
    utr = rt[(rt["mutants"] == "19T") & (rt["region"] == "3UTR")]
    assert len(utr) == 1
    assert utr.iloc[0]["region_type"] == "non-coding"
    assert pd.isna(utr.iloc[0]["wt_aa"])
    assert utr.iloc[0]["mutation_id"] == "3UTR:C19T"
