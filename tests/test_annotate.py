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

    # The amino-acid label too -- genotypeID here, genotypeName in the R. Only
    # the column name differs; the values must still match. It was NOT compared
    # for a long time, and a divergence hid there: the port joined the
    # per-mutation names in the order the tokens appeared in the genotype
    # string, while the R joined them in genome-position order, so "13A_8T"
    # came out "R5S_D3V" against the R's "D3V_R5S". Only multi-mutation
    # genotypes were affected, which is why nothing noticed. It is a
    # user-visible label -- it names nodes in the network files -- so it is
    # worth holding to the golden.
    def labels(df, column):
        d = df[["genotype", column]].drop_duplicates()
        return {row["genotype"]: row[column]
                for _, row in d.iterrows()
                # reference cells have no substitutions; the two sides render
                # that as "" and NaN respectively, which is not a disagreement.
                if isinstance(row[column], str) and row[column]}

    assert labels(got, "genotypeID") == labels(exp, "genotypeName"), (
        f"genotype labels differ:\n got: {labels(got, 'genotypeID')}"
        f"\n exp: {labels(exp, 'genotypeName')}"
    )
    assert "genotypeName" not in got.columns, (
        "genotypeName is back; the column holds an amino-acid translation of "
        "the nucleotide genotype, not the genotype's name, and reading it as a "
        "name is what put an amino-acid-grouped frequency on a node")

    # Both frequencies. THE NAMES ARE SWAPPED RELATIVE TO THE R, deliberately,
    # so this comparison is crosswise:
    #
    #     R genoFreq   grouped on the amino-acid name  ->  port idFreq
    #     R haploFreq  grouped on the nucleotide one   ->  port genoFreq
    #
    # The numbers are unchanged; what changed is which one is called genoFreq,
    # and therefore which one _genotypeNodes.csv offers Cytoscape as the node
    # frequency. A node is one nucleotide genotype, so it has to be the
    # genotype-grouped figure. See annotate.py.
    def freqs(df, column):
        d = df[["genotype", column]].drop_duplicates()
        return {row["genotype"]: round(row[column], 9) for _, row in d.iterrows()}

    assert freqs(got, "genoFreq") == freqs(exp, "haploFreq"), (
        f"nucleotide-genotype frequencies differ:\n got: {freqs(got, 'genoFreq')}"
        f"\n exp: {freqs(exp, 'haploFreq')}")
    assert freqs(got, "idFreq") == freqs(exp, "genoFreq"), (
        f"amino-acid-ID frequencies differ:\n got: {freqs(got, 'idFreq')}"
        f"\n exp: {freqs(exp, 'genoFreq')}")
    assert "haploFreq" not in got.columns


def test_network_matches_r_golden(tmp_path):
    """Python network CSVs reproduce the R epistatic/genotype network output.

    Asks for self_edges=True explicitly. The default now drops a genotype's edge
    to itself from the epistatic network, since Cytoscape draws each as a loop on
    the node and no node is lost by removing them. That is a deliberate
    divergence from the R, so the R-fidelity proof has to opt back in rather than
    quietly weaken -- see test_self_edges_dropped_by_default below for the new
    default's own test.
    """
    import pathlib
    data = pathlib.Path(__file__).parent / "data" / "annot"
    fixture_csv = data / "filtConsensus.csv"
    ref = data / "reference.txt"
    epi_golden = data / "annot_epistaticNetwork.csv"
    for p in (fixture_csv, ref, epi_golden):
        if not p.exists():
            pytest.skip(f"missing {p.name}; run make_annot_fixtures.py + R analysis first.")

    out_prefix = str(tmp_path / "py")
    run(str(fixture_csv), str(ref), out_prefix, network=True, self_edges=True)

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
# Self-edges: dropped from the epistatic network, kept in the genotype network
# --------------------------------------------------------------------------- #
def _fixture(request):
    data = request.path.parent / "data" / "annot"
    for name in ("filtConsensus.csv", "reference.txt"):
        if not (data / name).exists():
            pytest.skip(f"missing {name}; run tests/make_annot_fixtures.py first.")
    return data


def test_self_edges_dropped_by_default(request, tmp_path):
    """The default epistatic network has no genotype linked to itself."""
    data = _fixture(request)
    out = str(tmp_path / "py")
    run(str(data / "filtConsensus.csv"), str(data / "reference.txt"),
        out, network=True)

    edges = pd.read_csv(f"{out}_epistaticNetwork.csv")
    loops = edges[edges["source"] == edges["target"]]
    assert loops.empty, f"self-edges survived: {loops.to_dict('records')}"


def test_dropping_self_edges_loses_no_genotype(request, tmp_path):
    """The reason dropping them is safe HERE: every genotype is still present.

    Each one is reachable by a step edge or by its edge to the reference, so the
    epistatic network keeps the same node set either way. If that ever stops
    being true, this fails rather than silently shrinking someone's figure.
    """
    data = _fixture(request)
    kept, dropped = (str(tmp_path / "kept"), str(tmp_path / "dropped"))
    run(str(data / "filtConsensus.csv"), str(data / "reference.txt"),
        kept, network=True, self_edges=True)
    run(str(data / "filtConsensus.csv"), str(data / "reference.txt"),
        dropped, network=True)

    def nodes(prefix):
        e = pd.read_csv(f"{prefix}_epistaticNetwork.csv")
        return set(e["source"]) | set(e["target"])

    assert nodes(kept) == nodes(dropped)


def test_genotype_network_keeps_its_self_edges(request, tmp_path):
    """all_entries must NOT drop them -- there they are load-bearing.

    A genotype sharing no mutation with any other appears in the genotype
    network only as its own self-edge. On this fixture that is 19T, 6G and the
    reference: dropping self-edges there would delete 3 of 5 genotypes.
    """
    data = _fixture(request)
    out = str(tmp_path / "py")
    run(str(data / "filtConsensus.csv"), str(data / "reference.txt"),
        out, network=True)                       # default: drop where safe

    edges = pd.read_csv(f"{out}_genotypeNetwork.csv")
    loops = edges[edges["source"] == edges["target"]]
    assert not loops.empty, "genotype network lost its self-edges"

    present = set(edges["source"]) | set(edges["target"])
    for singleton in ("19T", "6G", "reference"):
        assert singleton in present, (
            f"{singleton} has no non-self edge and vanished from the genotype "
            f"network")
