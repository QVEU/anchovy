"""
test_consensus.py -- tests for the consensus stage.

TWO KINDS OF GUARANTEE, deliberately separated (see consensus.py header):

1. reference=None reproduces the ORIGINAL's behavior -> GOLDEN tested against the
   frozen, hand-verified output. Proves the migration preserved what worked.

2. reference=<given> is NEW capability (the old else-branch was dead + broken) ->
   UNIT tested on small inputs with known-correct answers. Proves the fix works.

Plus pure unit tests for select_sequences (the coverage/gap filters the golden
fixture deliberately doesn't trip) and parse_consensus_fasta.
"""

from __future__ import annotations

import csv

import pytest

from anchovy.config import ConsensusConfig
from anchovy.consensus import (
    ConsensusRecord,
    parse_consensus_fasta,
    select_sequences,
    genotype_summary,
    run,
)


def _rec(cbc_id, seq, coverage, length=None):
    return ConsensusRecord(
        cbc_id=cbc_id, seq=seq, description=f"{cbc_id} ref coverage:{coverage}",
        coverage=coverage, length=length if length is not None else len(seq),
    )


# --------------------------------------------------------------------------- #
# genotype_summary -- reference=None (the golden-preserved original behavior)
# --------------------------------------------------------------------------- #
def test_genotype_summary_computed_reference_basic():
    # pos 2 varies (A/T/A); majority A -> middle gets "2T", others empty.
    seqs = ["TACG", "TTCG", "TACG"]
    genos, ref = genotype_summary(seqs)
    assert ref == "TACG"
    assert genos == ["", "2T", ""]


def test_genotype_summary_multiple_variant_sites():
    # Mirrors the golden fixture's structure: variant columns at 1-based 3 and 4.
    seqs = ["AAAA", "AATA", "AAAG"]
    genos, ref = genotype_summary(seqs)
    # col2 (0-based): A/A/A -> not variant. col2 idx2: A/T/A variant. idx3: A/A/G variant.
    assert ref == "AAAA"
    assert genos == ["", "3T", "4G"]


def test_genotype_summary_empty_input():
    genos, ref = genotype_summary([])
    assert genos == [] and ref == ""


def test_genotype_summary_no_variation():
    # All identical -> no variant sites -> all empty genotypes.
    genos, ref = genotype_summary(["ACGT", "ACGT", "ACGT"])
    assert genos == ["", "", ""] and ref == "ACGT"


# --------------------------------------------------------------------------- #
# genotype_summary -- reference=<given> (the NEW, previously-broken capability)
# --------------------------------------------------------------------------- #
def test_genotype_summary_supplied_reference_changes_calls():
    # Same sequences, but with an explicit reference the calls differ from the
    # computed-consensus case. Reference says position 1 should be 'G'.
    seqs = ["GACG", "GTCG", "GACG"]
    genos, ref = genotype_summary(seqs, reference="AACG")
    assert ref == "AACG"
    # Variant sites: pos1 (G/G/G) is NOT variant across seqs, so even though all
    # differ from ref 'A' there, no token is emitted (matches original semantics:
    # tokens only at columns that vary across the input). pos2: A/T/A is variant;
    # 'T' differs from ref 'A' -> "2T".
    assert genos == ["", "2T", ""]


def test_genotype_summary_supplied_reference_does_not_crash():
    # The whole point of the fix: this path used to NameError. Now it runs.
    seqs = ["AAA", "ATA"]
    genos, ref = genotype_summary(seqs, reference="AAA")
    assert ref == "AAA"
    assert genos == ["", "2T"]


def test_supplied_reference_and_computed_can_differ():
    # Demonstrate the two paths give different reference (the reason to support it).
    seqs = ["TT", "TT", "TA"]
    _, computed = genotype_summary(seqs)            # majority -> "TT"
    _, supplied = genotype_summary(seqs, reference="AA")
    assert computed == "TT"
    assert supplied == "AA"


# --------------------------------------------------------------------------- #
# select_sequences -- the filters the golden fixture doesn't exercise
# --------------------------------------------------------------------------- #
def test_select_drops_low_coverage():
    # depth_min default 10; coverage 5 should be dropped, 20 kept.
    recs = [_rec("low", "AAAAAA", 5), _rec("high", "AAAAAA", 20)]
    kept = select_sequences(recs, start=0, end=6)
    assert [r.cbc_id for r in kept] == ["high"]


def test_select_drops_too_many_gaps():
    # max_gaps default 3; a seq with 3 gaps in region is dropped (filter is < 3).
    recs = [
        _rec("gappy", "A--A--", 50),   # 4 gaps in [0,6) -> dropped
        _rec("clean", "AACAAA", 50),   # 0 gaps -> kept
    ]
    kept = select_sequences(recs, start=0, end=6)
    assert [r.cbc_id for r in kept] == ["clean"]


def test_select_trims_to_region():
    recs = [_rec("x", "AACCGGTT", 50)]
    kept = select_sequences(recs, start=2, end=6)
    assert kept[0].seq == "CCGG"


def test_select_gap_count_only_within_region():
    # Gaps OUTSIDE the region don't count toward the filter.
    recs = [_rec("edge", "---AAAA---", 50)]   # gaps at edges, none in [3,7)
    kept = select_sequences(recs, start=3, end=7)
    assert [r.cbc_id for r in kept] == ["edge"]
    assert kept[0].seq == "AAAA"


# --------------------------------------------------------------------------- #
# parse_consensus_fasta
# --------------------------------------------------------------------------- #
def test_parse_fasta_reads_records_and_metadata(tmp_path):
    fa = tmp_path / "x_allConsensus.fasta"
    fa.write_text(
        ">CELL1 ref coverage:37 length:6\nAACCGG\n"
        ">CELL2 ref coverage:12 length:6\nAATCGG\n"
    )
    recs = parse_consensus_fasta(fa)
    assert [r.cbc_id for r in recs] == ["CELL1", "CELL2"]
    assert recs[0].seq == "AACCGG"
    assert recs[0].coverage == 37.0
    assert recs[1].coverage == 12.0


# --------------------------------------------------------------------------- #
# GOLDEN: reference=None path reproduces the frozen, hand-verified output
# --------------------------------------------------------------------------- #
def test_consensus_matches_golden(data_dir, golden_dir, tmp_run_dir):
    """run() with reference=None reproduces the original ConsensusTool output."""
    fixture = data_dir / "consensus_test_allConsensus.fasta"
    ref_golden = golden_dir / "consensus_reference.txt"
    csv_golden = golden_dir / "consensus_genotypes.csv"
    for p in (fixture, ref_golden, csv_golden):
        if not p.exists():
            pytest.skip(f"missing {p.name}; run the fixture/golden generators first.")

    # Run into a temp prefix so we don't clobber anything.
    out_prefix = str(tmp_run_dir / "consensus_test")
    result = run(fasta=str(fixture), start=3, end=27, reference=None,
                 out_prefix=out_prefix)

    # 1. Reference matches the frozen consensus.
    expected_ref = ref_golden.read_text().strip()
    assert result["reference"] == expected_ref

    # 2. Per-CBC genotypes match the frozen table.
    expected = {}
    with open(csv_golden) as fh:
        for row in csv.DictReader(fh):
            expected[row["CBC_ID"]] = row["genotype"] or ""

    got = {r["CBC_ID"]: r["genotype"] for r in result["records"]}
    assert got == expected
