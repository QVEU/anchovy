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
    # computed-consensus case. Reference says position 1 should be 'A'.
    seqs = ["GACG", "GTCG", "GACG"]
    genos, ref = genotype_summary(seqs, reference="AACG")
    assert ref == "AACG"
    # pos1: every cell has 'G' where the reference says 'A'. That is a FIXED
    # DIFFERENCE and it IS reported, for all three cells -- see the test below
    # for why. pos2: A/T/A, so only the middle cell differs from ref 'A'.
    assert genos == ["1G", "1G_2T", "1G"]


def test_fixed_difference_from_supplied_reference_is_reported():
    """A mutation shared by EVERY cell must still be called against the reference.

    This is a deliberate behavior change from the original ConsensusTool, which
    chose variant sites by asking whether a column varied ACROSS THE CELLS. Under
    that rule a position where the whole population has moved off the reference
    never became a variant site, so every cell silently read as wild-type there
    -- the more uniformly a mutation had swept the population, the more certainly
    it was dropped. On a passaged or lab-adapted stock that can discard most of
    the real variants, so the reference is now what sites are judged against.
    """
    # Every cell carries 'G' at position 3 where the reference says 'T'.
    seqs = ["AAG", "AAG", "AAG"]
    genos, _ = genotype_summary(seqs, reference="AAT")
    assert genos == ["3G", "3G", "3G"]

    # The computed-reference path is unaffected: with no reference supplied the
    # consensus IS 'G' there, so there is nothing to report.
    genos_computed, ref = genotype_summary(seqs)
    assert ref == "AAG"
    assert genos_computed == ["", "", ""]


def test_widened_predicate_is_a_no_op_for_the_computed_reference():
    """The safety property behind the change above.

    Choosing sites by "any cell differs from the reference" instead of "the column
    varies" is the SAME SET whenever the reference is the computed consensus,
    because the consensus is always one of the characters present. That is what
    lets the supplied-reference behavior change without disturbing any frozen
    golden. Checked here on alignments that mix agreement, variation and gaps.
    """
    for seqs in (["ACGT", "ACGT", "ACGT"],          # no variation at all
                 ["ACGT", "AGGT", "ACGA"],          # several variant columns
                 ["A-GT", "ACGT", "A--T"],          # gap-bearing columns
                 ["AAAA"]):                          # single sequence
        genos, ref = genotype_summary(seqs)
        # Under either rule a cell can only be called where it differs from the
        # consensus, so every token must correspond to a real difference.
        for seq, geno in zip(seqs, genos):
            positions = [int(t[:-1]) for t in geno.split("_") if t]
            assert all(seq[i - 1] != ref[i - 1] for i in positions)
            # ...and nothing that differs from the consensus was missed.
            expected = {i + 1 for i in range(len(seq)) if seq[i] != ref[i]}
            assert set(positions) == expected


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


# --------------------------------------------------------------------------- #
# select_sequences -- breadth / depth-where-called (the unconflated filters)
# --------------------------------------------------------------------------- #
# The case that motivated them: `coverage` is (summed depth over covered
# positions) / (full length), so a deep-but-partial cell scores LOWER than a
# shallow-but-complete one. Each test below pins one half of that apart.
def test_breadth_and_depth_called_are_off_by_default():
    # Half-covered but deep. Nothing is set, so only depth_min applies and the
    # record survives exactly as it did before these filters existed.
    recs = [_rec("partial", "AAAAA-----", 20)]
    assert [r.cbc_id for r in select_sequences(recs, config=ConsensusConfig(
        depth_min=10))] == ["partial"]


def test_min_breadth_drops_narrow_cells():
    cfg = ConsensusConfig(depth_min=0, min_breadth=0.5)
    recs = [
        _rec("narrow", "AAA-------", 20),   # 30% called -> dropped
        _rec("wide", "AAAAAA----", 20),     # 60% called -> kept
    ]
    assert [r.cbc_id for r in select_sequences(recs, config=cfg)] == ["wide"]


def test_min_depth_called_undoes_the_breadth_penalty():
    """The whole point: judge a cell on depth where it called, not on average.

    Both cells have the same genome-wide `coverage` of 16, so depth_min alone
    cannot tell them apart. Factoring breadth out separates them cleanly:
      deep    16 * 10/4 = 40x over the 4 positions it called
      shallow 16 * 10/10 = 16x across all 10
    """
    cfg = ConsensusConfig(depth_min=0, min_depth_called=20)
    recs = [
        _rec("deep", "AAAA------", 16),
        _rec("shallow", "AAAAAAAAAA", 16),
    ]
    assert [r.cbc_id for r in select_sequences(recs, config=cfg)] == ["deep"]


def test_depth_min_would_have_kept_exactly_the_wrong_one():
    # The inverse of the test above, proving the old filter is not merely
    # coarser but actively inverted on this pair: it keeps the 16x-everywhere
    # cell and discards the 40x-where-called one.
    recs = [
        _rec("deep", "AAAA------", 16),      # 40x where called
        _rec("shallow", "AAAAAAAAAA", 21),   # 21x everywhere
    ]
    kept = select_sequences(recs, config=ConsensusConfig(depth_min=20))
    assert [r.cbc_id for r in kept] == ["shallow"]


def test_filters_compose_and_report_their_own_tolls():
    cfg = ConsensusConfig(depth_min=5, min_breadth=0.5, min_depth_called=20)
    recs = [
        _rec("low", "AAAAAAAAAA", 4),        # under depth_min
        _rec("narrow", "AAA-------", 30),    # 30% called
        _rec("shallow", "AAAAAAAAAA", 10),   # 10x where called
        _rec("good", "AAAAAAAA--", 20),      # 80% called, 25x where called
    ]
    stats: dict = {}
    kept = select_sequences(recs, config=cfg, stats=stats)
    assert [r.cbc_id for r in kept] == ["good"]
    assert stats["dropped_low_depth"] == 1
    assert stats["dropped_narrow"] == 1
    assert stats["dropped_shallow_called"] == 1


def test_fully_gapped_record_does_not_divide_by_zero():
    # An all-gap consensus reaches the filter when depth_min is lowered, and
    # n_called is then 0. It must be dropped, not raise.
    cfg = ConsensusConfig(depth_min=0, min_breadth=0.1, min_depth_called=1)
    assert select_sequences([_rec("empty", "----------", 1)], config=cfg) == []


def test_min_breadth_is_measured_over_the_window_when_given():
    """Ragged flanks must not count against a cell that covers the core.

    The motivating case: amplicon reads never reach the extreme ends, so every
    real cell carries gap flanks. Judged whole-genome this cell is 50% covered
    and fails; judged over the core it named, it is complete and passes.
    """
    cfg = ConsensusConfig(depth_min=0, min_breadth=1.0)
    rec = _rec("flanked", "-----AAAAAAAAAA-----", 20)
    assert select_sequences([rec], config=cfg) == []
    kept = select_sequences([rec], start=5, end=15, config=cfg, trim=False)
    assert [r.cbc_id for r in kept] == ["flanked"]


def test_window_breadth_still_rejects_a_hole_in_the_core():
    # The window must not become a rubber stamp: a gap INSIDE it still fails,
    # which is the whole point of asking for full coverage of the core.
    cfg = ConsensusConfig(depth_min=0, min_breadth=1.0)
    rec = _rec("holed", "-----AAAA-AAAAA-----", 20)
    assert select_sequences([rec], start=5, end=15, config=cfg, trim=False) == []


def test_depth_called_stays_whole_sequence_under_a_window():
    # depth_called cannot be windowed (the per-position depths are gone), so a
    # window must not silently change its denominator. Coverage 10 over a
    # 20-col record with 10 called positions is 20x called, window or not.
    cfg = ConsensusConfig(depth_min=0, min_depth_called=20)
    rec = _rec("half", "-----AAAAAAAAAA-----", 10)
    assert [r.cbc_id for r in select_sequences(
        [rec], start=5, end=15, config=cfg, trim=False)] == ["half"]


# --------------------------------------------------------------------------- #
# run() -- consensus sequences that are not reference length
# --------------------------------------------------------------------------- #
# sam2consensus appends called insertions as extra columns, so a cell can come
# out LONGER than the reference. Genotyping compares by index, so such a cell is
# shifted from the insertion onward. It used to be caught only by a check
# against sequences[0] -- i.e. by luck -- and then vanish into annotate's
# hypermutant cap.
def _write_fasta(tmp_path, entries):
    p = tmp_path / "x_allConsensus.fasta"
    p.write_text("".join(
        f">{name} ref coverage:{cov} length:{len(seq)}\n{seq}\n"
        for name, seq, cov in entries))
    return p


def test_run_drops_insertion_shifted_cells(tmp_path, capsys):
    ref = "ACGTACGTAC"
    fasta = _write_fasta(tmp_path, [
        ("normal", "ACGTACGTAC", 50),
        ("inserted", "ACGTTACGTAC", 50),      # 11 nt: one inserted column
        ("normal2", "ACGTACGTAT", 50),
    ])
    result = run(str(fasta), reference=ref, trim=False,
                 out_prefix=str(tmp_path / "out"),
                 config=ConsensusConfig(depth_min=10))
    assert [r["CBC_ID"] for r in result["records"]] == ["normal", "normal2"]
    assert result["stats"]["dropped_length_mismatch"] == 1
    assert "not 10 nt" in capsys.readouterr().out


def test_run_survives_a_shifted_cell_sorting_first(tmp_path):
    """The old check looked at sequences[0], so this ordering used to raise."""
    ref = "ACGTACGTAC"
    fasta = _write_fasta(tmp_path, [
        ("inserted", "ACGTTACGTAC", 50),      # longer, and FIRST
        ("normal", "ACGTACGTAC", 50),
    ])
    result = run(str(fasta), reference=ref, trim=False,
                 out_prefix=str(tmp_path / "out"),
                 config=ConsensusConfig(depth_min=10))
    assert [r["CBC_ID"] for r in result["records"]] == ["normal"]


def test_run_uses_modal_length_when_no_reference_supplied(tmp_path):
    # With a computed reference there is no external length to trust, so the
    # majority length is the only defensible answer.
    fasta = _write_fasta(tmp_path, [
        ("a", "ACGTACGTAC", 50),
        ("b", "ACGTACGTAT", 50),
        ("odd", "ACGTTACGTAC", 50),
    ])
    result = run(str(fasta), trim=False, out_prefix=str(tmp_path / "out"),
                 config=ConsensusConfig(depth_min=10))
    assert [r["CBC_ID"] for r in result["records"]] == ["a", "b"]
    assert result["stats"]["dropped_length_mismatch"] == 1
