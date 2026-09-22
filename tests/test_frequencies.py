"""
Tests for frequencies.py -- allele frequencies with per-position denominators.

THE CENTRAL PROPERTY, and the reason this stage exists separately from the
genotype network: a cell that covers part of the reference must count at the
positions it covers and be absent elsewhere. No cell should be discarded for
partial coverage, and no cell should sit in a denominator at a position it has
no data for.
"""

from __future__ import annotations

import pytest

from anchovy.consensus import ConsensusRecord
from anchovy.frequencies import (
    cell_allele_rows,
    counts_path_for,
    population_rows,
    read_counts,
    run,
)
from anchovy.schema import AlleleFrequencyColumns as F


def _rec(cbc_id, seq):
    return ConsensusRecord(cbc_id=cbc_id, seq=seq, description=cbc_id,
                           coverage=50, length=len(seq))


# --------------------------------------------------------------------------- #
# population_rows -- per-position denominators
# --------------------------------------------------------------------------- #
def test_gapped_cells_leave_the_denominator_where_they_have_no_data():
    """A partial cell counts where it covers and vanishes where it does not."""
    ref = "AAAA"
    recs = [
        _rec("full_ref", "AAAA"),
        _rec("full_alt", "AATA"),
        _rec("partial", "AA--"),   # no data at positions 3 and 4
    ]
    rows = {r[F.POSITION]: r for r in population_rows(recs, ref)}
    # Position 3: only two cells could be counted, so the denominator is 2 and
    # the one alt call is 1/2 -- not 1/3, which would silently dilute it.
    assert rows[3][F.CELLS_CALLED] == 2
    assert rows[3][F.CELLS_ALT] == 1
    assert rows[3][F.FREQ_CELLS] == 0.5


def test_ambiguous_calls_are_in_neither_numerator_nor_denominator():
    # An IUPAC code means the reads disagreed, so the cell has no confident
    # call. Counting it as reference would understate the frequency.
    ref = "AA"
    recs = [_rec("alt", "AT"), _rec("ref", "AA"), _rec("mixed", "AY")]
    rows = {r[F.POSITION]: r for r in population_rows(recs, ref)}
    assert rows[2][F.CELLS_CALLED] == 2
    assert rows[2][F.FREQ_CELLS] == 0.5


def test_subset_columns_report_the_network_cells_separately():
    """The subset is reported, not assumed representative.

    Here the strict subset is enriched for the variant relative to the whole
    population -- exactly the bias these columns exist to make visible.
    """
    ref = "AA"
    recs = [_rec("a", "AT"), _rec("b", "AT"), _rec("c", "AA"), _rec("d", "AA")]
    rows = population_rows(recs, ref, subset_ids={"a", "c"})
    row = rows[0]
    assert (row[F.CELLS_ALT], row[F.CELLS_CALLED], row[F.FREQ_CELLS]) == (2, 4, 0.5)
    assert row[F.CELLS_ALT_SUBSET] == 1
    assert row[F.CELLS_CALLED_SUBSET] == 2


def test_no_row_for_a_position_nobody_calls():
    recs = [_rec("a", "A-"), _rec("b", "A-")]
    assert [r[F.POSITION] for r in population_rows(recs, "AA")] == []


def test_min_cells_filters_singletons():
    ref = "AA"
    recs = [_rec("a", "AT"), _rec("b", "AA"), _rec("c", "AA")]
    assert population_rows(recs, ref, min_cells=2) == []
    assert len(population_rows(recs, ref, min_cells=1)) == 1


# --------------------------------------------------------------------------- #
# cell vote vs read weighting -- the reason both are reported
# --------------------------------------------------------------------------- #
def test_cell_votes_and_read_sums_can_disagree():
    """One deep cell can carry the read frequency while losing the cell vote.

    This is why the population table reports both: reads within a cell are
    amplification copies, so summing them weights cells by depth. Here a single
    variant cell sequenced deeply outweighs two reference cells by reads while
    being outvoted 2:1 by cells.
    """
    ref = "A"
    recs = [_rec("deep_alt", "T"), _rec("ref1", "A"), _rec("ref2", "A")]
    read_totals = {1: {"A": 20, "C": 0, "G": 0, "T": 80, "N": 0, "gap": 0}}
    row = population_rows(recs, ref, read_totals=read_totals)[0]
    assert row[F.FREQ_CELLS] == pytest.approx(1 / 3)
    assert row[F.FREQ_READS] == pytest.approx(0.8)


def test_read_columns_are_blank_without_counts():
    # The stage stays runnable against results produced before --counts existed.
    row = population_rows([_rec("a", "T"), _rec("b", "A")], "A")[0]
    assert row[F.READS_ALT] == "" and row[F.FREQ_READS] == ""


# --------------------------------------------------------------------------- #
# per-cell pileup rows
# --------------------------------------------------------------------------- #
def _counts(**cols):
    base = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "gap": 0}
    base.update(cols)
    return base


def test_cell_rows_report_within_cell_mixtures_the_consensus_hides():
    # 2 of 8 reads carry the variant: the cell's own consensus calls reference,
    # so this mixture is invisible downstream of sam2consensus.
    rows = cell_allele_rows("cell", {1: _counts(A=6, T=2)}, "A")
    assert len(rows) == 1
    assert rows[0]["reads"] == 2 and rows[0]["depth"] == 8
    assert rows[0]["freq"] == 0.25


def test_singleton_alt_reads_are_dropped_by_default():
    # Every covered position carries error reads; min_alt_reads=2 is what keeps
    # the table to variants rather than the error spectrum.
    assert cell_allele_rows("cell", {1: _counts(A=9, T=1)}, "A") == []
    assert len(cell_allele_rows("cell", {1: _counts(A=9, T=1)}, "A",
                                min_alt_reads=1)) == 1


def test_depth_counts_deletions_and_ns():
    # A read carrying a deletion is still a read at that position, so it belongs
    # in the denominator -- matching how sam2consensus computes coverage.
    rows = cell_allele_rows("cell", {1: _counts(A=5, T=3, N=1, gap=1)}, "A")
    assert rows[0]["depth"] == 10
    assert rows[0]["freq"] == 0.3


def test_reference_allele_never_gets_a_row():
    assert cell_allele_rows("cell", {1: _counts(A=10)}, "A") == []


# --------------------------------------------------------------------------- #
# IO
# --------------------------------------------------------------------------- #
def test_read_counts_round_trip(tmp_path):
    p = tmp_path / "c.tsv"
    p.write_text("position\tA\tC\tG\tT\tN\tgap\n7\t3\t0\t1\t0\t0\t2\n")
    assert read_counts(p) == {7: {"A": 3, "C": 0, "G": 1, "T": 0,
                                  "N": 0, "gap": 2}}


def test_counts_path_is_derived_from_the_consensus_id(tmp_path):
    # sam2consensus writes "{ref}__{prefix}_counts.tsv" and headers as
    # "{prefix}|c{threshold}", so the id locates the file.
    (tmp_path / "ref__cellA_BC_counts.tsv").write_text("position\tA\tC\tG\tT\tN\tgap\n")
    assert counts_path_for("cellA_BC|c50", tmp_path, "ref") is not None
    assert counts_path_for("missing|c50", tmp_path, "ref") is None


def test_run_warns_on_a_mismatched_reference(tmp_path):
    """A shifted reference reports everything as mutant instead of failing."""
    fasta = tmp_path / "x_allConsensus.fasta"
    fasta.write_text(">a ref coverage:50 length:4\nACGT\n")
    with pytest.warns(UserWarning, match="Positions are compared by index"):
        run(str(fasta), reference="ACG", out_prefix=str(tmp_path / "o"))
