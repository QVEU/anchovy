"""
test_fasta.py -- unit + golden tests for the fasta stage.

The split (build_cell_fastas pure, write_cell_fastas thin) pays off here: the
core logic is tested with tiny in-memory DataFrames -- no files -- and the golden
test asserts on the returned dict, chaining directly off the extract stage's
frozen CSV. That chaining is the first place two stages connect through the
shared schema, which is the whole point of schema.py.
"""

from __future__ import annotations

import pandas as pd
import pytest

from anchovy.config import FastaConfig
from anchovy.fasta import build_cell_fastas, write_cell_fastas
from anchovy.schema import AnchovyColumns


def _mini_df(rows):
    """Build a minimal anchovy-shaped frame from (cbc, read, seq) tuples."""
    return pd.DataFrame(
        rows, columns=[AnchovyColumns.CBC, AnchovyColumns.READ, AnchovyColumns.MAPPED_SEQ]
    )


# --------------------------------------------------------------------------- #
# build_cell_fastas -- the pure filter/format logic
# --------------------------------------------------------------------------- #
def test_barcode_below_threshold_is_dropped():
    # 2 reads for one barcode, threshold 5 -> nothing emitted.
    df = _mini_df([("BC1", "r1", "AAA"), ("BC1", "r2", "CCC")])
    result = build_cell_fastas(df, FastaConfig(min_reads_per_cbc=5))
    assert result == {}


def test_barcode_at_threshold_is_kept():
    # Exactly threshold reads -> kept (filter is >=, matching the original).
    rows = [("BC1", f"r{i}", "AAA") for i in range(5)]
    result = build_cell_fastas(_mini_df(rows), FastaConfig(min_reads_per_cbc=5))
    assert set(result.keys()) == {"BC1"}


def test_fasta_record_format():
    # Records are ">{read}\n{seq}", joined by newline, in row order.
    rows = [("BC1", "r1", "AAA"), ("BC1", "r2", "CCC"), ("BC1", "r3", "GGG")]
    result = build_cell_fastas(_mini_df(rows), FastaConfig(min_reads_per_cbc=3))
    assert result["BC1"] == ">r1\nAAA\n>r2\nCCC\n>r3\nGGG"


def test_multiple_barcodes_filtered_independently():
    # BC1 has 3 reads (kept at threshold 3), BC2 has 1 (dropped).
    rows = [("BC1", "a", "AA"), ("BC1", "b", "AA"), ("BC1", "c", "AA"),
            ("BC2", "d", "CC")]
    result = build_cell_fastas(_mini_df(rows), FastaConfig(min_reads_per_cbc=3))
    assert set(result.keys()) == {"BC1"}


def test_barcode_whitespace_is_stripped():
    # Barcodes are stripped in the output keys (matches j.strip() in original).
    rows = [(" BC1 ", f"r{i}", "AA") for i in range(3)]
    result = build_cell_fastas(_mini_df(rows), FastaConfig(min_reads_per_cbc=3))
    assert set(result.keys()) == {"BC1"}


def test_default_threshold_is_five():
    # FastaConfig() default reproduces the original literal 5.
    assert FastaConfig().min_reads_per_cbc == 5


# --------------------------------------------------------------------------- #
# write_cell_fastas -- the thin I/O layer
# --------------------------------------------------------------------------- #
def test_write_creates_one_file_per_barcode(tmp_path):
    fastas = {"BC1": ">r1\nAAA", "BC2": ">r2\nCCC"}
    written = write_cell_fastas(fastas, tmp_path)
    assert [p.name for p in written] == ["BC1.fa", "BC2.fa"]
    assert (tmp_path / "BC1.fa").read_text() == ">r1\nAAA"


# --------------------------------------------------------------------------- #
# Golden test -- chains off the frozen extract output
# --------------------------------------------------------------------------- #
def test_fasta_matches_golden(golden_dir):
    """build_cell_fastas on the frozen extract CSV produces expected cells.

    Input is the extract stage's golden CSV (the stages chaining through the
    shared schema). We assert on the per-cell FASTA dict. On first run, if the
    fasta golden doesn't exist yet, we generate it from current behavior and
    skip -- the same freeze-then-guard pattern as the extract golden.
    """
    extract_csv = golden_dir / "extract_anchovy.csv"
    if not extract_csv.exists():
        pytest.skip("extract golden missing; generate it first.")

    df = pd.read_csv(extract_csv)
    result = build_cell_fastas(df)  # default config (min_reads_per_cbc=5)

    golden_path = golden_dir / "fasta_cells.json"

    import json
    if not golden_path.exists():
        # Freeze current behavior on first run, then skip with a clear message.
        golden_path.write_text(json.dumps(result, indent=2, sort_keys=True))
        pytest.skip(
            f"Wrote initial fasta golden to {golden_path}. "
            "Inspect it, commit it, then re-run to enable the guard."
        )

    expected = json.loads(golden_path.read_text())
    assert result == expected
