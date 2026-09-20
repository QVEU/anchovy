"""
test_reference_io.py -- reading a reference sequence, and refusing a wrong one.

THE BUG THIS GUARDS
-------------------
Both the consensus and annotate stages used to read their reference with a bare
`read().strip()`. Handed a FASTA -- the format every reference genome actually
ships in -- that spliced the '>' header line into the sequence. Nothing raised:
every genome coordinate shifted by the length of the header, so the run produced
a full set of confidently wrong calls. A reference that is merely too short at
least failed loudly; this failed silently, which is worse.

Two guards, doing different jobs:
  * io.read_reference_sequence -- parse FASTA properly at the boundary, so a bad
    sequence is never constructed in the first place.
  * genotype_summary's length check -- catch a reference that does not line up
    with the alignment, whatever produced it.
"""

from __future__ import annotations

import csv
import json

import pandas as pd
import pytest

from anchovy import annotate
from anchovy.cli import main
from anchovy.consensus import genotype_summary
from anchovy.io import read_reference_sequence


# --------------------------------------------------------------------------- #
# read_reference_sequence
# --------------------------------------------------------------------------- #
def test_reads_fasta_without_the_header(tmp_path):
    """THE REGRESSION: the header must not end up in the sequence."""
    p = tmp_path / "ref.fasta"
    p.write_text(">testref some description here\nACGTACGT\n")
    seq = read_reference_sequence(p)
    assert seq == "ACGTACGT"
    assert ">" not in seq and "\n" not in seq


def test_reads_wrapped_fasta_as_one_sequence(tmp_path):
    """Line-wrapped FASTA (how real reference genomes are distributed)."""
    p = tmp_path / "ref.fasta"
    p.write_text(">ref\nACGT\nACGT\nGGGG\n")
    assert read_reference_sequence(p) == "ACGTACGTGGGG"


def test_reads_raw_sequence_file(tmp_path):
    """The existing format must keep working -- it is what the pipeline writes."""
    p = tmp_path / "ref.txt"
    p.write_text("ACGTACGT\n")
    assert read_reference_sequence(p) == "ACGTACGT"


def test_raw_and_fasta_of_the_same_sequence_agree(tmp_path):
    (tmp_path / "a.txt").write_text("ACGTACGTGGGG\n")
    (tmp_path / "b.fasta").write_text(">ref\nACGTACGT\nGGGG\n")
    assert (read_reference_sequence(tmp_path / "a.txt")
            == read_reference_sequence(tmp_path / "b.fasta"))


def test_case_is_preserved(tmp_path):
    """Callers normalize; the reader does not decide for them."""
    p = tmp_path / "ref.fasta"
    p.write_text(">ref\nacgtACGT\n")
    assert read_reference_sequence(p) == "acgtACGT"


def test_multi_record_fasta_is_refused(tmp_path):
    """A segmented genome must fail loudly, not silently use the first record."""
    p = tmp_path / "segments.fasta"
    p.write_text(">seg1\nACGT\n>seg2\nTTTT\n")
    with pytest.raises(ValueError, match="found 2 FASTA records"):
        read_reference_sequence(p)


def test_empty_file_is_refused(tmp_path):
    p = tmp_path / "empty.txt"
    p.write_text("   \n")
    with pytest.raises(ValueError, match="no reference sequence"):
        read_reference_sequence(p)


def test_header_only_fasta_is_refused(tmp_path):
    p = tmp_path / "headeronly.fasta"
    p.write_text(">ref\n")
    with pytest.raises(ValueError, match="no reference sequence"):
        read_reference_sequence(p)


def test_real_template_fixture_reads_at_full_length(data_dir):
    """The actual fixture FASTA, read straight through with no hand-stripping."""
    seq = read_reference_sequence(data_dir / "mapping" / "template.fasta")
    assert len(seq) == 600
    assert set(seq) <= set("ACGT")


# --------------------------------------------------------------------------- #
# The length check in genotype_summary
# --------------------------------------------------------------------------- #
def test_mismatched_reference_length_is_refused():
    """Positions are compared by index, so a wrong-length reference is fatal.

    This is the shape the old FASTA bug produced: a reference longer than the
    alignment by the width of its header line.
    """
    with pytest.raises(ValueError, match="does not match the alignment length"):
        genotype_summary(["AAAT", "AAAA"], reference=">testref AAAT")

    with pytest.raises(ValueError, match="does not match the alignment length"):
        genotype_summary(["AAAT", "AAAA"], reference="AAA")      # too short


def test_matching_reference_length_is_accepted():
    genos, ref = genotype_summary(["AAAT", "AAAA"], reference="AAAA")
    assert ref == "AAAA" and genos == ["4T", ""]


def test_computed_reference_is_unaffected_by_the_length_check():
    """No reference supplied -> nothing to validate, original path untouched."""
    genos, ref = genotype_summary(["AAAT", "AAAA"])
    assert ref == "AAAA" and genos == ["4T", ""]


# --------------------------------------------------------------------------- #
# End to end through the CLI boundary
# --------------------------------------------------------------------------- #
def test_consensus_cli_accepts_a_fasta_reference(data_dir, tmp_run_dir):
    """`anchovy consensus --reference template.fasta` must just work.

    Before the fix this ran without error and produced wrong calls, because the
    header shifted every coordinate. Now the reference is parsed properly, so the
    cell carrying the planted variant is the one reported as mutant.
    """
    expected = json.loads(
        (data_dir / "mapping" / "mapping_expected.json").read_text())
    all_consensus = tmp_run_dir / "test_allConsensus.fasta"
    all_consensus.write_text("".join(
        f">{cell}_BC|c50 ref coverage:5.0 length:{len(seq)}\n{seq}\n"
        for cell, seq in expected["consensus"].items()))

    out_prefix = str(tmp_run_dir / "out")
    rc = main([
        "consensus", str(all_consensus), "100", "500", "--whole-reference",
        "--depth-min", "1",
        "--reference", str(data_dir / "mapping" / "template.fasta"),
        "--out-prefix", out_prefix,
    ])
    assert rc == 0

    with open(f"{out_prefix}_filtConsensus.csv") as fh:
        rows = {r["CBC_ID"]: r["genotype"] for r in csv.DictReader(fh)}

    # cellA matches the template, cellB carries the planted variant -- so with the
    # template as reference it is cellB that is called, at the genome position.
    cell_a = next(k for k in rows if k.startswith("cellA"))
    cell_b = next(k for k in rows if k.startswith("cellB"))
    assert rows[cell_a] == ""
    assert rows[cell_b] == f"{expected['variant_genome_pos']}{expected['variant_alt_base']}"


def test_annotate_accepts_a_fasta_reference(data_dir, tmp_path):
    """The annotate stage has the same hazard and the same fix."""
    expected = json.loads(
        (data_dir / "mapping" / "mapping_expected.json").read_text())
    pos = expected["variant_genome_pos"]
    cons = tmp_path / "filtConsensus.csv"
    pd.DataFrame([{"CBC_ID": "cellB", "genotype": f"{pos}{expected['variant_alt_base']}",
                   "sequence": "", "description": "coverage:50"}]).to_csv(cons, index=False)

    result = annotate.run(
        filt_consensus_csv=str(cons),
        reference_file=str(data_dir / "mapping" / "template.fasta"),
        out_prefix=str(tmp_path / "out"), network=False,
        gff=str(data_dir / "mapping" / "regions.gff3"))

    row = result["regions"].iloc[0]
    # The WT base is read from the real template, not from a header-shifted one.
    assert row["wt_base"] == expected["variant_ref_base"]
    assert row["residue"] == expected["variant_residue"]
    assert row["region"] == expected["variant_region"]
