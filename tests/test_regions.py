"""
test_regions.py -- unit tests for the region model, GFF3 parser, and both-strand
codon annotation.

The codon math here was verified numerically (against real Dengue coordinates and
a synthetic minus-strand gene translated independently) before the code was
written; these tests lock in those verified answers.
"""

from __future__ import annotations

import pytest

from anchovy.regions import (
    Region, parse_gff3, validate_regions, classify,
    annotate_in_region, annotate_mutation,
)


# --------------------------------------------------------------------------- #
# GFF3 parsing
# --------------------------------------------------------------------------- #
GFF3 = """\
##gff-version 3
# a comment line
ref\tanchovy\tfive_prime_UTR\t1\t96\t.\t+\t.\tID=5UTR;Name=5UTR
ref\tanchovy\tCDS\t97\t10272\t.\t+\t0\tID=poly;Name=polyprotein
ref\tanchovy\tmature_protein_region\t7570\t10269\t.\t+\t0\tID=ns5;Name=NS5
ref\tanchovy\tthree_prime_UTR\t10273\t10727\t.\t+\t.\tID=3UTR;Name=3UTR
gene\tanchovy\tgene\t1\t10727\t.\t+\t.\tID=ignored_gene
"""


def test_parse_gff3_reads_expected_regions(tmp_path):
    p = tmp_path / "x.gff3"
    p.write_text(GFF3)
    regions = parse_gff3(str(p))
    names = [r.name for r in regions]
    # gene row ignored (unrecognized type); the other four kept.
    assert names == ["5UTR", "polyprotein", "NS5", "3UTR"]
    poly = next(r for r in regions if r.name == "polyprotein")
    assert poly.coding and poly.start == 97 and poly.end == 10272 and poly.phase == 0
    utr = next(r for r in regions if r.name == "5UTR")
    assert not utr.coding


def test_parse_gff3_rejects_bad_strand(tmp_path):
    p = tmp_path / "bad.gff3"
    p.write_text("ref\ta\tCDS\t1\t9\t.\tX\t0\tID=x\n")
    with pytest.raises(ValueError, match="strand"):
        parse_gff3(str(p))


def test_parse_gff3_rejects_start_after_end(tmp_path):
    p = tmp_path / "bad2.gff3"
    p.write_text("ref\ta\tCDS\t99\t9\t.\t+\t0\tID=x\n")
    with pytest.raises(ValueError, match="start"):
        parse_gff3(str(p))


# --------------------------------------------------------------------------- #
# validity checks
# --------------------------------------------------------------------------- #
def test_validate_flags_non_multiple_of_three():
    # A coding region of length 10 (not divisible by 3) should warn.
    regions = [Region("bad", 1, 10, "+", 0, coding=True)]
    warnings = validate_regions(regions)
    assert any("multiple of 3" in w for w in warnings)


def test_validate_passes_clean_cds():
    # Dengue polyprotein: (10272-97+1) = 10176, divisible by 3 -> no warning.
    regions = [Region("polyprotein", 97, 10272, "+", 0, coding=True)]
    assert validate_regions(regions) == []


# --------------------------------------------------------------------------- #
# classification / overlap (long-format)
# --------------------------------------------------------------------------- #
def test_classify_returns_all_containing_regions():
    regions = [
        Region("polyprotein", 97, 10272, "+", 0, coding=True),
        Region("NS5", 7570, 10269, "+", 0, coding=True),
        Region("5UTR", 1, 96, "+", 0, coding=False),
    ]
    # genome 7573 is inside both polyprotein and NS5
    hits = classify(regions, 7573)
    assert {r.name for r in hits} == {"polyprotein", "NS5"}
    # genome 50 only in 5UTR
    assert [r.name for r in classify(regions, 50)] == ["5UTR"]


def test_annotate_mutation_long_format_two_rows():
    # Overlapping regions -> two rows for one mutation, each in its own frame.
    # Use a short reference so codon extraction works; only frames matter here.
    ref = "A" * 10300  # dummy; residue numbers are what we check
    regions = [
        Region("polyprotein", 97, 10272, "+", 0, coding=True),
        Region("NS5", 7570, 10269, "+", 0, coding=True),
    ]
    rows = annotate_mutation(regions, 7573, "A", "G", ref)
    assert len(rows) == 2
    by_region = {r["region"]: r for r in rows}
    # verified earlier: polyprotein residue 2493, NS5 residue 2
    assert by_region["polyprotein"]["residue"] == 2493
    assert by_region["NS5"]["residue"] == 2
    # mutation_id is genome/forward-strand in both
    assert by_region["polyprotein"]["mutation_id"] == "polyprotein:A7573G"
    assert by_region["NS5"]["mutation_id"] == "NS5:A7573G"


def test_annotate_intergenic_when_no_region():
    rows = annotate_mutation([], 500, "C", "T", "N" * 1000)
    assert len(rows) == 1
    assert rows[0]["region"] == "intergenic"
    assert rows[0]["mutation_id"] == "intergenic:C500T"


# --------------------------------------------------------------------------- #
# plus-strand coding annotation (real amino-acid calls)
# --------------------------------------------------------------------------- #
def test_plus_strand_amino_acid_call():
    # Reference where CDS starts at position 1, codon 1 = ATG (M), codon 2 = AAA (K).
    ref = "ATGAAAGATTTT"   # M K D F
    region = Region("cds", 1, 12, "+", 0, coding=True)
    # mutate position 6 (3rd base of codon 2, AAA): A->G => AAG still K (synonymous)
    row = annotate_in_region(region, 6, "A", "G", ref)
    assert row["wt_aa"] == "K" and row["mut_aa"] == "K"
    assert row["sub_class"] == "Syn" and row["residue"] == 2
    # mutate position 5 (2nd base of codon 2): A->T => ATA = I (non-syn)
    row2 = annotate_in_region(region, 5, "A", "T", ref)
    assert row2["wt_aa"] == "K" and row2["mut_aa"] == "I"
    assert row2["sub_class"] == "Non-Syn"


# --------------------------------------------------------------------------- #
# minus-strand coding annotation (the verified synthetic case)
# --------------------------------------------------------------------------- #
def test_minus_strand_amino_acid_call():
    # Verified synthetic gene: reference "CCCCCTTAATCTTTCAT", minus-strand CDS
    # genome 6..17, whose mRNA (revcomp) is ATGAAAGATTAA -> M K D *.
    # A genome T->C at position 14 (mRNA A->G) changes residue 2 from K to E.
    ref = "CCCCCTTAATCTTTCAT"
    region = Region("cds", 6, 17, "-", 0, coding=True)
    row = annotate_in_region(region, 14, "T", "C", ref)
    assert row["residue"] == 2
    assert row["wt_aa"] == "K" and row["mut_aa"] == "E"
    assert row["sub_class"] == "Non-Syn"
    # nucleotide id stays genome/forward-strand
    assert row["mutation_id"] == "cds:T14C"
    # codon columns are the mRNA (revcomp) codons
    assert row["codon_wt"] == "AAA" and row["codon_mut"] == "GAA"


def test_minus_strand_residue_one_at_high_coordinate():
    # residue 1 sits at the feature END on minus strand (genome 17 here).
    ref = "CCCCCTTAATCTTTCAT"
    region = Region("cds", 6, 17, "-", 0, coding=True)
    row = annotate_in_region(region, 17, "T", "C", ref)   # start codon region
    assert row["residue"] == 1


# --------------------------------------------------------------------------- #
# non-coding annotation
# --------------------------------------------------------------------------- #
def test_noncoding_has_no_amino_acid():
    region = Region("5UTR", 1, 96, "+", 0, coding=False)
    row = annotate_in_region(region, 50, "A", "G", "N" * 100)
    assert row["region_type"] == "non-coding"
    assert row["wt_aa"] is None and row["mut_aa"] is None and row["residue"] is None
    assert row["mutation_id"] == "5UTR:A50G"
