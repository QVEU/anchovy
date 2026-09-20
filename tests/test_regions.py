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


# A GFF3 mixing strands. Every other fixture in the suite is plus-strand only,
# so without this the parser's strand handling is never exercised from a FILE --
# the minus-strand tests below it build Region objects by hand, which skips
# parse_gff3 entirely. Strand and phase must come from each feature's own
# columns, never inferred from neighbours or nesting.
MIXED_STRAND_GFF3 = """\
##gff-version 3
ref\tanchovy\tCDS\t6\t17\t.\t-\t0\tID=rev;Name=revgene
ref\tanchovy\tCDS\t20\t31\t.\t+\t0\tID=fwd;Name=fwdgene
ref\tanchovy\tfive_prime_UTR\t1\t5\t.\t-\t.\tID=5UTR;Name=5UTR
"""


def test_parse_gff3_reads_strand_and_phase_per_feature(tmp_path):
    p = tmp_path / "mixed.gff3"
    p.write_text(MIXED_STRAND_GFF3)
    by_name = {r.name: r for r in parse_gff3(str(p))}

    assert by_name["revgene"].strand == "-"
    assert by_name["fwdgene"].strand == "+"
    # A minus-strand non-coding feature keeps its strand too, even though strand
    # does not affect how it is annotated.
    assert by_name["5UTR"].strand == "-" and not by_name["5UTR"].coding
    # No feature borrowed a neighbour's strand.
    assert [by_name[n].strand for n in ("revgene", "fwdgene", "5UTR")] == ["-", "+", "-"]


def test_parse_gff3_accepts_minus_strand_phase(tmp_path):
    """Phase is read for a minus-strand CDS as it is for a plus-strand one."""
    p = tmp_path / "phase.gff3"
    p.write_text("##gff-version 3\nref\ta\tCDS\t6\t17\t.\t-\t2\tID=x;Name=x\n")
    region = parse_gff3(str(p))[0]
    assert region.strand == "-" and region.phase == 2


def test_minus_strand_annotation_through_a_parsed_gff(tmp_path):
    """The verified minus-strand case, driven from a FILE rather than by hand.

    Same synthetic gene as test_minus_strand_amino_acid_call below: reference
    "CCCCCTTAATCTTTCAT", minus-strand CDS at genome 6-17, whose mRNA (the
    reverse complement) is ATGAAAGATTAA -> M K D *. A genome T->C at 14 is an
    mRNA A->G, turning residue 2 from K to E.

    The expected values were re-derived independently for this test by
    retranslating the whole mutated gene with Bio.Seq, not by reusing the codon
    helper under test: the only changed residue is 2, K->E.

    This closes the gap where minus-strand codon math was covered but only via
    hand-built Region objects, so parse_gff3 never produced a minus-strand
    feature that anything then annotated.
    """
    reference = "CCCCCTTAATCTTTCAT"
    p = tmp_path / "minus.gff3"
    p.write_text("##gff-version 3\nref\ta\tCDS\t6\t17\t.\t-\t0\tID=cds;Name=cds\n")

    region = parse_gff3(str(p))[0]
    assert validate_regions([region]) == []          # 12 nt, a clean 4 codons

    row = annotate_in_region(region, 14, "T", "C", reference)
    assert row["residue"] == 2
    assert row["wt_aa"] == "K" and row["mut_aa"] == "E"
    assert row["codon_wt"] == "AAA" and row["codon_mut"] == "GAA"
    assert row["sub_class"] == "Non-Syn"
    # The nucleotide identity stays genome/forward-strand even though the amino
    # acid columns are read in the minus-strand direction.
    assert row["mutation_id"] == "cds:T14C"

    # Parsing the feature must give exactly what building it by hand gives.
    hand_built = annotate_in_region(
        Region("cds", 6, 17, "-", 0, coding=True), 14, "T", "C", reference)
    assert row == hand_built


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


# --------------------------------------------------------------------------- #
# The example in the README must actually work
# --------------------------------------------------------------------------- #
def test_readme_gff3_examples_parse():
    """Every GFF3 example in the README parses, with tabs intact.

    Users are told to copy these and edit the numbers, so an example that has
    rotted -- or, far more likely, has had its tabs silently converted to spaces
    by an editor -- is a real bug in the documentation. GFF3 is tab-separated and
    tabs are invisible on screen, so nothing but a parse will catch it.
    """
    import re
    from pathlib import Path

    readme = Path(__file__).resolve().parent.parent / "README.md"
    if not readme.exists():                       # pragma: no cover
        pytest.skip("README.md not found")

    blocks = re.findall(r"```\n(##gff-version 3\n.*?)```", readme.read_text(),
                        re.DOTALL)
    assert blocks, "no GFF3 example found in README.md -- did the section move?"

    for i, block in enumerate(blocks):
        feature_lines = [ln for ln in block.splitlines()
                         if ln and not ln.startswith("#")]
        assert feature_lines, f"README GFF3 example {i} has no feature rows"
        for line in feature_lines:
            assert "\t" in line, (
                f"README GFF3 example {i} lost its tabs -- this row is "
                f"space-separated and will not parse: {line!r}")

        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / f"readme_{i}.gff3"
            path.write_text(block)
            regions = parse_gff3(str(path))

        assert regions, f"README GFF3 example {i} parsed to zero regions"
        # A documented example should not itself trip the coordinate sanity
        # check -- it is what people copy as a starting point.
        assert validate_regions(regions) == [], (
            f"README GFF3 example {i} triggers a validation warning")
