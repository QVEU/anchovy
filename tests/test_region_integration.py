"""
test_region_integration.py -- tests for region awareness wired THROUGH the
pipeline, as opposed to test_regions.py which unit-tests the region model alone.

Three things are guarded here, in increasing scope:

1. WHOLE-REFERENCE MODE in the consensus stage (consensus.run(trim=False)):
   positions stay genome coordinates, the gap filter doesn't nuke every cell,
   and gap columns never produce substitution calls.

2. THE BACKWARD-COMPATIBILITY CLAIM: annotate.run(gff=...) must reproduce the
   legacy frame-1 amino-acid calls exactly when the GFF describes a single
   + strand CDS starting at genome position 1 with phase 0. If region awareness
   changed those calls, it would be a silent behavior change, not an addition.

3. END TO END on the mapping fixture: whole-reference consensus -> region-aware
   annotate, asserting the hand-verified prediction recorded in
   mapping_expected.json (genome 201 -> polyprotein residue 18).
"""

from __future__ import annotations

import json

import pandas as pd
import pytest

from anchovy import annotate, consensus
from anchovy.config import ConsensusConfig
from anchovy.consensus import ConsensusRecord, genotype_summary, select_sequences


def _rec(cbc_id, seq, coverage):
    return ConsensusRecord(cbc_id=cbc_id, seq=seq,
                           description=f"{cbc_id} ref coverage:{coverage}",
                           coverage=coverage, length=len(seq))


# --------------------------------------------------------------------------- #
# 1. Whole-reference mode in the consensus stage
# --------------------------------------------------------------------------- #
def test_whole_reference_keeps_full_length_and_genome_positions():
    """trim=False keeps sequences whole, so token positions ARE genome coords."""
    recs = [_rec("a", "AAAAAAAAAA", 50), _rec("b", "AAAATAAAAA", 50)]
    kept = select_sequences(recs, start=2, end=8, config=None, trim=False)
    assert [len(r.seq) for r in kept] == [10, 10]          # nothing was cut

    genos, ref = genotype_summary([r.seq for r in kept], window=(2, 8),
                                  skip_gaps=True)
    # The variant sits at 0-based index 4 -> genome position 5, NOT renumbered
    # relative to the window start.
    assert genos == ["", "5T"]
    assert len(ref) == 10


def test_legacy_trim_still_renumbers_within_region():
    """The contrast case: trim=True keeps the original renumbering behavior."""
    recs = [_rec("a", "AAAAAAAAAA", 50), _rec("b", "AAAATAAAAA", 50)]
    kept = select_sequences(recs, start=2, end=8)          # trim defaults True
    assert [len(r.seq) for r in kept] == [6, 6]            # cut to the region

    genos, _ = genotype_summary([r.seq for r in kept])
    # Same variant, now numbered from 1 within the trimmed region: index 4 of the
    # full sequence is index 2 of [2:8], so "3T".
    assert genos == ["", "3T"]


def test_gap_filter_skipped_when_no_window_given():
    """The bug that made whole-reference mode unusable.

    sam2consensus gap-fills uncovered flanks, so a whole-genome gap count is
    dominated by them and the max_gaps filter drops every cell. With no window
    there is nothing to judge gappiness against, so the filter must not run.
    """
    gappy = "-" * 20 + "ACGTACGTAC" + "-" * 20     # 40 gaps, way over max_gaps=3
    recs = [_rec("a", gappy, 50), _rec("b", gappy, 50)]

    assert select_sequences(recs, trim=False) != []        # no window -> kept
    # ...and with a window over the covered core, still kept (no gaps in there).
    assert select_sequences(recs, start=20, end=30, trim=False) != []
    # ...but a window over the ragged flank correctly drops them.
    assert select_sequences(recs, start=0, end=20, trim=False) == []


def test_no_variant_called_at_a_gap_position():
    """A coverage difference is not a substitution, in either direction."""
    # Column 0: cell 'b' has no coverage ('-'). Column 4: the reference (majority)
    # is a gap. Neither is a real substitution; only column 2 is.
    seqs = ["AAAAA", "-AGA-"]
    genos, _ = genotype_summary(seqs, reference="AAAA-", skip_gaps=True)
    assert genos == ["", "3G"]

    # Without skip_gaps (the legacy path) both spurious calls come back, which is
    # exactly why the flag exists and why it defaults to off.
    genos_legacy, _ = genotype_summary(seqs, reference="AAAA-")
    assert genos_legacy == ["5A", "1-_3G"]


# --------------------------------------------------------------------------- #
# 2. Region-aware annotation reproduces the legacy calls on a single CDS
# --------------------------------------------------------------------------- #
def _write_single_cds_gff(path, length):
    """A GFF3 describing one + strand CDS over the whole reference, phase 0.

    This is the region model that is EQUIVALENT to the legacy assumption, so it
    is the case where the two paths must agree.
    """
    path.write_text(
        "##gff-version 3\n"
        f"ref\tanchovy\tCDS\t1\t{length - (length % 3)}\t.\t+\t0\tID=cds;Name=cds\n")
    return str(path)


def test_annotate_region_matches_legacy_on_single_cds(tmp_path):
    """The backward-compatibility claim, checked call by call.

    A GFF3 whose only feature is a frame-1 CDS starting at position 1 encodes
    exactly the legacy assumption, so every amino-acid call must be identical.
    """
    reference = ("ATGAAAGATTTTCCCGGGAAATTTCCCTAAATGGGGCCCAAATTTGGGCCC" * 3)[:150]
    ref_file = tmp_path / "reference.txt"
    ref_file.write_text(reference)
    gff = _write_single_cds_gff(tmp_path / "single.gff3", len(reference))

    # A genotype table with several mutations spread across codon offsets.
    tokens = ["5A", "6G", "7C", "23T", "100G", "149A"]
    cons = tmp_path / "filtConsensus.csv"
    pd.DataFrame([
        {"CBC_ID": f"cell{i}", "genotype": tok, "sequence": reference,
         "description": "coverage:50"}
        for i, tok in enumerate(tokens)
    ]).to_csv(cons, index=False)

    legacy = annotate.run(str(cons), str(ref_file), str(tmp_path / "legacy"),
                          network=False)["annot"]
    region = annotate.run(str(cons), str(ref_file), str(tmp_path / "region"),
                          network=False, gff=gff)["annot"]

    cols = ["ref.codon", "ref.resPos", "ref.AA",
            "mut.codon", "mut.resPos", "mut.AA", "subName", "subClass"]
    left = legacy.set_index("mutants")[cols].sort_index()
    right = region.set_index("mutants")[cols].sort_index()
    pd.testing.assert_frame_equal(left, right, check_dtype=False)

    # And the calls are real, not all-NaN agreement.
    assert left["subClass"].isin(["Syn", "Non-Syn"]).all()


def test_region_aware_writes_long_format_for_overlapping_regions(tmp_path):
    """One mutation inside two nested regions yields one row per region."""
    reference = "ATG" + "AAA" * 40          # 123 nt, frame-1 CDS
    ref_file = tmp_path / "reference.txt"
    ref_file.write_text(reference)
    gff = tmp_path / "nested.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "ref\ta\tCDS\t1\t123\t.\t+\t0\tID=poly;Name=polyprotein\n"
        "ref\ta\tmature_protein_region\t61\t123\t.\t+\t0\tID=np;Name=NSP\n")

    cons = tmp_path / "filtConsensus.csv"
    pd.DataFrame([{"CBC_ID": "c1", "genotype": "64G", "sequence": reference,
                   "description": "coverage:50"}]).to_csv(cons, index=False)

    result = annotate.run(str(cons), str(ref_file), str(tmp_path / "out"),
                          network=False, gff=str(gff))
    regions_tbl = result["regions"]

    assert len(regions_tbl) == 2
    by_region = regions_tbl.set_index("region")
    # Position 64 is residue 22 of the polyprotein but residue 2 of the NSP --
    # each region numbered in its OWN frame of reference.
    assert by_region.loc["polyprotein", "residue"] == 22
    assert by_region.loc["NSP", "residue"] == 2
    # The nucleotide identity is genome-anchored and identical in both rows.
    assert set(by_region["genome_pos"]) == {64}

    # The legacy table is back-filled from the PRIMARY region, which is the first
    # one listed in the GFF3 (polyprotein), NOT the first alphabetically (NSP).
    assert result["annot"].iloc[0]["ref.resPos"] == 22


def test_region_aware_annotates_non_coding_without_amino_acids(tmp_path):
    """A UTR mutation gets a row with no amino-acid columns, not a bogus codon."""
    reference = "ACGT" * 50
    ref_file = tmp_path / "reference.txt"
    ref_file.write_text(reference)
    gff = tmp_path / "utr.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "ref\ta\tfive_prime_UTR\t1\t60\t.\t+\t.\tID=5UTR;Name=5UTR\n"
        "ref\ta\tCDS\t61\t180\t.\t+\t0\tID=cds;Name=polyprotein\n")

    cons = tmp_path / "filtConsensus.csv"
    pd.DataFrame([{"CBC_ID": "c1", "genotype": "10G", "sequence": reference,
                   "description": "coverage:50"}]).to_csv(cons, index=False)

    result = annotate.run(str(cons), str(ref_file), str(tmp_path / "out"),
                          network=False, gff=str(gff))
    row = result["regions"].iloc[0]
    assert row["region"] == "5UTR" and row["region_type"] == "non-coding"
    # reference is "ACGT" repeated, so genome position 10 is a 'C'.
    assert row["mutation_id"] == "5UTR:C10G"
    assert pd.isna(row["wt_aa"]) and pd.isna(row["residue"])

    # The legacy table still gets a usable name for it rather than a blank.
    assert result["annot"].iloc[0]["subName"] == "5UTR:C10G"


def test_region_annotations_csv_is_written_only_with_a_gff(tmp_path):
    reference = "ATG" + "AAA" * 20
    ref_file = tmp_path / "reference.txt"
    ref_file.write_text(reference)
    cons = tmp_path / "filtConsensus.csv"
    pd.DataFrame([{"CBC_ID": "c1", "genotype": "5G", "sequence": reference,
                   "description": "coverage:50"}]).to_csv(cons, index=False)

    legacy = annotate.run(str(cons), str(ref_file), str(tmp_path / "legacy"),
                          network=False)
    assert "regions" not in legacy["written"]
    assert not (tmp_path / "legacy_regionAnnotations.csv").exists()

    gff = tmp_path / "g.gff3"
    gff.write_text("##gff-version 3\nref\ta\tCDS\t1\t63\t.\t+\t0\tID=c;Name=c\n")
    annotate.run(str(cons), str(ref_file), str(tmp_path / "region"),
                 network=False, gff=str(gff))
    assert (tmp_path / "region_regionAnnotations.csv").exists()


# --------------------------------------------------------------------------- #
# 3. End to end on the mapping fixture
# --------------------------------------------------------------------------- #
def test_end_to_end_region_aware_on_mapping_fixture(data_dir, tmp_run_dir):
    """Whole-reference consensus -> region-aware annotate, against the prediction.

    The fixture plants TWO variants in different cells, so one pass exercises both
    halves of the region model:
      cellB  genome 201, inside the CDS   -> coding, frame-correct residue
      cellC  genome 121, inside the 5'UTR -> non-coding, NO amino-acid columns
    cellA matches the template and must come out clean.

    Starts from the sam2consensus output recorded in mapping_expected.json (itself
    verified against the real tool), so this runs without minimap2 installed while
    still exercising the real consensus and annotate code paths.
    """
    mapping = data_dir / "mapping"
    expected = json.loads((mapping / "mapping_expected.json").read_text())
    gff = mapping / "regions.gff3"
    if not gff.exists():
        pytest.skip("regions.gff3 missing; run tests/make_mapping_fixtures.py")

    # Rebuild the merged _allConsensus.fasta the workflow would produce.
    all_consensus = tmp_run_dir / "test_allConsensus.fasta"
    all_consensus.write_text("".join(
        f">{cell}_BC|c50 ref coverage:5.0 length:{len(seq)}\n{seq}\n"
        for cell, seq in expected["consensus"].items()))

    out_prefix = str(tmp_run_dir / "test")
    cons_result = consensus.run(
        fasta=str(all_consensus), start=100, end=500, trim=False,
        config=ConsensusConfig(depth_min=1),
        out_prefix=out_prefix)

    # The gap filter did not eat any cell, and the reference stayed whole.
    assert len(cons_result["records"]) == 3
    assert len(cons_result["reference"]) == expected["template_len"]

    # Each variant is called in its own cell, at its planted GENOME position --
    # not renumbered relative to the window or the ORF. With three cells the
    # majority at each column is the template base, so the computed reference
    # agrees with the template and the calls point the right way round.
    by_cell = {r["CBC_ID"].split("_")[0]: r["genotype"]
               for r in cons_result["records"]}
    assert by_cell["cellA"] == ""
    assert by_cell["cellB"] == (f"{expected['variant_genome_pos']}"
                                f"{expected['variant_alt_base']}")
    assert by_cell["cellC"] == (f"{expected['utr_variant_genome_pos']}"
                                f"{expected['utr_variant_alt_base']}")

    annot_result = annotate.run(
        filt_consensus_csv=f"{out_prefix}_filtConsensus.csv",
        reference_file=f"{out_prefix}_consensus_reference.txt",
        out_prefix=out_prefix, network=False, gff=str(gff))

    regions_tbl = annot_result["regions"].set_index("mutation_id")
    assert len(regions_tbl) == 2

    # --- the CODING variant ------------------------------------------------- #
    coding_id = (f"{expected['variant_region']}:{expected['variant_ref_base']}"
                 f"{expected['variant_genome_pos']}{expected['variant_alt_base']}")
    coding = regions_tbl.loc[coding_id]
    assert coding["region_type"] == "coding"
    assert coding["genome_pos"] == expected["variant_genome_pos"]
    # THE POINT: residue numbered from the CDS start (150), not from genome 1.
    assert coding["residue"] == expected["variant_residue"]
    assert coding["sub_class"] == "Non-Syn"

    # --- the NON-CODING variant --------------------------------------------- #
    # The behavior this half of the fixture exists for: a UTR mutation is
    # REPORTED (it is a real change, and the old annotator could not describe it
    # at all) but carries no amino-acid interpretation, rather than being
    # silently translated as though the UTR were coding.
    utr_id = (f"{expected['utr_variant_region']}:{expected['utr_variant_ref_base']}"
              f"{expected['utr_variant_genome_pos']}{expected['utr_variant_alt_base']}")
    utr = regions_tbl.loc[utr_id]
    assert utr["region_type"] == "non-coding"
    assert utr["genome_pos"] == expected["utr_variant_genome_pos"]
    for column in ("residue", "wt_aa", "mut_aa", "codon_wt", "codon_mut",
                   "sub_class"):
        assert pd.isna(utr[column]), f"{column} should be empty for a UTR variant"

    # --- the legacy table is still written, back-filled per primary region --- #
    legacy_tbl = annot_result["annot"].dropna(subset=["pos"]).set_index("mutants")
    coding_tok = f"{expected['variant_genome_pos']}{expected['variant_alt_base']}"
    utr_tok = f"{expected['utr_variant_genome_pos']}{expected['utr_variant_alt_base']}"
    assert legacy_tbl.loc[coding_tok, "ref.resPos"] == expected["variant_residue"]
    # The UTR mutation still gets a usable, genome-anchored name rather than a
    # blank that would collapse distinct genotypes together.
    assert legacy_tbl.loc[utr_tok, "subName"] == utr_id
    assert legacy_tbl.loc[utr_tok, "subClass"] == "non-coding"

    # --- and the coding call differs from the old frame-1 annotator ---------- #
    # which is the whole reason this feature exists.
    legacy_only = annotate.run(
        filt_consensus_csv=f"{out_prefix}_filtConsensus.csv",
        reference_file=f"{out_prefix}_consensus_reference.txt",
        out_prefix=str(tmp_run_dir / "legacy"), network=False)
    legacy_res = (legacy_only["annot"].dropna(subset=["pos"])
                  .set_index("mutants").loc[coding_tok, "ref.resPos"])
    assert legacy_res != expected["variant_residue"]


def test_mixed_region_types_in_one_run(tmp_path):
    """Coding, non-coding, phase lead-in and intergenic rows in a single table.

    REGRESSION: building one DataFrame from rows whose amino-acid fields are
    sometimes absent coerces those Nones to NaN, and `NaN is not None` is True.
    The legacy back-fill therefore has to test with pd.notna(); with a bare
    `is not None` check, the phase lead-in row below reaches int(NaN) and raises.
    """
    reference = "ACGT" * 60                    # 240 nt
    ref_file = tmp_path / "reference.txt"
    ref_file.write_text(reference)
    gff = tmp_path / "mixed.gff3"
    # CDS phase 2 => the reading frame starts at 62 + 2 = 64, so genome 62 and 63
    # sit in the lead-in and belong to no complete codon.
    gff.write_text(
        "##gff-version 3\n"
        "r\ta\tfive_prime_UTR\t1\t60\t.\t+\t.\tName=5UTR\n"
        "r\ta\tCDS\t62\t181\t.\t+\t2\tName=cds\n")

    cons = tmp_path / "filtConsensus.csv"
    pd.DataFrame([{"CBC_ID": "c1", "genotype": "10G_62T_70G_200A",
                   "sequence": reference,
                   "description": "coverage:50"}]).to_csv(cons, index=False)

    with pytest.warns(UserWarning, match="multiple of 3"):
        result = annotate.run(str(cons), str(ref_file), str(tmp_path / "out"),
                              network=False, gff=str(gff))

    by_id = result["regions"].set_index("mutation_id")
    assert by_id.loc["5UTR:C10G", "region_type"] == "non-coding"
    assert by_id.loc["cds:C62T", "sub_class"] == "non-coding-frame"
    assert pd.isna(by_id.loc["cds:C62T", "residue"])
    # genome 70 is frame position 70 - 64 + 1 = 7 -> residue 3
    assert by_id.loc["cds:C70G", "residue"] == 3
    assert by_id.loc["intergenic:T200A", "region"] == "intergenic"

    legacy = result["annot"].dropna(subset=["pos"]).set_index("mutants")
    assert legacy.loc["70G", "subName"] == "R3G"
    assert legacy.loc["10G", "subClass"] == "non-coding"
    assert legacy.loc["200A", "subName"] == "intergenic:T200A"


def test_validate_warnings_surface_from_run(tmp_path):
    """A suspicious region definition warns rather than silently mis-annotating."""
    reference = "ACGT" * 30
    ref_file = tmp_path / "reference.txt"
    ref_file.write_text(reference)
    gff = tmp_path / "bad.gff3"
    gff.write_text("##gff-version 3\nr\ta\tCDS\t1\t100\t.\t+\t0\tName=oops\n")
    cons = tmp_path / "filtConsensus.csv"
    pd.DataFrame([{"CBC_ID": "c1", "genotype": "10G", "sequence": reference,
                   "description": "coverage:50"}]).to_csv(cons, index=False)

    with pytest.warns(UserWarning, match="not a multiple of 3"):
        annotate.run(str(cons), str(ref_file), str(tmp_path / "out"),
                     network=False, gff=str(gff))


def test_minus_strand_end_to_end_through_annotate_run(tmp_path):
    """A minus-strand CDS driven all the way through annotate.run(gff=...).

    The codon math and the parser are covered separately in test_regions.py; this
    checks the whole path holds together for a minus-strand feature -- GFF3 on
    disk, genotype tokens in, region table and back-filled legacy columns out.

    Same verified gene: reference "CCCCCTTAATCTTTCAT", CDS at genome 6-17 on the
    minus strand, mRNA ATGAAAGATTAA (M K D *). Genome T14C is residue 2 K->E,
    confirmed by independently retranslating the mutated gene.
    """
    reference = "CCCCCTTAATCTTTCAT"
    ref_file = tmp_path / "reference.txt"
    ref_file.write_text(reference)
    gff = tmp_path / "minus.gff3"
    gff.write_text("##gff-version 3\n"
                   "ref\ta\tCDS\t6\t17\t.\t-\t0\tID=cds;Name=revgene\n")

    cons = tmp_path / "filtConsensus.csv"
    pd.DataFrame([{"CBC_ID": "c1", "genotype": "14C", "sequence": reference,
                   "description": "coverage:50"}]).to_csv(cons, index=False)

    result = annotate.run(str(cons), str(ref_file), str(tmp_path / "out"),
                          network=False, gff=str(gff))

    row = result["regions"].iloc[0]
    assert row["strand"] == "-"
    assert row["residue"] == 2
    assert row["wt_aa"] == "K" and row["mut_aa"] == "E"
    assert row["sub_class"] == "Non-Syn"
    # Genome-anchored id, forward-strand bases, even on a minus-strand feature.
    assert row["mutation_id"] == "revgene:T14C"
    assert row["wt_base"] == "T" and row["mut_base"] == "C"

    # The legacy back-fill names the residue in the FEATURE's frame, so a
    # minus-strand gene reads the same way a plus-strand one does.
    assert result["annot"].iloc[0]["subName"] == "K2E"

    # The written file carries it too.
    written = pd.read_csv(f"{tmp_path / 'out'}_regionAnnotations.csv")
    assert written.iloc[0]["mutation_id"] == "revgene:T14C"
    assert written.iloc[0]["strand"] == "-"
