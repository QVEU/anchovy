"""
test_genbank_to_gff3.py -- the GenBank -> GFF3 converter used by the worked example.

WHY IT IS TESTED rather than left as example scaffolding: it exists specifically
to stop region coordinates being copied by hand, because a coordinate typo in a
GFF3 does not crash anything -- it silently renumbers every amino acid
downstream. A converter with that job has to be right, and its output has to be
something anchovy's own parser accepts.

The fixture is a miniature picornavirus laid out like EV-A71: 5'UTR, one
polyprotein CDS, mature peptides tiling it, 3'UTR.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from anchovy.regions import annotate_mutation, parse_gff3, validate_regions

SCRIPT = Path(__file__).resolve().parent.parent / "examples" / "eva71_sra" / "genbank_to_gff3.py"


@pytest.fixture(scope="module")
def converter():
    """Import the example script as a module (it is not part of the package)."""
    if not SCRIPT.exists():
        pytest.skip(f"missing {SCRIPT}")
    spec = importlib.util.spec_from_file_location("genbank_to_gff3", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules["genbank_to_gff3"] = module
    spec.loader.exec_module(module)
    return module


GENOME_LEN = 900
CDS_START, CDS_END = 101, 799          # 699 nt = 233 whole codons


def _feature(kind, start, end, strand=1, **qualifiers):
    return SeqFeature(SimpleLocation(start - 1, end, strand=strand),
                      type=kind, qualifiers=qualifiers)


@pytest.fixture
def picornavirus_gb(tmp_path):
    """A structurally faithful miniature picornavirus record."""
    record = SeqRecord(
        Seq("".join("ACGTTGCA"[i % 8] for i in range(GENOME_LEN))),
        id="AF000000.1", name="TESTEV", description="synthetic picornavirus")
    record.annotations["molecule_type"] = "ss-RNA"
    record.features = [
        _feature("5'UTR", 1, CDS_START - 1),
        _feature("CDS", CDS_START, CDS_END, gene=["polyprotein"], codon_start=["1"]),
        _feature("mat_peptide", 101, 310, product=["VP4"]),
        _feature("mat_peptide", 311, 520, product=["VP2"]),
        _feature("mat_peptide", 521, CDS_END, product=["VP1"]),
        _feature("3'UTR", CDS_END + 1, GENOME_LEN),
        # Must be ignored -- they carry no region anchovy annotates against.
        _feature("source", 1, GENOME_LEN, organism=["synthetic"]),
        _feature("gene", 1, GENOME_LEN, gene=["ignored"]),
    ]
    path = tmp_path / "record.gb"
    SeqIO.write(record, path, "genbank")
    return path, str(record.seq)


def _convert(converter, gb_path, out_path):
    assert converter.main(["genbank_to_gff3", str(gb_path), str(out_path)]) == 0
    return out_path


# --------------------------------------------------------------------------- #
# What it produces
# --------------------------------------------------------------------------- #
def test_converts_the_features_anchovy_annotates_against(converter, picornavirus_gb, tmp_path):
    gb, _ = picornavirus_gb
    gff = _convert(converter, gb, tmp_path / "regions.gff3")
    regions = {r.name: r for r in parse_gff3(str(gff))}

    assert set(regions) == {"five_prime_UTR_1", "polyprotein", "VP4", "VP2",
                            "VP1", "three_prime_UTR_800"}
    assert regions["polyprotein"].coding
    assert regions["VP1"].coding                    # mat_peptide -> coding
    assert not regions["five_prime_UTR_1"].coding


def test_skips_features_that_are_not_regions(converter, picornavirus_gb, tmp_path):
    """source and gene must not become regions, or they would shadow the CDS."""
    gb, _ = picornavirus_gb
    gff = _convert(converter, gb, tmp_path / "regions.gff3")
    names = [r.name for r in parse_gff3(str(gff))]
    assert "ignored" not in names
    assert "synthetic" not in names


def test_coordinates_survive_the_conversion(converter, picornavirus_gb, tmp_path):
    """Biopython is 0-based half-open, GFF3 is 1-based inclusive.

    Getting this wrong by one would shift every residue number silently, which
    is the whole failure mode the converter exists to prevent.
    """
    gb, _ = picornavirus_gb
    gff = _convert(converter, gb, tmp_path / "regions.gff3")
    regions = {r.name: r for r in parse_gff3(str(gff))}

    assert (regions["polyprotein"].start, regions["polyprotein"].end) == (CDS_START, CDS_END)
    assert (regions["VP4"].start, regions["VP4"].end) == (101, 310)
    assert (regions["three_prime_UTR_800"].start,
            regions["three_prime_UTR_800"].end) == (CDS_END + 1, GENOME_LEN)


def test_output_passes_anchovys_own_validation(converter, picornavirus_gb, tmp_path):
    """No coding region should come out a non-multiple of three."""
    gb, _ = picornavirus_gb
    gff = _convert(converter, gb, tmp_path / "regions.gff3")
    assert validate_regions(parse_gff3(str(gff))) == []


def test_overlapping_regions_are_usable_end_to_end(converter, picornavirus_gb, tmp_path):
    """The payoff: a mutation in VP1 is numbered in VP1 AND in the polyprotein."""
    gb, sequence = picornavirus_gb
    gff = _convert(converter, gb, tmp_path / "regions.gff3")
    regions = parse_gff3(str(gff))

    pos = 601                                        # inside VP1 (521-799)
    rows = {r["region"]: r for r in
            annotate_mutation(regions, pos, sequence[pos - 1], "A", sequence)}

    assert set(rows) == {"polyprotein", "VP1"}
    # polyprotein frame starts at 101: (601-101+1) = 501 -> residue 167
    assert rows["polyprotein"]["residue"] == 167
    # VP1 frame starts at 521: (601-521+1) = 81 -> residue 27
    assert rows["VP1"]["residue"] == 27


# --------------------------------------------------------------------------- #
# codon_start -> phase, and hostile input
# --------------------------------------------------------------------------- #
def test_codon_start_becomes_phase(converter, tmp_path):
    """GenBank codon_start is 1-based; GFF3 phase counts bases to SKIP.

    They differ by one, which is an easy off-by-one to inherit silently.
    """
    record = SeqRecord(Seq("ACGT" * 60), id="X.1", description="x")
    record.annotations["molecule_type"] = "ss-RNA"
    record.features = [_feature("CDS", 10, 72, gene=["p"], codon_start=["3"])]
    gb = tmp_path / "phase.gb"
    SeqIO.write(record, gb, "genbank")

    gff = _convert(converter, gb, tmp_path / "phase.gff3")
    assert parse_gff3(str(gff))[0].phase == 2        # codon_start 3 -> phase 2


def test_a_name_containing_gff3_delimiters_is_sanitized(converter, tmp_path):
    """';' and '=' would corrupt column 9 and silently split the attributes."""
    record = SeqRecord(Seq("ACGT" * 60), id="X.1", description="x")
    record.annotations["molecule_type"] = "ss-RNA"
    record.features = [_feature("CDS", 1, 60, product=["VP1; capsid=protein"])]
    gb = tmp_path / "nasty.gb"
    SeqIO.write(record, gb, "genbank")

    gff = _convert(converter, gb, tmp_path / "nasty.gff3")
    name = parse_gff3(str(gff))[0].name
    assert ";" not in name and "=" not in name
    assert name == "VP1, capsid-protein"


def test_minus_strand_feature_keeps_its_strand(converter, tmp_path):
    record = SeqRecord(Seq("ACGT" * 60), id="X.1", description="x")
    record.annotations["molecule_type"] = "ss-RNA"
    record.features = [_feature("CDS", 10, 72, strand=-1, gene=["rev"])]
    gb = tmp_path / "minus.gb"
    SeqIO.write(record, gb, "genbank")

    gff = _convert(converter, gb, tmp_path / "minus.gff3")
    assert parse_gff3(str(gff))[0].strand == "-"


# --------------------------------------------------------------------------- #
# Failure modes report clearly instead of raising a stack trace
# --------------------------------------------------------------------------- #
def test_unparseable_download_reports_clearly(converter, tmp_path, capsys):
    """A failed efetch can return an HTML error page with a 200 status."""
    bad = tmp_path / "bad.gb"
    bad.write_text("<html><body>Error: cannot get document</body></html>\n")
    assert converter.main(["genbank_to_gff3", str(bad), str(tmp_path / "o.gff3")]) == 1
    err = capsys.readouterr().err
    assert "not a readable GenBank record" in err or "no GenBank records" in err


def test_record_without_annotation_reports_clearly(converter, tmp_path, capsys):
    """A bare sequence with no features cannot drive region-aware annotation."""
    record = SeqRecord(Seq("ACGT" * 60), id="X.1", description="x")
    record.annotations["molecule_type"] = "ss-RNA"
    gb = tmp_path / "bare.gb"
    SeqIO.write(record, gb, "genbank")

    assert converter.main(["genbank_to_gff3", str(gb), str(tmp_path / "o.gff3")]) == 1
    assert "no CDS/mat_peptide/UTR features" in capsys.readouterr().err


def test_multi_record_file_is_refused(converter, tmp_path, capsys):
    """The pipeline assumes a single reference; a segmented genome must fail loudly."""
    records = []
    for name in ("A.1", "B.1"):
        r = SeqRecord(Seq("ACGT" * 60), id=name, description="x")
        r.annotations["molecule_type"] = "ss-RNA"
        r.features = [_feature("CDS", 1, 60, gene=["p"])]
        records.append(r)
    gb = tmp_path / "multi.gb"
    SeqIO.write(records, gb, "genbank")

    assert converter.main(["genbank_to_gff3", str(gb), str(tmp_path / "o.gff3")]) == 1
    assert "single reference" in capsys.readouterr().err
