"""
Tests for the workflow's input contract: a folder of FASTQs.

WHAT THESE REPLACED. The pipeline used to need a script run before it --
fetch.sh mapped the reads to a SAM and run_cluster.sh derived the sample name
to hand over -- and the earlier tests here guarded that handover, because the
two sides disagreed the moment a run was pointed at different reads: the reads
were mapped under one name, the workflow looked for another, found the previous
run's SAM still in place, and reported it up to date without failing.

Mapping is a pipeline stage now and the sample name is derived once, in the
package, from the FASTQ itself. There is no handover left to get wrong, so the
tests moved to the thing that replaced it.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from anchovy.io import discover_fastqs, sample_name_from_fastq

REPO = Path(__file__).resolve().parent.parent
SNAKEFILE = (REPO / "workflow" / "Snakefile").read_text()


# --------------------------------------------------------------------------- #
# What a sample is called
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("fastq,expected", [
    # The dataset the example is built around: '.ccs' is part of the name, and
    # is what config_cluster.yaml's outputs have always been called.
    ("/data/Sample_5/5_EVA71_6h_P5.ccs.fastq", "5_EVA71_6h_P5.ccs"),
    ("SRR28178313.fastq.gz", "SRR28178313"),
    ("reads.fq", "reads"),
    ("reads.fq.gz", "reads"),
    ("/a/b/sample_A.fastq", "sample_A"),
    ("/a/b/no_extension", "no_extension"),
])
def test_sample_name_from_fastq(fastq, expected):
    assert sample_name_from_fastq(fastq) == expected


# --------------------------------------------------------------------------- #
# Discovery
# --------------------------------------------------------------------------- #
def _fq(d: Path, name: str) -> Path:
    p = d / name
    if name.endswith(".gz"):
        with gzip.open(p, "wt") as fh:
            fh.write("@r\nACGT\n+\nIIII\n")
    else:
        p.write_text("@r\nACGT\n+\nIIII\n")
    return p


def test_discovers_every_flavour_of_fastq(tmp_path):
    for name in ("a.fastq", "b.fastq.gz", "c.fq", "d.fq.gz"):
        _fq(tmp_path, name)
    assert sorted(discover_fastqs(tmp_path)) == ["a", "b", "c", "d"]


def test_ignores_everything_that_is_not_a_fastq(tmp_path):
    # input_dir doubles as the place the reference and GFF may sit, so
    # discovery has to be indifferent to them.
    _fq(tmp_path, "reads.fastq")
    (tmp_path / "AF304458.fasta").write_text(">ref\nACGT\n")
    (tmp_path / "AF304458.gff3").write_text("##gff-version 3\n")
    (tmp_path / "notes.txt").write_text("hello")
    assert list(discover_fastqs(tmp_path)) == ["reads"]


def test_does_not_recurse(tmp_path):
    # A nested layout is far likelier to be a mistake about what input_dir
    # means than an intent to run a whole tree, so it is not silently included.
    nested = tmp_path / "sub"
    nested.mkdir()
    _fq(nested, "deep.fastq")
    _fq(tmp_path, "top.fastq")
    assert list(discover_fastqs(tmp_path)) == ["top"]


def test_colliding_sample_names_are_an_error(tmp_path):
    """Their results would overwrite each other rather than collide."""
    _fq(tmp_path, "reads.fastq")
    _fq(tmp_path, "reads.fq.gz")
    with pytest.raises(ValueError, match="reduce to the sample name"):
        discover_fastqs(tmp_path)


def test_an_empty_directory_says_what_it_looked_for(tmp_path):
    with pytest.raises(ValueError, match="no FASTQs in input_dir"):
        discover_fastqs(tmp_path)


# --------------------------------------------------------------------------- #
# The workflow wiring these replaced a script for
# --------------------------------------------------------------------------- #
def test_mapping_is_a_pipeline_stage():
    assert "rule map_reads:" in SNAKEFILE, (
        "mapping left the workflow; running your own data needs a script again")
    assert "minimap2 -ax" in SNAKEFILE


def test_the_workflow_discovers_samples_rather_than_being_told_one():
    assert "discover_fastqs" in SNAKEFILE
    assert 'config["sample"]' not in SNAKEFILE
    assert 'config["data_dir"]' not in SNAKEFILE


def test_the_report_is_a_default_target():
    assert "rule report:" in SNAKEFILE
    assert "REPORT" in SNAKEFILE.split("rule all:")[1].split("rule ")[0]


def test_the_whitelist_is_fetched_and_verified():
    assert "rule get_whitelist:" in SNAKEFILE
    # The count check is the point: a wrong whitelist does not fail on its own.
    assert "barcodes, got" in SNAKEFILE


def test_fetch_no_longer_maps():
    fetch = (REPO / "examples" / "eva71_sra" / "fetch.sh").read_text()
    assert "minimap2" not in fetch, (
        "fetch.sh maps again; that is the workflow's job and the two disagreed "
        "about the sample name when both did it")
