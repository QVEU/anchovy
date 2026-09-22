"""
Tests for the example runner scripts.

THE BUG THESE GUARD. fetch.sh names the SAM it writes after the FASTQ, while
the workflow reads `sample` from the config. Nothing connected the two, so
pointing READS at a new FASTQ mapped it to its own name and then ran the
workflow against the PREVIOUS sample's SAM -- which was still sitting in
data_dir, so snakemake reported it up to date. The expensive mapping was
discarded and the results were the old sample's under the new sample's name.

Nothing failed, which is what made it worth a test rather than a comment.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

EXAMPLES = Path(__file__).resolve().parent.parent / "examples" / "eva71_sra"


def _sample_name(fastq: str) -> str:
    """Call the shell function the scripts actually use."""
    r = subprocess.run(
        ["bash", "-c",
         f'. "{EXAMPLES}/sample_name.sh"; sample_name_from_fastq "$1"',
         "_", fastq],
        capture_output=True, text=True, check=True)
    return r.stdout


@pytest.mark.parametrize("fastq,expected", [
    # The dataset this example is built around: the name keeps '.ccs', which is
    # what config_cluster.yaml's `sample` has always been, so the default run
    # is unchanged by the fix.
    ("/data/Sample_5/5_EVA71_6h_P5.ccs.fastq", "5_EVA71_6h_P5.ccs"),
    ("SRR28178313.fastq.gz", "SRR28178313"),
    ("reads.fq", "reads"),
    ("reads.fq.gz", "reads"),
    ("/a/b/sample_A.fastq", "sample_A"),
    ("/a/b/no_extension", "no_extension"),
])
def test_sample_name_derivation(fastq, expected):
    assert _sample_name(fastq) == expected


def test_run_cluster_passes_the_sample_to_the_workflow():
    """The fix itself: without this the workflow silently reads the old SAM."""
    script = (EXAMPLES / "run_cluster.sh").read_text()
    assert "--config sample=" in script, (
        "run_cluster.sh no longer overrides `sample` on the snakemake command "
        "line; a run with a new READS will read the previous sample's SAM and "
        "report it up to date")


def test_run_cluster_and_fetch_share_one_derivation():
    """Two copies of this rule is how they drifted apart in the first place."""
    for name in ("run_cluster.sh", "fetch.sh"):
        assert "sample_name.sh" in (EXAMPLES / name).read_text(), (
            f"{name} no longer sources the shared sample-name derivation")


def test_run_cluster_hands_the_same_name_to_fetch():
    # Passing it down means fetch.sh cannot pick a different name even if the
    # two derivations somehow diverge again.
    script = (EXAMPLES / "run_cluster.sh").read_text()
    assert 'SAMPLE="$SAMPLE"' in script


@pytest.mark.parametrize("name", ["run_cluster.sh", "fetch.sh", "sample_name.sh"])
def test_scripts_parse(name):
    subprocess.run(["bash", "-n", str(EXAMPLES / name)], check=True)
