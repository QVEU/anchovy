"""
Runs the shipped example end to end.

workflow/config.yaml is the command the README opens with, and the thing it is
supposed to prove -- that an install works -- is exactly what a structural check
on the fixture cannot prove. This runs it: minimap2, extract, per-cell
consensus, genotypes, networks and the frequency tables, against a fixture small
enough to finish in seconds.

SKIPPED when minimap2 or snakemake are absent, which is a dev container without
the conda environment. CI builds from environment.yml and has both, so the
example is verified there on every push.

The report is switched off: it needs R, it is the one stage with no assertable
output beyond "the file exists", and tests/test_report_rmd.py covers the
document itself. `--config report=false` is passed on the command line rather
than being a second config file, which also keeps the Snakefile honest about
coercing a string "false" back to a boolean.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
CONFIG = REPO / "workflow" / "config.yaml"

pytestmark = pytest.mark.skipif(
    not (shutil.which("minimap2") and shutil.which("snakemake")),
    reason="needs minimap2 and snakemake (conda env from environment.yml)")


@pytest.fixture(scope="module")
def results(tmp_path_factory) -> Path:
    """Run the example once, into a scratch results directory.

    results_dir keeps the run out of tests/data/fastqs, so a test run never
    leaves output beside the committed fixture and two runs cannot collide.
    """
    out = tmp_path_factory.mktemp("example_run")
    proc = subprocess.run(
        ["snakemake", "-s", str(REPO / "workflow" / "Snakefile"),
         "--configfile", str(CONFIG), "--cores", "2",
         "--config", "report=false", f"results_dir={out}"],
        cwd=REPO, capture_output=True, text=True)
    assert proc.returncode == 0, (
        f"the example in the README failed.\n--- stderr ---\n{proc.stderr[-4000:]}")
    return out / "example"


def test_report_false_on_the_command_line_is_honoured(results):
    """A string "false" is truthy in Python; the Snakefile has to coerce it."""
    assert not list(results.glob("*_report.html")), (
        "report=false was ignored, so every --config boolean override is being "
        "read as its opposite")


def test_every_cell_survives_and_is_named_by_its_barcode(results):
    consensus = (results / "example_allConsensus.fasta").read_text()
    names = [l.split()[0][1:] for l in consensus.splitlines() if l.startswith(">")]
    assert len(names) == 6, f"expected 6 cells, got {len(names)}: {names}"
    whitelist = set((REPO / "tests/data/fastqs/whitelist.txt").read_text().split())
    assert {n.split("_")[0] for n in names} == whitelist


def test_the_genotypes_are_the_ones_that_were_planted(results):
    annot = pd.read_csv(results / "example_annot_v3.csv")
    assert set(annot.genotype) == {"reference", "121C", "201A", "121C_201A"}
    # Three clean cells, one per single mutant, one double. Anything else means
    # a read was assigned to the wrong cell or a variant was called at the wrong
    # coordinate -- both of which produce a plausible-looking table.
    assert annot.genotype.value_counts()["reference"] == 3


def test_a_coding_and_a_non_coding_call_come_out_of_one_gff(results):
    regions = pd.read_csv(results / "example_regionAnnotations.csv")
    by_type = dict(zip(regions.region_type, regions.mutation_id))
    assert set(by_type) == {"coding", "non-coding"}, (
        "the region model is the reason to supply a GFF3; the example has to "
        "show both kinds of call")

    coding = regions[regions.region_type == "coding"].iloc[0]
    # Residue 18, not 67: the CDS starts at genome 150, so frame read from the
    # GFF is the only way to number this correctly. This is the specific thing
    # the pre-GFF annotator got wrong.
    assert (coding.genome_pos, coding.residue) == (201, 18)
    assert (coding.wt_aa, coding.mut_aa, coding.sub_class) == ("C", "S", "Non-Syn")

    non_coding = regions[regions.region_type == "non-coding"].iloc[0]
    assert non_coding.region == "5UTR"
    assert pd.isna(non_coding.residue), (
        "a UTR variant was given a residue number, i.e. translated as though "
        "it were coding")


def test_the_networks_have_an_edge_to_draw(results):
    nodes = pd.read_csv(results / "example_genotypeNodes.csv")
    assert set(nodes.genotype) == {"reference", "121C", "201A", "121C_201A"}
    assert nodes.set_index("genotype").nCells["reference"] == 3

    epi = pd.read_csv(results / "example_epistaticNetwork.csv")
    pairs = set(zip(epi.source, epi.target))
    assert ("121C", "121C_201A") in pairs or ("121C_201A", "121C") in pairs, (
        "the double mutant is isolated; the example shows nothing about "
        "epistasis")


def test_the_frequency_tables_use_a_per_position_denominator(results):
    pop = pd.read_csv(results / "example_alleleFrequencies.csv").set_index("position")
    assert sorted(pop.index) == [121, 201]
    for position in (121, 201):
        row = pop.loc[position]
        # Two of six cells carry each variant, and every cell covers every
        # position here, so the cell and read views must agree exactly. They
        # diverge only when depth differs between cells, which is the bias the
        # two columns exist to expose.
        assert (row.cells_alt, row.cells_called) == (2, 6)
        assert (row.reads_alt, row.reads_depth) == (16, 48)
        assert row.freq_cells == pytest.approx(row.freq_reads)

    per_cell = pd.read_csv(results / "example_cellAlleleFreq.csv")
    assert len(per_cell) == 4          # one single + one single + one double
    assert (per_cell.freq == 1.0).all(), (
        "every planted variant is on every read of its cell, so a within-cell "
        "frequency below 1.0 means reads leaked between cells")
