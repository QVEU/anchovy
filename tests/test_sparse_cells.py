"""
test_sparse_cells.py -- cells with too little coverage must not kill the run.

THE FAILURE THIS GUARDS, which only appeared on real data. sam2consensus
collects the sequences it intends to write into a dict, skipping any whose bases
are all gaps, then writes one file per entry. A cell where NO position reaches
--min-depth therefore leaves that dict empty, writes nothing at all, and still
exits 0. The workflow declares that file as a required output, so one sparse
cell aborted the entire run with MissingOutputException.

On real 10X data most barcodes carry few reads, so this is the common case, not
an edge case. It never showed up on the mapping fixture because every fixture
cell has 6 reads at every position -- comfortably over the default min-depth of
5. A fixture that is uniformly well covered cannot exercise the path where
coverage runs out.

Dropping an under-covered cell is correct; that is what the depth filter is for.
The bug was only that "no consensus" is signalled by an absent file. The
workflow now touches the output, so the cell becomes an empty file and simply
does not appear downstream.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from anchovy.config import ConsensusConfig
from anchovy.consensus import parse_consensus_fasta, run

SAM2CONSENSUS = Path(__file__).resolve().parent.parent / "workflow" / "scripts" / "sam2consensus.py"


def _template(data_dir) -> str:
    path = data_dir / "mapping" / "template.fasta"
    if not path.exists():
        pytest.skip("mapping fixture missing; run tests/make_mapping_fixtures.py")
    return "".join(l.strip() for l in path.read_text().splitlines()
                   if not l.startswith(">"))


def _write_sam(path: Path, template: str, n_reads: int,
               start: int = 101, length: int = 100) -> None:
    """A SAM with n_reads identical alignments, i.e. uniform depth n_reads."""
    seq = template[start - 1:start - 1 + length]
    lines = [f"@SQ\tSN:testref\tLN:{len(template)}\n"]
    lines += [f"r{i}\t0\ttestref\t{start}\t60\t{length}M\t*\t0\t0\t"
              f"{seq}\t{'I' * length}\n" for i in range(n_reads)]
    path.write_text("".join(lines))


def _run_sam2consensus(sam: Path, outdir: Path, min_depth: int):
    outdir.mkdir(parents=True, exist_ok=True)
    return subprocess.run(
        [sys.executable, str(SAM2CONSENSUS), "-c", "0.5", "-m", str(min_depth),
         "--outfolder", f"{outdir}/", "-i", str(sam)],
        # check=False deliberately: the exit code is part of what these tests
        # assert, since the bug is that sam2consensus succeeds while producing
        # nothing.
        capture_output=True, text=True, check=False)


# --------------------------------------------------------------------------- #
# The upstream behavior we compensate for
# --------------------------------------------------------------------------- #
def test_sam2consensus_writes_nothing_when_depth_is_unmet(data_dir, tmp_path):
    """Characterizes the vendored tool: no output file, and exit code 0.

    If a future sam2consensus starts emitting an empty file itself, this fails
    and the `touch` in the workflow can be reconsidered. Until then it is what
    makes that touch necessary.
    """
    template = _template(data_dir)
    sam = tmp_path / "sparse_BC.sam"
    _write_sam(sam, template, n_reads=2)

    result = _run_sam2consensus(sam, tmp_path / "out", min_depth=5)

    assert result.returncode == 0, result.stderr
    assert list((tmp_path / "out").iterdir()) == [], (
        "sam2consensus produced a file; the workflow's touch may no longer be "
        "needed -- re-check before removing it")


def test_sam2consensus_writes_a_file_when_depth_is_met(data_dir, tmp_path):
    """The contrast case, so the test above is about DEPTH and not about setup."""
    template = _template(data_dir)
    sam = tmp_path / "covered_BC.sam"
    _write_sam(sam, template, n_reads=6)

    result = _run_sam2consensus(sam, tmp_path / "out", min_depth=5)

    assert result.returncode == 0, result.stderr
    produced = [p.name for p in (tmp_path / "out").iterdir()]
    assert produced == ["testref__covered_BC.fasta"]


# --------------------------------------------------------------------------- #
# Downstream tolerance of what the workflow leaves behind
# --------------------------------------------------------------------------- #
def test_empty_consensus_file_parses_to_no_records(tmp_path):
    """What `touch` leaves: an empty file, which must read as zero cells."""
    empty = tmp_path / "empty.fasta"
    empty.write_text("")
    assert parse_consensus_fasta(empty) == []


def test_merged_fasta_tolerates_cells_that_produced_nothing(tmp_path):
    """The real shape after a merge: sparse cells contribute no bytes at all.

    `cat` of an empty file adds nothing, so the merged FASTA simply holds fewer
    cells than there were barcodes. The consensus stage must treat that as
    normal rather than as malformed input.
    """
    merged = tmp_path / "test_allConsensus.fasta"
    # Two covered cells; three sparse ones contributed empty files to the cat.
    merged.write_text(
        ">cellA ref coverage:9 length:8\nACGTACGT\n"
        ">cellB ref coverage:9 length:8\nACGTACGA\n")

    result = run(fasta=str(merged), start=None, end=None, trim=False,
                 config=ConsensusConfig(depth_min=1),
                 out_prefix=str(tmp_path / "out"))

    assert [r["CBC_ID"] for r in result["records"]] == ["cellA", "cellB"]
    # cellB differs from cellA at the last base; with two cells the computed
    # reference breaks the tie alphabetically, so one of them carries the call.
    assert sum(1 for r in result["records"] if r["genotype"]) == 1


def test_a_merge_of_only_empty_files_does_not_crash(tmp_path):
    """Every cell sparse -- degenerate, but it must fail cleanly or empty."""
    merged = tmp_path / "test_allConsensus.fasta"
    merged.write_text("")

    with pytest.warns(UserWarning, match="passed filtering"):
        result = run(fasta=str(merged), start=None, end=None, trim=False,
                     config=ConsensusConfig(depth_min=1),
                     out_prefix=str(tmp_path / "out"))
    assert result["records"] == []


# --------------------------------------------------------------------------- #
# The workflow rule actually carries the fix
# --------------------------------------------------------------------------- #
def test_workflow_touches_the_consensus_output():
    """Guards the `touch` itself -- removing it reintroduces the abort."""
    snakefile = (Path(__file__).resolve().parent.parent
                 / "workflow" / "Snakefile").read_text()
    rule = snakefile.split("rule cell_consensus:")[1].split("\nrule ")[0]
    # `{output}` rather than `{output.fasta}`: the rule gained a second, optional
    # output (the pileup counts for the frequencies stage), and a sparse cell
    # writes neither file. Touching the whole output set covers both, and keeps
    # covering any output added later -- naming one field would silently stop
    # protecting the others.
    assert "touch {output}" in rule, (
        "cell_consensus no longer touches its outputs; a cell whose coverage "
        "never reaches --min-depth will abort the whole run")


# --------------------------------------------------------------------------- #
# When EVERY cell is sparse
# --------------------------------------------------------------------------- #
def test_annotate_reports_clearly_when_no_cells_survived(tmp_path):
    """Zero cells is a failed run, and must say so rather than crash obscurely.

    Fixing the missing-output abort moved this failure downstream: the workflow
    now reaches annotate, which used to die with KeyError: 'mutants' several
    stages after the real problem. It is a reachable state, not a curiosity --
    on a small slice of a real run the reads spread thinly over thousands of
    barcodes, so every cell can fall under the depth thresholds.
    """
    import pandas as pd

    from anchovy import annotate

    reference = tmp_path / "reference.txt"
    reference.write_text("ACGT" * 30)
    empty = tmp_path / "filtConsensus.csv"
    pd.DataFrame(columns=["CBC_ID", "genotype", "sequence",
                          "description"]).to_csv(empty, index=False)

    with pytest.raises(ValueError) as excinfo:
        annotate.run(str(empty), str(reference), str(tmp_path / "out"))

    message = str(excinfo.value)
    assert "no cells" in message
    # The message has to name the knobs, or it is just a nicer crash.
    assert "cons_min_depth" in message and "depth_min" in message


def test_haplo_analysis_survives_an_empty_table():
    """The pure function is total: an empty input gives an empty, well-shaped
    table rather than a frame with no columns at all."""
    import pandas as pd

    from anchovy.annotate import haplo_analysis

    out = haplo_analysis(pd.DataFrame(columns=["CBC_ID", "genotype"]))
    assert out.empty
    for column in ("mutants", "pos", "base", "CBC_ID", "genotype"):
        assert column in out.columns


def test_one_surviving_cell_is_enough(tmp_path):
    """The boundary next to the empty case: a single cell must still work."""
    import pandas as pd

    from anchovy import annotate

    reference = tmp_path / "reference.txt"
    reference.write_text("ATG" + "AAA" * 20)
    cons = tmp_path / "filtConsensus.csv"
    pd.DataFrame([{"CBC_ID": "only", "genotype": "5G", "sequence": "",
                   "description": "coverage:50"}]).to_csv(cons, index=False)

    result = annotate.run(str(cons), str(reference), str(tmp_path / "out"),
                          network=True)
    assert len(result["annot"]) == 1
    assert len(result["nodes"]) == 1


# --------------------------------------------------------------------------- #
# A run that finds no variants at all
# --------------------------------------------------------------------------- #
def test_annotate_survives_a_run_with_no_variants(tmp_path):
    """Zero variants is a RESULT, not an error, and must not crash.

    Taken from a real run: one cell cleared the depth filter, so the computed
    reference was that cell and nothing could differ from it. The annotation
    table was then built from an empty list of rows, giving a frame with no
    columns, and the merge on ["pos", "base"] died with KeyError: 'pos'.

    Same shape as the KeyError: 'mutants' above -- an empty intermediate losing
    its columns -- but in the other branch, so fixing one did not fix this.
    """
    import pandas as pd

    from anchovy import annotate

    reference = tmp_path / "reference.txt"
    reference.write_text("ATG" + "AAA" * 20)
    cons = tmp_path / "filtConsensus.csv"
    # One cell, empty genotype: exactly what a single-cell run produces.
    pd.DataFrame([{"CBC_ID": "only_BC|c50", "genotype": "", "sequence": "",
                   "description": "coverage:24.33"}]).to_csv(cons, index=False)

    result = annotate.run(str(cons), str(reference), str(tmp_path / "out"),
                          network=True)

    assert len(result["annot"]) == 1
    assert result["annot"].iloc[0]["genotype"] == "reference"
    # The node table still describes the one cell, honestly.
    nodes = result["nodes"]
    assert len(nodes) == 1
    assert nodes.iloc[0]["genotype"] == "reference"
    assert nodes.iloc[0]["nCells"] == 1
    assert nodes.iloc[0]["nMutations"] == 0


def test_no_variants_works_in_region_aware_mode_too(tmp_path):
    """The region-aware branch has its own empty path; check it as well."""
    import pandas as pd

    from anchovy import annotate

    reference = tmp_path / "reference.txt"
    reference.write_text("ACGT" * 50)
    gff = tmp_path / "r.gff3"
    gff.write_text("##gff-version 3\nr\ta\tCDS\t1\t60\t.\t+\t0\tName=cds\n")
    cons = tmp_path / "filtConsensus.csv"
    pd.DataFrame([{"CBC_ID": "only", "genotype": "", "sequence": "",
                   "description": "coverage:24"}]).to_csv(cons, index=False)

    result = annotate.run(str(cons), str(reference), str(tmp_path / "out"),
                          network=True, gff=str(gff))
    assert result["regions"].empty
    assert len(result["annot"]) == 1


def test_one_cell_with_a_computed_reference_warns(tmp_path):
    """The vacuity is worth saying out loud.

    With a computed reference the per-column consensus IS the single cell, so
    an empty genotype means "nothing to compare against", not "matches the
    virus". Those read identically in the output file.
    """
    from anchovy.consensus import run as consensus_run

    merged = tmp_path / "one_allConsensus.fasta"
    merged.write_text(">onecell ref coverage:24 length:8\nACGTACGT\n")

    with pytest.warns(UserWarning, match="nothing to compare against"):
        result = consensus_run(fasta=str(merged), trim=False,
                               config=ConsensusConfig(depth_min=1),
                               out_prefix=str(tmp_path / "out"))
    assert result["records"][0]["genotype"] == ""


def test_a_supplied_reference_makes_one_cell_meaningful(tmp_path):
    """...and the escape hatch the warning points at actually works.

    Against a supplied reference a single cell CAN carry variants, so no
    warning and a real genotype.
    """
    import warnings as _warnings

    from anchovy.consensus import run as consensus_run

    merged = tmp_path / "one_allConsensus.fasta"
    merged.write_text(">onecell ref coverage:24 length:8\nACGTACGT\n")

    with _warnings.catch_warnings():
        _warnings.simplefilter("error")          # any warning fails the test
        result = consensus_run(fasta=str(merged), trim=False,
                               reference="ACGTACGA",
                               config=ConsensusConfig(depth_min=1),
                               out_prefix=str(tmp_path / "out"))
    assert result["records"][0]["genotype"] == "8T"


# --------------------------------------------------------------------------- #
# Calling against a supplied reference rather than the crowd
# --------------------------------------------------------------------------- #
def test_a_single_cell_yields_variants_against_a_supplied_reference(tmp_path):
    """The fix for the single-cell dead end, end to end.

    A computed reference is the consensus ACROSS cells, so one cell means the
    reference IS that cell and no genotype can be non-empty -- regardless of
    what the cell actually carries. Against the genome, the same cell reports
    its variants.

    This is the difference between a run that "succeeded" with nothing in it and
    one that answers the question, so it is worth holding both halves.
    """
    from anchovy.consensus import run as consensus_run

    template = "".join("ACGTTGCA"[i % 8] for i in range(300))
    variant_pos = 150                                  # 1-based
    cell = (template[:variant_pos - 1]
            + ("A" if template[variant_pos - 1] != "A" else "C")
            + template[variant_pos:])

    merged = tmp_path / "one_allConsensus.fasta"
    merged.write_text(f">cell ref coverage:24 length:{len(cell)}\n{cell}\n")
    ref_fasta = tmp_path / "genome.fasta"
    ref_fasta.write_text(f">genome\n{template}\n")

    # Computed: nothing to compare against, so nothing is reported.
    with pytest.warns(UserWarning, match="nothing to compare against"):
        computed = consensus_run(fasta=str(merged), trim=False,
                                 config=ConsensusConfig(depth_min=1),
                                 out_prefix=str(tmp_path / "computed"))
    assert computed["records"][0]["genotype"] == ""

    # Supplied: the variant is found, at its genome position.
    from anchovy.io import read_reference_sequence
    supplied = consensus_run(fasta=str(merged), trim=False,
                             reference=read_reference_sequence(ref_fasta),
                             config=ConsensusConfig(depth_min=1),
                             out_prefix=str(tmp_path / "supplied"))
    genotype = supplied["records"][0]["genotype"]
    assert genotype == f"{variant_pos}{cell[variant_pos - 1]}"


def test_workflow_passes_the_reference_when_configured():
    """The Snakefile must forward `reference`, and track it as an input.

    As a params string alone it would be untracked, so changing the reference
    would leave stale genotypes in place without Snakemake noticing.
    """
    snakefile = (Path(__file__).resolve().parent.parent
                 / "workflow" / "Snakefile").read_text()
    rule = snakefile.split("rule consensus:")[1].split("\nrule ")[0]
    assert "--reference" in rule, "consensus rule does not forward `reference`"
    assert '"reference": REFERENCE' in rule, (
        "the reference is not declared as an input, so Snakemake cannot tell "
        "that changing it invalidates the genotypes")
