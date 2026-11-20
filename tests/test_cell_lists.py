"""
test_cell_lists.py -- the master barcode list, custom cell lists against it,
and exact matching as the default.

Extract matches against the full whitelist; which barcodes are real cells is
applied afterwards from Cell Ranger / Seurat lists, which change between
analyses. The snakemake dry runs need snakemake and are skipped without it.
"""

from __future__ import annotations

import gzip
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from anchovy import cell_lists
from anchovy.cli import build_parser, main
from anchovy.config import ExtractConfig
from anchovy.schema import AnchovyColumns

REPO = Path(__file__).resolve().parent.parent
SNAKEFILE = REPO / "workflow" / "Snakefile"
CONFIG = REPO / "workflow" / "config.yaml"

needs_snakemake = pytest.mark.skipif(
    not shutil.which("snakemake"), reason="needs snakemake")

# Overhangs at barcode bases 8-9: GT = OB1, CA = OB2.
BC_OB1_A = "AAAAAAAGTAAAAAAA"
BC_OB1_B = "CCCCCCCGTCCCCCCC"
BC_OB2_A = "GGGGGGGCAGGGGGGG"
BC_NONE = "TTTTTTTTTTTTTTTT"


def _anchovy_csv(path: Path, rows):
    """An anchovy CSV with (CBC, UMI) rows; the other columns are filler."""
    frame = pd.DataFrame({c: ["x"] * len(rows) for c in AnchovyColumns.ORDER})
    frame[AnchovyColumns.CBC] = [r[0] for r in rows]
    frame[AnchovyColumns.UMI] = [r[1] for r in rows]
    frame.to_csv(path)
    return path


@pytest.fixture
def run_csv(tmp_path):
    return _anchovy_csv(tmp_path / "chan_anchovy.csv", [
        (BC_OB1_A, "U1"), (BC_OB1_A, "U1"), (BC_OB1_A, "U2"),
        (BC_OB1_B, "U3"),
        (BC_OB2_A, "U4"), (BC_OB2_A, "U5"),
        (BC_NONE, "U6"),
    ])


# --------------------------------------------------------------------------- #
# Exact matches are the default
# --------------------------------------------------------------------------- #
def test_exact_matching_is_the_default():
    assert ExtractConfig().max_barcode_errors == 0


def _extract_errors(argv):
    args = build_parser().parse_args(["extract", "x.sam", "wl.txt", *argv])
    from anchovy import extract
    captured = {}
    original = extract.run_to_csv
    extract.run_to_csv = lambda **kw: captured.update(kw) or 0
    try:
        args.func(args)
    finally:
        extract.run_to_csv = original
    return captured["config"].max_barcode_errors


def test_the_cli_keeps_the_mismatch_option():
    assert _extract_errors([]) == 0
    assert _extract_errors(["--max-barcode-errors", "1"]) == 1
    # The original unbounded nearest-match, spelled out now that 0 is default.
    assert _extract_errors(["--max-barcode-errors", "none"]) is None
    assert _extract_errors(["--max-barcode-errors", "None"]) is None
    with pytest.raises(SystemExit):
        _extract_errors(["--max-barcode-errors", "-1"])


def test_a_non_exact_barcode_is_dropped_by_default(data_dir, tmp_path):
    """A whitelist entry one substitution from the reads' barcode no longer
    claims them unless mismatches are asked for."""
    from anchovy import extract
    from make_golden import SIGNATURE

    original = (data_dir / "whitelist.txt").read_text().split()
    changed = [b[:-1] + ("A" if b[-1] != "A" else "C") for b in original]
    wl = tmp_path / "wl.txt"
    wl.write_text("\n".join(changed) + "\n")
    exact = extract.run(sam=str(data_dir / "test.sam"), whitelist=str(wl),
                        signature=SIGNATURE,
                        config=ExtractConfig(min_exact_rate=0, nthreads=1))
    assert len(exact) == 0
    one = extract.run(sam=str(data_dir / "test.sam"), whitelist=str(wl),
                      signature=SIGNATURE,
                      config=ExtractConfig(min_exact_rate=0, nthreads=1,
                                           max_barcode_errors=1))
    assert len(one) > 0


# --------------------------------------------------------------------------- #
# The master list
# --------------------------------------------------------------------------- #
def test_barcode_table_counts_reads_and_distinct_umis(run_csv):
    table = cell_lists.barcode_table(str(run_csv))
    assert list(table.columns) == ["CBC", "reads", "umis"]
    row = table.set_index("CBC").loc[BC_OB1_A]
    assert (row.reads, row.umis) == (3, 2)
    assert table.CBC.iloc[0] == BC_OB1_A, "most reads first"
    assert table.reads.sum() == 7


def test_an_ocm_channel_list_names_each_barcodes_sample(run_csv):
    from anchovy.ocm import resolve_mapping
    ocm = resolve_mapping({"S1": "OB1", "S2": "OB2"})
    table = cell_lists.barcode_table(str(run_csv), ocm=ocm).set_index("CBC")
    assert table.loc[BC_OB1_A, "ocm_sample"] == "S1"
    assert table.loc[BC_OB2_A, "ocm_sample"] == "S2"
    assert table.loc[BC_NONE, "ocm_sample"] == ""


def test_the_barcodes_command_writes_the_master_list(run_csv, tmp_path):
    out = tmp_path / "chan_barcodes.csv"
    assert main(["barcodes", str(run_csv), str(out),
                 "--ocm", "S1=OB1", "--ocm", "S2=CA"]) == 0
    table = pd.read_csv(out, keep_default_na=False)
    assert set(table.ocm_sample) == {"S1", "S2", ""}


# --------------------------------------------------------------------------- #
# Custom cell lists
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("content", [
    f"{BC_OB1_A}\n{BC_OB2_A}\n",
    f"{BC_OB1_A}-1\n{BC_OB2_A}-1\n",                     # Cell Ranger
    f'"x"\n"P10_5_{BC_OB1_A}-1"\n"P10_5_{BC_OB2_A}-1"\n',  # Seurat colnames
    f"barcode,cluster\n{BC_OB1_A}_1,3\n{BC_OB2_A}_1,4\n",  # Seurat CSV
])
def test_read_cell_list_accepts_what_cell_ranger_and_seurat_write(
        tmp_path, content):
    path = tmp_path / "cells.csv"
    path.write_text(content)
    assert cell_lists.read_cell_list(str(path)) == {BC_OB1_A, BC_OB2_A}


def test_read_cell_list_reads_gzip(tmp_path):
    path = tmp_path / "barcodes.tsv.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(f"{BC_OB1_A}-1\n")
    assert cell_lists.read_cell_list(str(path)) == {BC_OB1_A}


def test_read_cell_list_refuses_a_file_that_is_not_one(tmp_path):
    path = tmp_path / "genes.tsv"
    path.write_text("ENSG0001\tGAPDH\nENSG0002\tACTB\n")
    with pytest.raises(ValueError, match="no 16 nt barcode"):
        cell_lists.read_cell_list(str(path))


def test_match_cells_keeps_only_listed_barcodes(run_csv, tmp_path, capsys):
    master = tmp_path / "chan_barcodes.csv"
    cell_lists.write_barcode_table(str(run_csv), str(master))
    cells = tmp_path / "seurat_res08.txt"
    cells.write_text(f"{BC_OB1_A}-1\n{BC_OB2_A}-1\nACGTACGTACGTACGT-1\n")

    assert main(["match-cells", str(master), str(cells)]) == 0
    out = tmp_path / "chan_seurat_res08.csv"
    matched = pd.read_csv(out)
    assert set(matched.CBC) == {BC_OB1_A, BC_OB2_A}
    assert matched.reads.sum() == 5
    printed = capsys.readouterr().out
    assert "2 of 3 listed cells" in printed


def test_label_for_strips_the_extensions():
    assert cell_lists.label_for("/x/barcodes.tsv.gz") == "barcodes"
    assert cell_lists.label_for("res 0.8 cells.csv") == "res_0.8_cells"


# --------------------------------------------------------------------------- #
# The workflow
# --------------------------------------------------------------------------- #
def _dry_run(tmp_path, extra_yaml: str):
    extra = tmp_path / "extra.yaml"
    extra.write_text(extra_yaml)
    return subprocess.run(
        ["snakemake", "-s", str(SNAKEFILE), "-n", "-p", "--cores", "1",
         "--configfile", str(CONFIG), str(extra),
         "--config", f"results_dir={tmp_path / 'results'}", "report=false"],
        cwd=REPO, capture_output=True, text=True, timeout=120)


@needs_snakemake
def test_every_run_writes_a_master_list(tmp_path):
    proc = _dry_run(tmp_path, "{}\n")
    assert proc.returncode == 0, proc.stderr[-3000:]
    assert "anchovy barcodes" in proc.stdout
    assert "example_barcodes.csv" in proc.stdout


@needs_snakemake
def test_cell_lists_are_matched_against_the_master_list(tmp_path):
    cells = tmp_path / "res08.txt"
    cells.write_text("AAACCCAAGAAACACT-1\n")
    proc = _dry_run(tmp_path, f"cell_lists:\n  example: {{seurat: {cells}}}\n")
    assert proc.returncode == 0, proc.stderr[-3000:]
    assert "anchovy match-cells" in proc.stdout
    assert "cell_lists/example_seurat.csv" in proc.stdout


@needs_snakemake
def test_a_cell_list_for_an_unknown_sample_fails_before_any_job(tmp_path):
    proc = _dry_run(tmp_path, "cell_lists:\n  nobody: /nowhere.txt\n")
    assert proc.returncode != 0
    assert "not a sample" in proc.stderr + proc.stdout
