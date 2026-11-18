"""
test_ocm.py -- splitting one On-Chip Multiplexing channel into its samples.

OCM writes each sample's identity into the cell barcode: a 2 nt overhang at
offset 7 of the 16 nt barcode (Cell Ranger chemistry_defs.json, every *-OCM
chemistry). These tests build a channel's anchovy CSV with known overhangs and
check that the split, the table read back from Cell Ranger's per-sample calls,
and the workflow wiring all agree on which read belongs to which sample.

The snakemake dry runs need snakemake and are skipped without it (CI has it).
"""

from __future__ import annotations

import gzip
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from anchovy import ocm
from anchovy.cli import main
from anchovy.schema import AnchovyColumns as C

REPO = Path(__file__).resolve().parent.parent
SNAKEFILE = REPO / "workflow" / "Snakefile"
CONFIG = REPO / "workflow" / "config.yaml"

needs_snakemake = pytest.mark.skipif(
    not shutil.which("snakemake"), reason="needs snakemake")


def _barcode(overhang: str, n: int) -> str:
    """A 16 nt barcode carrying `overhang` at positions 8-9."""
    stem = format(n, "014b").replace("0", "A").replace("1", "T")
    return stem[:7] + overhang + stem[7:]


def _channel_csv(path: Path, layout: dict[str, int]) -> pd.DataFrame:
    rows = []
    i = 0
    for over, n in layout.items():
        for _ in range(n):
            rows.append({
                C.CBC: _barcode(over, i), C.MIN_DISTANCE: 12, C.BARCODE_ID: i,
                C.READ_POSITION: 0, C.CLIP_LENGTH: 100, C.RECON_QUERY: "Q",
                C.MATCH_SEQ: "M", C.READ: f"read{i}", C.MAPPED_SEQ: "ACGT",
                C.UMI: "U" * 12,
            })
            i += 1
    df = pd.DataFrame(rows)[C.ORDER]
    df.to_csv(path)
    return df


def test_the_overhang_is_barcode_positions_8_and_9():
    assert ocm.overhang("AAAAAAACGTTTTTTT") == "CG"
    assert ocm.OVERHANG_OFFSET == 7 and ocm.OVERHANG_LEN == 2


def test_mapping_accepts_lists_and_strings():
    got = ocm.validate_mapping({"a": ["ac", "AG"], "b": "TT,TA"})
    assert got == {"a": {"AC", "AG"}, "b": {"TT", "TA"}}


@pytest.mark.parametrize("mapping, message", [
    ({}, "empty"),
    ({"a": []}, "no overhangs"),
    ({"a": ["ACG"]}, "not an overhang"),
    ({"a": ["AN"]}, "not an overhang"),
    ({"a": ["AC"], "b": ["AC"]}, "both"),
])
def test_a_bad_mapping_is_refused(mapping, message):
    with pytest.raises(ValueError, match=message):
        ocm.validate_mapping(mapping)


def test_split_keeps_exactly_one_samples_reads(tmp_path):
    csv = tmp_path / "chan_anchovy.csv"
    df = _channel_csv(csv, {"AC": 5, "AG": 3, "TT": 4, "GG": 2})

    out = tmp_path / "s1_anchovy.csv"
    kept, total = ocm.split_csv(str(csv), {"AC", "AG"}, str(out))
    assert (kept, total) == (8, 14)

    got = pd.read_csv(out, index_col=0)
    assert list(got.columns) == C.ORDER
    assert list(got.index) == list(range(8))
    assert set(got[C.CBC].str[7:9]) == {"AC", "AG"}
    expected = df[df[C.CBC].str[7:9].isin({"AC", "AG"})]
    assert list(got[C.READ]) == list(expected[C.READ])


def test_split_streams_across_chunks_without_losing_rows(tmp_path):
    csv = tmp_path / "chan_anchovy.csv"
    _channel_csv(csv, {"AC": 7, "TT": 6})
    out = tmp_path / "s_anchovy.csv"
    kept, total = ocm.split_csv(str(csv), {"AC"}, str(out), chunksize=3)
    got = pd.read_csv(out, index_col=0)
    assert (kept, total, len(got)) == (7, 13, 7)
    assert list(got.index) == list(range(7))


def test_a_sample_with_no_reads_still_gets_a_readable_table(tmp_path):
    csv = tmp_path / "chan_anchovy.csv"
    _channel_csv(csv, {"AC": 3})
    out = tmp_path / "empty_anchovy.csv"
    assert ocm.split_csv(str(csv), {"GG"}, str(out)) == (0, 3)
    got = pd.read_csv(out, index_col=0)
    assert got.empty and list(got.columns) == C.ORDER


def _per_sample_outs(root: Path, layout: dict[str, dict[str, int]]) -> Path:
    i = 0
    for sample, overs in layout.items():
        d = root / sample / "count" / "sample_filtered_feature_bc_matrix"
        d.mkdir(parents=True)
        lines = []
        for over, n in overs.items():
            for _ in range(n):
                lines.append(_barcode(over, i) + "-1")
                i += 1
        with gzip.open(d / "barcodes.tsv.gz", "wt") as handle:
            handle.write("\n".join(lines) + "\n")
    return root


def test_the_table_is_read_back_from_cellranger_per_sample_calls(tmp_path):
    root = _per_sample_outs(tmp_path / "per_sample_outs", {
        "P10_5": {"AC": 300, "AG": 290, "TT": 1},   # TT under 1%: noise
        "P10_8": {"GA": 400, "GC": 410},
    })
    assert ocm.table_from_cellranger(str(root)) == {
        "P10_5": ["AC", "AG"], "P10_8": ["GA", "GC"]}


def test_overlapping_samples_mean_this_is_not_ocm(tmp_path):
    root = _per_sample_outs(tmp_path / "per_sample_outs", {
        "A": {"AC": 100, "GG": 100}, "B": {"GG": 100, "TT": 100}})
    with pytest.raises(ValueError, match="both"):
        ocm.table_from_cellranger(str(root))


def test_cli_table_then_split(tmp_path, capsys):
    root = _per_sample_outs(tmp_path / "per_sample_outs", {
        "S1": {"AC": 50}, "S2": {"TT": 50}})
    assert main(["ocm-table", str(root)]) == 0
    out = capsys.readouterr().out
    assert "ocm:" in out and "S1: [AC]" in out and "S2: [TT]" in out

    csv = tmp_path / "chan_anchovy.csv"
    _channel_csv(csv, {"AC": 2, "TT": 5})
    dest = tmp_path / "S2_anchovy.csv"
    assert main(["ocm-split", str(csv), str(dest), "--overhangs", "TT",
                 "--sample", "S2"]) == 0
    assert "S2: kept 5 of 7 reads" in capsys.readouterr().out
    assert len(pd.read_csv(dest, index_col=0)) == 5


def test_the_workflow_routes_ocm_samples_through_the_split():
    text = SNAKEFILE.read_text()
    split = text.split("rule ocm_split:")[1].split("\nrule ")[0]
    extract = text.split("rule extract:")[1].split("\nrule ")[0]
    assert "anchovy ocm-split" in split
    assert "sample=OCM_SAMPLE_RE" in split
    assert "sample=NOT_OCM_SAMPLE_RE" in extract


# --------------------------------------------------------------------------- #
# snakemake dry runs
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
def test_an_ocm_channel_is_extracted_once_and_split_per_sample(tmp_path):
    # The example FASTQ is named "example"; treat it as a two-sample channel.
    proc = _dry_run(tmp_path, "ocm:\n  example:\n    s_one: [AC]\n"
                              "    s_two: [GT]\n")
    out = proc.stdout + proc.stderr
    assert proc.returncode == 0, out[-3000:]
    assert out.count("anchovy extract") == 1
    assert "--overhangs AC --sample s_one" in out
    assert "--overhangs GT --sample s_two" in out
    assert "results/s_one/s_one_annot_v3.csv" in out
    assert "example_annot_v3.csv" not in out


@needs_snakemake
def test_an_ocm_key_that_is_not_a_fastq_fails_up_front(tmp_path):
    proc = _dry_run(tmp_path, "ocm:\n  nosuch:\n    a: [AC]\n")
    assert proc.returncode != 0
    assert "no such FASTQ" in proc.stdout + proc.stderr
