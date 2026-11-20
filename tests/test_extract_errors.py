"""
test_extract_errors.py -- failures that say why, and a scan that says what it
will cost.

A run that failed the exact-match floor reported only Snakemake's "exit status
1", with the explanation lost in a scheduler file. A rerun against the full
5' v3 whitelist with no `max_barcode_errors` then sat in the exhaustive scan
for a day with nothing on screen. These pin the fixes:

- the floor is checked before the expensive assignment, not after it;
- the scan prints its projected cost before starting;
- the CLI turns a stage's ValueError into a short marked message, exit 1;
- every anchovy stage in the workflow writes a log, and a misspelt config key
  is named instead of silently ignored.

The snakemake dry run needs snakemake and is skipped without it (CI has it).
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from anchovy import extract
from anchovy.cli import ERROR_MARKER, main
from anchovy.config import ExtractConfig

REPO = Path(__file__).resolve().parent.parent
SNAKEFILE = REPO / "workflow" / "Snakefile"
CONFIG = REPO / "workflow" / "config.yaml"

needs_snakemake = pytest.mark.skipif(
    not shutil.which("snakemake"), reason="needs snakemake")

SIGNATURE = "CTACACGACGCTCTTCCGATCT" + "N" * 26 + "TTTCTTATAT"


def _whitelist(n=400, seed=3):
    rng = np.random.default_rng(seed)
    return pd.DataFrame({"CBC": ["".join(c) for c in
                                 rng.choice(list("ACGT"), size=(n, 16))]})


def _chunk(barcodes):
    rows = []
    for i, bc in enumerate(barcodes):
        block = SIGNATURE[:22] + bc + "ACGTACGTAC" + SIGNATURE[-10:]
        rows.append({"read": f"r{i}", "seq": block + "A" * 40, "minPos": 0,
                     "matchseq": block, "readLen": len(block) + 40,
                     "offset": len(block), "minD": 10})
    return pd.DataFrame(rows)


class _PoolThatMustNotRun:
    def imap(self, *args, **kwargs):
        raise AssertionError("barcode assignment ran before the floor check")


def test_the_floor_fails_before_any_read_is_assigned():
    wl = _whitelist()
    config = ExtractConfig(nthreads=1, min_exact_rate=0.05)
    state = extract._BarcodeState(SIGNATURE, wl, config)
    off_list = ["T" * 16] * 50                     # none in the whitelist
    with pytest.raises(ValueError, match="EXACT entry"):
        extract._assign_chunk(_chunk(off_list), state, _PoolThatMustNotRun(),
                              extract._Tally(), first=True)


def test_later_chunks_are_not_rechecked():
    wl = _whitelist()
    config = ExtractConfig(nthreads=1, min_exact_rate=0.05,
                           max_barcode_errors=0)
    state = extract._BarcodeState(SIGNATURE, wl, config)

    class _Pool:
        def imap(self, func, payloads, chunksize=1):
            return [None for _ in payloads]

    rows = extract._assign_chunk(_chunk(["T" * 16] * 5), state, _Pool(),
                                 extract._Tally(), first=False)
    assert rows == []


def test_the_scan_prints_its_cost_before_it_starts(capsys):
    extract._warn_scan_cost(n_scan=90_000, n_reads=100_000,
                            whitelist_size=3_686_400, nthreads=16)
    out = capsys.readouterr().out
    assert "WARNING" in out and "max_barcode_errors" in out
    assert "hours" in out


def test_a_cheap_scan_says_nothing(capsys):
    extract._warn_scan_cost(n_scan=10, n_reads=100, whitelist_size=400,
                            nthreads=1)
    assert capsys.readouterr().out == ""


def test_a_stage_error_is_one_marked_message_and_exit_1(tmp_path, capsys):
    sam = tmp_path / "x.sam"
    sam.write_text("@HD\tVN:1.6\n")
    wl = tmp_path / "wl.txt"
    wl.write_text("AAAACCCCGGGGTTTT\n")
    rc = main(["extract", str(sam), str(wl), "--signature", "ACGTNNNN",
               "-o", str(tmp_path / "out.csv")])
    err = capsys.readouterr().err
    assert rc == 1
    assert f"{ERROR_MARKER}: anchovy extract stopped." in err
    assert "Traceback" not in err


def test_an_unexpected_error_keeps_its_traceback(monkeypatch, capsys):
    def boom(**kwargs):
        raise KeyError(0)

    monkeypatch.setattr(extract, "run_to_csv", boom)
    rc = main(["extract", "x.sam", "wl.txt"])
    err = capsys.readouterr().err
    assert rc == 1
    assert "Traceback" in err and "KeyError" in err
    assert f"{ERROR_MARKER}: anchovy extract crashed" in err


def test_every_anchovy_stage_writes_a_log():
    text = SNAKEFILE.read_text()
    for rule in ("rule map_reads:", "rule extract:", "checkpoint fasta:",
                 "rule ocm_split:", "rule consensus:", "rule annotate:",
                 "rule frequencies:"):
        body = text.split(rule)[1].split("\nrule ")[0].split("\ncheckpoint ")[0]
        assert "log:" in body, rule
        assert "{log}" in body, rule


def test_the_workflow_looks_for_the_cli_error_marker():
    assert f'_ERROR_MARKER = "{ERROR_MARKER}"' in SNAKEFILE.read_text()


@needs_snakemake
def test_a_misspelt_key_is_named(tmp_path):
    extra = tmp_path / "extra.yaml"
    extra.write_text("mismatches: 0\n")
    proc = subprocess.run(
        ["snakemake", "-s", str(SNAKEFILE), "-n", "--cores", "1",
         "--configfile", str(CONFIG), str(extra),
         "--config", f"results_dir={tmp_path / 'results'}", "report=false"],
        cwd=REPO, capture_output=True, text=True, timeout=120)
    out = proc.stdout + proc.stderr
    assert proc.returncode == 0, out[-3000:]
    assert "`mismatches` is not used" in out
    assert "max_barcode_errors" in out
    # The example config's own keys are all real.
    assert out.count("is not used by this workflow") == 1
