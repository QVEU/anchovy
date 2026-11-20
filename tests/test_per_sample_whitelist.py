"""
test_per_sample_whitelist.py -- a reachable min_exact_rate, and the per-sample
whitelist mapping that `cell_lists` replaced.

Samples with their own Cell Ranger cell calls each need their own list, and a
list that covers only some of a library's cells can sit below the exact-match
floor legitimately. Before this, `whitelist` took one path for every sample and
the floor's own error message pointed at a `min_exact_rate` setting that
neither the CLI nor the workflow could reach.

The snakemake dry runs need snakemake and are skipped without it (CI has it);
the rest run everywhere.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from anchovy import extract
from anchovy.cli import build_parser
from anchovy.config import ExtractConfig

from make_golden import SIGNATURE

REPO = Path(__file__).resolve().parent.parent
SNAKEFILE = REPO / "workflow" / "Snakefile"
CONFIG = REPO / "workflow" / "config.yaml"

needs_snakemake = pytest.mark.skipif(
    not shutil.which("snakemake"), reason="needs snakemake")


# --------------------------------------------------------------------------- #
# min_exact_rate is reachable
# --------------------------------------------------------------------------- #
def _extract_config(argv):
    """Run the CLI's extract handler up to the point it builds its config."""
    args = build_parser().parse_args(["extract", "x.sam", "wl.txt", *argv])
    captured = {}

    def fake_run_to_csv(**kwargs):
        captured.update(kwargs)
        return 0

    original = extract.run_to_csv
    extract.run_to_csv = fake_run_to_csv
    try:
        args.func(args)
    finally:
        extract.run_to_csv = original
    return captured["config"]


def test_min_exact_rate_flag_reaches_the_config():
    assert _extract_config(["--min-exact-rate", "0.01"]).min_exact_rate == 0.01
    assert _extract_config(["--min-exact-rate", "0"]).min_exact_rate == 0


def test_min_exact_rate_defaults_when_the_flag_is_absent():
    assert (_extract_config([]).min_exact_rate
            == ExtractConfig().min_exact_rate)


def test_the_extract_rule_passes_min_exact_rate_and_the_whitelist():
    text = SNAKEFILE.read_text()
    rule = text.split("rule extract:")[1].split("\nrule ")[0]
    assert "--min-exact-rate" in rule
    assert "whitelist=whitelist_for" in rule


def test_the_floor_error_names_a_setting_that_exists(data_dir, tmp_path):
    # A list none of the fixture's reads carry: every read misses exactly.
    wl = tmp_path / "other.txt"
    wl.write_text("\n".join(["ACGTACGTACGTACGT", "TTTTGGGGCCCCAAAA"]) + "\n")
    with pytest.raises(ValueError) as err:
        extract.run(sam=str(data_dir / "test.sam"), whitelist=str(wl),
                    signature=SIGNATURE)
    message = str(err.value)
    assert "--min-exact-rate" in message
    assert "`min_exact_rate`" in message

    # And the setting it names does let the run through.
    df = extract.run(sam=str(data_dir / "test.sam"), whitelist=str(wl),
                     signature=SIGNATURE,
                     config=ExtractConfig(min_exact_rate=0,
                                          max_barcode_errors=2))
    assert len(df) == 0


# --------------------------------------------------------------------------- #
# whitelist as a sample -> path mapping
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
def test_a_per_sample_mapping_is_refused_and_points_at_cell_lists(tmp_path):
    """Matching against called cells was replaced by matching the full list
    and applying `cell_lists` afterwards. An old config must say so, not run."""
    wl = tmp_path / "example_cells.txt"
    shutil.copy(REPO / "tests/data/fastqs/whitelist.txt", wl)
    proc = _dry_run(tmp_path, f"whitelist:\n  example: {wl}\n")
    assert proc.returncode != 0
    assert "cell_lists" in proc.stderr + proc.stdout


@needs_snakemake
def test_min_exact_rate_reaches_the_extract_command(tmp_path):
    proc = _dry_run(tmp_path, "min_exact_rate: 0.01\n")
    assert proc.returncode == 0, proc.stderr[-3000:]
    assert "--min-exact-rate 0.01" in proc.stdout
