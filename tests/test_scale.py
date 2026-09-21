"""Scaling properties that only bite on real data.

Every one of these passes trivially on the test fixture and fails on a real
run, which is exactly why they are pinned here. The fixture has three cells;
a 10X run has hundreds of thousands, and three code paths that are invisible
at three cells become either a multi-day stall or a hard crash at that scale.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from anchovy.config import FastaConfig
from anchovy.fasta import build_cell_fastas
from anchovy.schema import AnchovyColumns as C

SNAKEFILE = Path(__file__).resolve().parents[1] / "workflow" / "Snakefile"


def _frame(n_barcodes: int, reads_per: int = 8) -> pd.DataFrame:
    barcodes = [f"BC{i:012d}" for i in range(n_barcodes)]
    rows = n_barcodes * reads_per
    return pd.DataFrame({
        C.CBC: np.repeat(barcodes, reads_per),
        C.READ: [f"r{i}" for i in range(rows)],
        C.MAPPED_SEQ: ["ACGT" * 20] * rows,
    })


def _reference_impl(df: pd.DataFrame, config: FastaConfig) -> dict[str, str]:
    """The original per-barcode-rescan implementation, kept as the oracle.

    build_cell_fastas was rewritten for speed, so the thing to defend is that
    the rewrite did not change a single byte of output. This is what it used
    to do.
    """
    out: dict[str, str] = {}
    for barcode in np.unique(df[C.CBC]):
        cell = df[df[C.CBC] == barcode]
        if len(cell) >= config.min_reads_per_cbc:
            out[str(barcode).strip()] = "\n".join(
                ">{}\n{}".format(row[C.READ], row[C.MAPPED_SEQ])
                for _, row in cell.iterrows())
    return out


# --------------------------------------------------------------------------- #
# build_cell_fastas: same answer, without the quadratic rescan
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("seed", range(4))
def test_grouping_matches_the_original_on_ragged_input(seed):
    """Byte-identical output, including key order, on messy realistic input.

    Ragged read counts straddling the threshold, shuffled rows, and
    whitespace-padded barcodes -- the three things that could make a groupby
    diverge from the original loop.
    """
    rng = np.random.default_rng(seed)
    barcodes, reads, seqs = [], [], []
    for i in range(int(rng.integers(20, 120))):
        pad = " " if i % 7 == 0 else ""
        for j in range(int(rng.integers(1, 12))):
            barcodes.append(f"{pad}BC{rng.integers(0, 50):04d}{pad}")
            reads.append(f"read{i}_{j}")
            seqs.append("".join(rng.choice(list("ACGT"), 30)))
    df = pd.DataFrame({C.CBC: barcodes, C.READ: reads, C.MAPPED_SEQ: seqs})
    df = df.sample(frac=1, random_state=seed).reset_index(drop=True)

    config = FastaConfig()
    expected = _reference_impl(df, config)
    actual = build_cell_fastas(df, config)

    assert actual == expected
    assert list(actual) == list(expected), "barcode order must stay sorted"


def test_grouping_handles_empty_and_all_thin_input():
    config = FastaConfig()
    empty = pd.DataFrame({C.CBC: [], C.READ: [], C.MAPPED_SEQ: []})
    assert build_cell_fastas(empty, config) == {}

    thin = pd.DataFrame({C.CBC: ["A", "A"], C.READ: ["r1", "r2"],
                         C.MAPPED_SEQ: ["AC", "GT"]})
    assert build_cell_fastas(thin, config) == {}


def test_barcode_grouping_scales_subquadratically():
    """Doubling the cells must not roughly quadruple the time.

    The original rescanned the whole table once per barcode, so cost grew as
    barcodes x reads. Timing is noisy, so this asserts only the shape -- a
    genuine quadratic comes in near 4.0 and the linear replacement near 2.0,
    which a threshold of 3.0 separates comfortably.
    """
    import time

    config = FastaConfig()
    small, large = _frame(2000), _frame(4000)

    build_cell_fastas(small, config)          # warm pandas' import-time caches

    t0 = time.perf_counter()
    build_cell_fastas(small, config)
    small_elapsed = time.perf_counter() - t0

    t0 = time.perf_counter()
    build_cell_fastas(large, config)
    large_elapsed = time.perf_counter() - t0

    if small_elapsed < 1e-3:                  # too fast to time meaningfully
        pytest.skip("runtime below timer resolution")
    assert large_elapsed / small_elapsed < 3.0, (
        f"doubling cells multiplied runtime by "
        f"{large_elapsed / small_elapsed:.1f}x, which looks quadratic again")


# --------------------------------------------------------------------------- #
# merge_consensus: the gather must not put every path on one command line
# --------------------------------------------------------------------------- #
def test_merge_consensus_does_not_shell_out_to_cat():
    """`cat {input}` dies with E2BIG once the cell count is large.

    It is the gather step, so it fails only after every per-cell job has run.
    """
    text = SNAKEFILE.read_text()
    rule = text.split("rule merge_consensus:")[1].split("rule ")[0]
    assert "cat {input}" not in rule, (
        "merge_consensus must not expand every cell path onto one command line")
    assert "copyfileobj" in rule, "expected the streaming concatenation"


def test_streaming_concatenation_survives_an_argv_sized_file_list(tmp_path):
    """The rule's own body, on a list too long to have passed to `cat`."""
    parts = []
    for i in range(40000):
        p = tmp_path / f"AF304458.1__{i:016d}_BC.fasta"
        p.write_text(f">cell{i}\nACGT\n")
        parts.append(str(p))

    argv_bytes = sum(len(p) + 1 for p in parts)
    assert argv_bytes > os.sysconf("SC_ARG_MAX") * 0.9, (
        "fixture no longer approaches ARG_MAX; raise the file count")

    out = tmp_path / "merged.fasta"
    with open(out, "wb") as out_fh:                 # the rule body, verbatim
        for path in parts:
            with open(path, "rb") as in_fh:
                shutil.copyfileobj(in_fh, out_fh)

    assert out.read_text().count(">") == 40000

    # And confirm the old form really would have failed on this same list.
    with pytest.raises(OSError):
        subprocess.run(["cat", *parts], stdout=subprocess.DEVNULL, check=True)


# --------------------------------------------------------------------------- #
# The knobs that control how much work a cluster run creates
# --------------------------------------------------------------------------- #
def test_fasta_checkpoint_exposes_min_reads():
    """min_reads decides how many per-cell jobs exist, so config must reach it."""
    text = SNAKEFILE.read_text()
    rule = text.split("checkpoint fasta:")[1].split("rule ")[0]
    assert "min_reads" in rule
    assert "--min-reads" in rule


def test_per_cell_rules_share_a_group():
    """Both per-cell rules must be in one group, for --group-components."""
    text = SNAKEFILE.read_text()
    for rule_name in ("rule map_cell:", "rule cell_consensus:"):
        body = text.split(rule_name)[1].split("\nrule ")[0]
        assert re.search(r'group:\s*\n\s*"cell"', body), (
            f"{rule_name} is not in the 'cell' group")


# --------------------------------------------------------------------------- #
# AnnotationConfig: the hypermutation cap must actually be reachable
# --------------------------------------------------------------------------- #
def _cells(n_mutations_per_cell: dict[str, int]) -> pd.DataFrame:
    """A filtConsensus-shaped frame where each cell carries N mutations."""
    rows = []
    for cbc, n in n_mutations_per_cell.items():
        genotype = "_".join(f"{i + 1}A" for i in range(n))
        rows.append({"CBC_ID": cbc, "genotype": genotype})
    return pd.DataFrame(rows)


def test_max_mutations_per_cell_is_honoured():
    """The field existed but nothing read it; a literal 200 was used instead."""
    from anchovy.annotate import haplo_analysis
    from anchovy.config import AnnotationConfig

    cons = _cells({"quiet": 2, "loud": 10})

    kept_by_default = set(haplo_analysis(cons)["CBC_ID"])
    assert {"quiet", "loud"} <= kept_by_default, "default 200 should keep both"

    capped = haplo_analysis(cons, config=AnnotationConfig(max_mutations_per_cell=5))
    mutation_rows = capped[capped["mutants"] != ""]
    assert "quiet" in set(mutation_rows["CBC_ID"])
    assert "loud" not in set(mutation_rows["CBC_ID"]), (
        "a cell above the cap must be dropped; the setting is not wired through")


def test_max_mutations_bound_is_exclusive():
    """R used `< 200`, so a cell AT the threshold is dropped. Preserved."""
    from anchovy.annotate import haplo_analysis
    from anchovy.config import AnnotationConfig

    cons = _cells({"at_threshold": 5, "below": 4})
    out = haplo_analysis(cons, config=AnnotationConfig(max_mutations_per_cell=5))
    with_mutations = set(out[out["mutants"] != ""]["CBC_ID"])

    assert "below" in with_mutations
    assert "at_threshold" not in with_mutations


def test_annotation_config_has_no_settings_that_control_nothing():
    """plot_haplotypes and build_network described behavior that did not exist."""
    from anchovy.config import AnnotationConfig

    fields = set(AnnotationConfig.__dataclass_fields__)
    assert "plot_haplotypes" not in fields, "the port dropped plotting"
    assert "build_network" not in fields, "duplicated annotate.run(network=)"
    assert fields == {"max_mutations_per_cell"}


# --------------------------------------------------------------------------- #
# extract: shared read-only state must not travel in every payload
# --------------------------------------------------------------------------- #
def test_assign_payloads_carry_no_bulk_objects():
    """The whitelist must not be pickled once per chunk.

    Pool.map pickles payloads, so a 737,280-barcode whitelist inside each one
    made a single payload 185 MB. Across the ~64 chunks that 402k reads over 16
    processes produces, that is ~12 GB serialized in each direction -- to let
    the worker do one positional lookup. The payload is now the per-read fields
    only, with the shared state published once by the pool initializer.
    """
    import inspect
    import pickle

    from anchovy import extract

    source = inspect.getsource(extract.assign_barcodes)
    payload_block = source.split("payloads = [")[1].split("]")[0]
    for leaked in ("whitelist", "barcode_blocks"):
        assert leaked not in payload_block, (
            f"{leaked} is back in the per-read payload; it belongs in "
            f"_assign_init, or every chunk ships a copy of it")

    assert "initializer=_assign_init" in source, (
        "the pool must publish shared state via its initializer")

    # A payload must stay small enough that per-chunk pickling is irrelevant.
    payload = ("read1", "ACGT" * 500, 5, "x" * 58, 2000, 0)
    assert len(pickle.dumps(payload)) < 10_000


def test_assign_worker_reads_its_state_from_the_initializer(tmp_path):
    """_assign_worker must work once _assign_init has populated the state."""
    import numpy as np
    import pandas as pd

    from anchovy import extract
    from anchovy.barcodes import build_barcode_query_blocks

    signature = "CTACACGACGCTCTTCCGATCT" + "N" * 26 + "TTTCTTATAT"
    whitelist = pd.DataFrame({"CBC": ["AAACCCAAGAAACACT", "AAACCCAAGAAACCAT"]})
    blocks = build_barcode_query_blocks(signature, whitelist.CBC)

    extract._assign_init(whitelist, blocks, 38, 10)

    # A read whose matched block carries the first barcode exactly.
    matchseq = str(blocks[0]).replace("N" * 10, "TGTGTTATCT")
    row = ("read0", "G" * 20 + matchseq, 20, matchseq, 20 + len(matchseq), 20)
    result = extract._assign_worker(row)

    assert result[0] == "AAACCCAAGAAACACT"   # CBC
    assert result[7] == "read0"               # read id
    assert isinstance(np.int64(result[1]), np.integer) or isinstance(result[1], int)
