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
    from anchovy.barcodes import (build_barcode_index,
                                  build_barcode_query_blocks)

    signature = "CTACACGACGCTCTTCCGATCT" + "N" * 26 + "TTTCTTATAT"
    whitelist = pd.DataFrame({"CBC": ["AAACCCAAGAAACACT", "AAACCCAAGAAACCAT"]})
    blocks = build_barcode_query_blocks(signature, whitelist.CBC)

    index = build_barcode_index(whitelist.CBC)
    extract._assign_init(whitelist, blocks, index, 38, 10)

    # A read whose matched block carries the first barcode exactly.
    matchseq = str(blocks[0]).replace("N" * 10, "TGTGTTATCT")
    row = ("read0", "G" * 20 + matchseq, 20, matchseq, 20 + len(matchseq), 20)
    result = extract._assign_worker(row)

    assert result[0] == "AAACCCAAGAAACACT"   # CBC
    assert result[7] == "read0"               # read id
    assert isinstance(np.int64(result[1]), np.integer) or isinstance(result[1], int)


# --------------------------------------------------------------------------- #
# read_sam: one streaming pass, memory bounded by survivors
# --------------------------------------------------------------------------- #
# read_sam had no direct test -- the extract golden covered it only end to end.
# These pin the properties the streaming rewrite had to preserve.
SAM_HEADER = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ref\tLN:1000\n"


def _sam(tmp_path, name, records, read_len=80):
    """Write a small SAM. Each record is (qname, flag, rname, cigar)."""
    seq, qual = "ACGT" * (read_len // 4), "I" * read_len
    lines = [SAM_HEADER]
    for qname, flag, rname, cigar in records:
        lines.append(f"{qname}\t{flag}\t{rname}\t1\t60\t{cigar}\t*\t0\t0\t"
                     f"{seq}\t{qual}\tNM:i:0\n")
    path = tmp_path / name
    path.write_text("".join(lines))
    return path


def test_read_sam_drops_unmapped_and_short(tmp_path, capsys):
    """Survivors are RNAME != '*' then length > min, in that order."""
    from anchovy.io import read_sam

    path = _sam(tmp_path, "mixed.sam", [
        ("keep1", 0, "ref", "80M"),
        ("unmapped", 4, "*", "*"),
        ("keep2", 0, "ref", "80M"),
    ])
    df = read_sam(str(path), 10)
    assert list(df["read"]) == ["keep1", "keep2"]

    # Every record is still counted, including the ones dropped.
    assert "Total Candidate Reads: 3" in capsys.readouterr().out

    # The length filter is strictly greater-than, as the original was.
    assert read_sam(str(path), 80).empty
    assert len(read_sam(str(path), 79)) == 2


def test_read_sam_keeps_sam_fields_as_strings(tmp_path):
    """Fields stay text, as the old split-based parse produced them.

    pysam exposes FLAG and POS as ints. Taking those natively would change the
    frame's dtypes without changing any value -- the kind of divergence that
    looks correct and breaks a downstream comparison.
    """
    from anchovy.io import read_sam
    from anchovy.schema import SamColumns

    df = read_sam(str(_sam(tmp_path, "types.sam", [("r", 0, "ref", "80M")])), 10)

    assert list(df.columns)[:len(SamColumns.ORDER)] == SamColumns.ORDER
    for column in SamColumns.ORDER:
        assert pd.api.types.is_string_dtype(df[column]), \
            f"{column} should stay text, got {df[column].dtype}"
        assert isinstance(df[column].iloc[0], str)
    # The two that pysam would hand back as ints if taken natively.
    assert df[SamColumns.FLAG].iloc[0] == "0"
    assert df[SamColumns.POS].iloc[0] == "1"
    for derived in (SamColumns.READ_LEN, SamColumns.CLIP_READ_LEN,
                    SamColumns.OFFSET, SamColumns.LENGTH):
        assert derived in df.columns


def test_read_sam_handles_an_unmapped_read_carrying_a_reference(tmp_path):
    """FLAG 4 with an RNAME set: legal SAM, and it used to abort the run.

    The old two passes disagreed here -- the text pass kept the row on
    RNAME != '*', the pysam pass dropped it on is_unmapped -- and zipping the
    two raised `ValueError: Length of values (N) does not match length of
    index (M)`. One pass cannot disagree with itself.
    """
    from anchovy.io import read_sam

    path = _sam(tmp_path, "placed.sam", [
        ("mapped", 0, "ref", "80M"),
        ("placed_unmapped", 4, "ref", "80M"),
    ])
    df = read_sam(str(path), 10)
    assert list(df["read"]) == ["mapped", "placed_unmapped"]


def test_read_sam_memory_scales_with_survivors_not_file_size(tmp_path):
    """The point of the rewrite: rejected reads are never retained.

    The old version listed every non-header line before filtering anything, so
    peak memory tracked the whole file. This asserts the shape of the fix --
    a file that is almost entirely rejects costs little more than its few
    survivors.
    """
    import tracemalloc

    from anchovy.io import read_sam

    survivors = [(f"keep{i}", 0, "ref", "400M") for i in range(20)]
    rejects = [(f"drop{i}", 4, "*", "*") for i in range(2000)]

    mostly_rejects = _sam(tmp_path, "rejects.sam", survivors + rejects,
                          read_len=400)
    just_survivors = _sam(tmp_path, "survivors.sam", survivors, read_len=400)

    assert mostly_rejects.stat().st_size > just_survivors.stat().st_size * 20

    def peak(path):
        tracemalloc.start()
        read_sam(str(path), 10)
        _, high = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        return high

    big, small = peak(mostly_rejects), peak(just_survivors)
    assert big < small * 3, (
        f"a file that is 99% rejects peaked at {big/1e6:.1f} MB against "
        f"{small/1e6:.1f} MB for its survivors alone -- rejected reads are "
        f"being retained again")


# --------------------------------------------------------------------------- #
# Barcode assignment: exact-match lookup instead of a 737,280-way scan
# --------------------------------------------------------------------------- #
# min_distance_block compares a read against EVERY whitelist block through
# np.vectorize, which is a Python loop: 0.47s per read against a v2 whitelist,
# or 52 core-hours for one real run. Every block differs from every other in
# only the 16 barcode characters, so an exact hit is a dict lookup.
TEST_SIGNATURE = "CTACACGACGCTCTTCCGATCT" + "N" * 26 + "TTTCTTATAT"


def _whitelist(n=400, seed=3):
    import numpy as np
    rng = np.random.default_rng(seed)
    return pd.DataFrame({"CBC": ["".join(c) for c in
                                 rng.choice(list("ACGT"), size=(n, 16))]})


def _read_block(barcode, umi="ACGTACGTAC"):
    return TEST_SIGNATURE[:22] + barcode + umi + TEST_SIGNATURE[-10:]


def test_exact_lookup_returns_what_the_full_scan_returns():
    """Same four-tuple, not merely the same barcode."""
    from anchovy.barcodes import (assign_barcode, build_barcode_index,
                                  build_barcode_query_blocks)

    wl = _whitelist()
    blocks = build_barcode_query_blocks(TEST_SIGNATURE, wl.CBC)
    index = build_barcode_index(wl.CBC)

    for i in (0, 7, 199, len(wl) - 1):
        block = _read_block(wl.CBC[i])
        assert (assign_barcode(block, blocks, wl, barcode_index=index)
                == assign_barcode(block, blocks, wl))


def test_exact_lookup_survives_damage_outside_the_barcode():
    """An error in the UMI or the 5' handle must not change the assignment."""
    from anchovy.barcodes import (assign_barcode, build_barcode_index,
                                  build_barcode_query_blocks)

    wl = _whitelist()
    blocks = build_barcode_query_blocks(TEST_SIGNATURE, wl.CBC)
    index = build_barcode_index(wl.CBC)

    intact = _read_block(wl.CBC[42])
    damaged_umi = intact[:38] + "TTTTTTTTTT" + intact[48:]
    damaged_handle = "G" + intact[1:]

    for block in (damaged_umi, damaged_handle):
        assert (assign_barcode(block, blocks, wl, barcode_index=index)
                == assign_barcode(block, blocks, wl))


def test_damaged_barcode_falls_back_to_the_scan():
    """No exact hit means the original path, not a wrong answer."""
    from anchovy.barcodes import (assign_barcode, build_barcode_index,
                                  build_barcode_query_blocks)

    wl = _whitelist()
    blocks = build_barcode_query_blocks(TEST_SIGNATURE, wl.CBC)
    index = build_barcode_index(wl.CBC)

    true_bc = wl.CBC[10]
    mutated = ("T" if true_bc[0] != "T" else "A") + true_bc[1:]
    assert mutated not in index, "fixture no longer exercises the fallback"

    block = _read_block(mutated)
    assert (assign_barcode(block, blocks, wl, barcode_index=index)
            == assign_barcode(block, blocks, wl))


def test_barcode_index_first_occurrence_wins():
    """Ties must break the way argmin breaks them, or the paths diverge."""
    from anchovy.barcodes import build_barcode_index

    duplicated = pd.DataFrame({"CBC": ["AAAACCCCGGGGTTTT", "ACGTACGTACGTACGT",
                                       "AAAACCCCGGGGTTTT"]})
    index = build_barcode_index(duplicated.CBC)
    assert index["AAAACCCCGGGGTTTT"] == 0, "argmin returns the first match"


def test_lookup_path_is_actually_taken():
    """Guards against the index being accepted and then quietly ignored."""
    from unittest.mock import patch

    from anchovy import barcodes
    from anchovy.barcodes import (assign_barcode, build_barcode_index,
                                  build_barcode_query_blocks)

    wl = _whitelist()
    blocks = build_barcode_query_blocks(TEST_SIGNATURE, wl.CBC)
    index = build_barcode_index(wl.CBC)

    with patch.object(barcodes, "min_distance_block",
                      side_effect=AssertionError("scanned despite an exact hit")):
        result = assign_barcode(_read_block(wl.CBC[5]), blocks, wl,
                                barcode_index=index)
    assert result[0] == wl.CBC[5]
    assert result[2] == 5


# --------------------------------------------------------------------------- #
# Bounding the barcode search, rather than filtering after it
# --------------------------------------------------------------------------- #
# Without a limit the search returns the NEAREST whitelist entry however far
# away it is, so a read whose barcode is unreadable still becomes a cell. The
# obvious fix -- assign, then discard by distance -- pays the full scan for
# every read it throws away. With a limit the question changes from "which is
# nearest" to "is anything within k errors", which the neighbourhood answers.
def test_bounded_and_unbounded_agree_when_bounded_answers():
    from anchovy.barcodes import (assign_barcode, build_barcode_index,
                                  build_barcode_query_blocks)

    wl = _whitelist()
    blocks = build_barcode_query_blocks(TEST_SIGNATURE, wl.CBC)
    index = build_barcode_index(wl.CBC)

    for i in (3, 44, 120):
        exact = _read_block(wl.CBC[i])
        one_off = _read_block(("T" if wl.CBC[i][0] != "T" else "A") + wl.CBC[i][1:])
        for block in (exact, one_off):
            bounded = assign_barcode(block, blocks, wl, barcode_index=index,
                                     max_barcode_errors=1)
            if bounded is not None:
                assert bounded == assign_barcode(block, blocks, wl,
                                                 barcode_index=index)


def test_a_barcode_beyond_the_limit_is_dropped_without_scanning():
    """The whole point: no exhaustive scan for a read that cannot pass."""
    from unittest.mock import patch

    from anchovy import barcodes
    from anchovy.barcodes import (assign_barcode, build_barcode_index,
                                  build_barcode_query_blocks)

    wl = _whitelist()
    blocks = build_barcode_query_blocks(TEST_SIGNATURE, wl.CBC)
    index = build_barcode_index(wl.CBC)

    junk = _read_block("TTTTTTTTTTTTTTTT")
    assert "TTTTTTTTTTTTTTTT" not in index

    with patch.object(barcodes, "min_distance_block",
                      side_effect=AssertionError("scanned a read it will drop")):
        assert assign_barcode(junk, blocks, wl, barcode_index=index,
                              max_barcode_errors=1) is None


def test_zero_errors_admits_only_exact_barcodes():
    from anchovy.barcodes import (assign_barcode, build_barcode_index,
                                  build_barcode_query_blocks)

    wl = _whitelist()
    blocks = build_barcode_query_blocks(TEST_SIGNATURE, wl.CBC)
    index = build_barcode_index(wl.CBC)

    assert assign_barcode(_read_block(wl.CBC[8]), blocks, wl,
                          barcode_index=index, max_barcode_errors=0)[0] == wl.CBC[8]

    one_off = ("T" if wl.CBC[8][0] != "T" else "A") + wl.CBC[8][1:]
    assert assign_barcode(_read_block(one_off), blocks, wl,
                          barcode_index=index, max_barcode_errors=0) is None


def test_unset_limit_leaves_the_original_behaviour():
    """No limit means the scan, and an answer for every read."""
    from anchovy.barcodes import (assign_barcode, build_barcode_index,
                                  build_barcode_query_blocks)

    wl = _whitelist()
    blocks = build_barcode_query_blocks(TEST_SIGNATURE, wl.CBC)
    index = build_barcode_index(wl.CBC)

    junk = _read_block("TTTTTTTTTTTTTTTT")
    assert assign_barcode(junk, blocks, wl, barcode_index=index) is not None
    assert assign_barcode(junk, blocks, wl) is not None


def test_substitution_neighbourhood_sizes_and_exclusivity():
    from anchovy.barcodes import _substitution_neighbours

    one = [n for _, n in _substitution_neighbours("A" * 16, 1)]
    assert len(one) == 16 * 3 == len(set(one))
    assert "A" * 16 not in one, "the unchanged barcode is not a neighbour"
    assert all(len(n) == 16 for n in one)


# --------------------------------------------------------------------------- #
# The extract CSV must be complete or absent, never half-written
# --------------------------------------------------------------------------- #
def test_anchovy_csv_is_written_atomically(tmp_path):
    """A crash mid-write must not leave a plausible-looking truncated file.

    to_csv streams in chunks straight to the destination, so a run killed
    partway -- OOM, wall clock -- left a partial CSV at exactly the name the
    fasta stage reads. Snakemake usually removes a failed job's outputs, but a
    SIGKILL gives it no opportunity.
    """
    from unittest.mock import patch

    from anchovy.io import write_anchovy_csv
    from anchovy.schema import AnchovyColumns

    frame = pd.DataFrame(
        [["AAACCCAAGAAACACT", 10, 0, 5, 0, "q", "m", "r1", "ACGT", "TTTT"]],
        columns=AnchovyColumns.ORDER)
    target = tmp_path / "sample_anchovy.csv"

    # A write that dies partway through, as an OOM kill would.
    real_to_csv = pd.DataFrame.to_csv

    def die_midway(self, path_or_buf=None, **kwargs):
        real_to_csv(self, path_or_buf, **kwargs)     # leave bytes behind
        raise MemoryError("killed mid-write")

    with patch.object(pd.DataFrame, "to_csv", die_midway):
        with pytest.raises(MemoryError):
            write_anchovy_csv(frame, str(target))

    assert not target.exists(), (
        "a failed write left a file where the next stage expects a complete one")
    assert not list(tmp_path.glob("*.partial")), "temp file was not cleaned up"

    # And the ordinary path still produces a readable, correctly ordered file.
    write_anchovy_csv(frame, str(target))
    assert target.exists()
    written = pd.read_csv(target)
    assert list(written.columns)[1:] == AnchovyColumns.ORDER
