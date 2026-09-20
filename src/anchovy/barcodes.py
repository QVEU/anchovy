"""
barcodes.py -- the barcode/UMI matching logic (pure, testable core).

WHY THIS MODULE EXISTS SEPARATELY
---------------------------------
In the original anchovy.py, the distance-matching math was tangled together with
multiprocessing Pool calls and DataFrame column plumbing. That made the actual
logic -- "given a read window and a query, find the closest block" -- impossible
to test without spinning up a process pool and building a DataFrame.

Here the pure functions live on their own: sequence in, result out, no pools, no
files, no globals. extract.py handles the parallelism and I/O and calls into these.
This is the same split we used for io.py's pure functions, and it's why they can
have fast, deterministic unit tests.

BEHAVIOR PRESERVATION
---------------------
These reproduce the original blockDist / cellMatch math exactly, including:
  - the fuzzy min-distance-block search,
  - the barcode-augmented query construction in cellIDPool,
  - the UMI slice matchseq[umi_start : len-umi_end_trim].
The one change is the pandas compatibility fix: the original did
CBCs.iloc[minPos][0] (which breaks on pandas 2.x label indexing); here we use
positional .iat[minPos, 0], which is what the old code MEANT.
"""

from __future__ import annotations

import numpy as np

from anchovy.schema import (
    SIGNATURE_BARCODE_LEN,
    SIGNATURE_NON_UMI_LEN,
    SIGNATURE_PREFIX_LEN,
    SIGNATURE_SUFFIX_LEN,
)
import Levenshtein


def min_distance_block(blocks: np.ndarray, target: str) -> tuple[int, int, str]:
    """Find the block with the smallest Levenshtein distance to `target`.

    Consolidates the vectorized-distance pattern the original repeated inline in
    both blockDist and cellMatch.

    Args:
        blocks: array of candidate substrings.
        target: the string to match against.

    Returns:
        (min_distance, min_index, matching_block)
    """
    distance = lambda block: Levenshtein.distance(block, target)
    vector_distance = np.vectorize(distance, otypes=[int])
    dist = vector_distance(blocks).astype(int)
    min_pos = int(dist.argmin())
    return int(dist[min_pos]), min_pos, blocks[min_pos]


def make_read_blocks(seq: str, query_len: int) -> np.ndarray:
    """Slice a read sequence into overlapping query-length blocks.

    Reproduces the original blockDist windowing:
        [seq[i : i+query_len] for i in range(max(1, len(seq) - query_len))]
    The max(1, ...) guard preserves the original's behavior on short sequences.
    """
    return np.array([
        seq[i:min(len(seq), i + query_len)]
        for i in range(max(1, len(seq) - query_len))
    ])


def best_query_match(seq: str, query: str) -> tuple[int, int, str]:
    """Find where `query` best matches within `seq` (the blockDist core).

    Returns (min_distance, min_position, matched_block). On an empty/degenerate
    sequence the vectorized search would raise; we preserve the original's
    no-hit fallback of (len(query), -1, "") instead of propagating the error.
    """
    query_len = len(query)
    blocks = make_read_blocks(seq, query_len)
    try:
        return min_distance_block(blocks, query)
    except Exception:
        # Original fallback: no hit -> distance == query length, pos -1, empty.
        return query_len, -1, ""


def validate_signature(query: str) -> None:
    """Check a 10X signature has the layout the rest of this module assumes.

    WHY THIS EXISTS. The slice points below are not arbitrary: they encode a
    fixed 22-base 5' handle, a 16-base barcode and a fixed 10-base 3' handle,
    with the UMI width derived as len(query) - 48. That derivation is what lets
    a single signature string carry the whole chemistry -- swap the 26-N v2
    signature for the 28-N v3 one and the UMI follows automatically, because
    only the UMI length differs between them.

    It holds only while the signature actually has that shape. A signature with,
    say, a 4-base prefix is not rejected by anything: it just shifts the slice
    points, and the "UMI" that comes out is mostly barcode. Nothing downstream
    can tell, because a wrong-but-consistent UMI still groups reads -- it just
    groups the wrong ones, collapsing distinct molecules or splitting one.

    So the assumption is checked once, loudly, rather than silently relied on.

    Raises:
        ValueError: if the signature does not have the expected structure.
    """
    query = query.upper()
    n_start = query.find("N")
    n_end = query.rfind("N")

    if n_start == -1:
        raise ValueError(
            f"10X signature has no N-run marking the barcode+UMI region: "
            f"{query!r}")

    prefix, suffix = query[:n_start], query[n_end + 1:]
    n_run = query[n_start:n_end + 1]

    # The span between the first and last N must be ALL Ns. Checking the prefix
    # and suffix instead would be vacuous: they are defined as the text outside
    # those two positions, so neither can contain an N by construction.
    if set(n_run) != {"N"}:
        raise ValueError(
            f"10X signature must be <constant prefix><N-run><constant suffix>, "
            f"but the Ns are not contiguous -- found "
            f"{sorted(set(n_run) - {'N'})} inside the N-run: {query!r}")

    problems = []
    if len(prefix) != SIGNATURE_PREFIX_LEN:
        problems.append(
            f"5' handle is {len(prefix)} nt, expected {SIGNATURE_PREFIX_LEN}")
    if len(suffix) != SIGNATURE_SUFFIX_LEN:
        problems.append(
            f"3' handle is {len(suffix)} nt, expected {SIGNATURE_SUFFIX_LEN}")
    if len(n_run) <= SIGNATURE_BARCODE_LEN:
        problems.append(
            f"N-run is {len(n_run)} nt, which leaves no UMI after the "
            f"{SIGNATURE_BARCODE_LEN} nt barcode")

    if problems:
        umi = len(query) - SIGNATURE_NON_UMI_LEN
        raise ValueError(
            "10X signature does not have the layout anchovy assumes:\n"
            "    " + "\n    ".join(problems) + "\n"
            f"  got:      {query!r} ({len(query)} nt)\n"
            f"  expected: {SIGNATURE_PREFIX_LEN} nt constant 5' handle, then "
            f"{SIGNATURE_BARCODE_LEN} nt barcode + UMI as Ns, then "
            f"{SIGNATURE_SUFFIX_LEN} nt constant 3' handle.\n"
            f"  For reference, 10X v2 is 26 Ns (10 nt UMI) and v3 is 28 "
            f"(12 nt UMI).\n"
            f"  Proceeding would slice the UMI at the wrong offsets and "
            f"silently mis-group reads (this signature would give a "
            f"{umi} nt UMI).")


def build_barcode_query_blocks(query: str, barcodes) -> np.ndarray:
    """Build one barcode-augmented query template per whitelist barcode.

    Reproduces cellIDPool's construction:
        query[0:22] + barcode + "N"*(len(query)-48) + query[len-10:len]
    i.e. the fixed 5' handle, then the barcode, then N-padding sized to absorb
    varying UMI lengths, then the fixed 3' handle. The magic slice points
    (22, 48, 10) are structural to the 10X signature layout.

    Args:
        query: the 10X signature string.
        barcodes: iterable of barcode strings (e.g. pdCBCs.CBC).
    """
    n = len(query)
    return np.array([
        query[0:SIGNATURE_PREFIX_LEN]
        + bc
        + "N" * (n - SIGNATURE_NON_UMI_LEN)          # UMI width, derived
        + query[n - SIGNATURE_SUFFIX_LEN:n]
        for bc in barcodes
    ])


def extract_umi(matchseq: str, umi_start: int, umi_end_trim: int) -> str:
    """Slice the UMI out of a matched signature block.

    Reproduces matchseq[38:(len(matchseq)-10)] with the offsets parameterized
    (config.extract.umi_start_offset / umi_end_trim) instead of hardcoded.
    """
    return matchseq[umi_start:len(matchseq) - umi_end_trim]


def assign_barcode(matchseq: str, barcode_blocks: np.ndarray, whitelist):
    """Assign a read's matched signature to its closest whitelist barcode.

    The pure core of the original cellMatch: given the read's matched signature
    block and the precomputed barcode-augmented blocks, find the nearest barcode.

    Args:
        matchseq: the read's matched signature block.
        barcode_blocks: output of build_barcode_query_blocks.
        whitelist: the barcode DataFrame (single 'CBC' column).

    Returns:
        (barcode_string, min_distance, min_position, matched_block)

    PANDAS COMPAT: the original used whitelist.iloc[minPos][0], which on pandas
    2.x tries label-based lookup of column 0 and raises KeyError. We use
    positional .iat[minPos, 0] -- the behavior the original intended.
    """
    min_d, min_pos, matchblock = min_distance_block(barcode_blocks, matchseq)
    barcode = whitelist.iat[min_pos, 0]
    return barcode, min_d, min_pos, matchblock
