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
        query[0:22] + bc + "N" * (n - 48) + query[n - 10:n]
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
