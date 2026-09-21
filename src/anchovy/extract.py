"""
extract.py -- the anchovy extract stage: mapped SAM -> per-read CBC/UMI table.

This is the migrated core of the original anchovy.py. It orchestrates:
  1. read the SAM (io.read_sam) and whitelist (io.read_whitelist)
  2. find where the 10X signature best matches in each read (poolBlocks)
  3. assign each read to its closest whitelist barcode + extract UMI (cellIDPool)

The pure matching math lives in barcodes.py; the I/O lives in io.py; the tunable
parameters live in config.py. This module holds only the glue and the parallelism.

PUBLIC API
----------
run(sam, whitelist, signature=..., config=...) -> pd.DataFrame
    Returns the anchovy table in canonical schema column order. This is what the
    golden test calls. A thin CLI wrapper (cli.py) adds file output.

DESIGN NOTES
------------
- The multiprocessing worker functions must be module-level (not closures) so
  they're picklable by multiprocessing. They take a single tuple argument, which
  is why the arguments are packed/unpacked -- the same reason the original did it.
- We pass column VALUES into the workers by name (via itertuples attributes),
  not by positional index. This is the core fix for the original's f[1]/f[10]/
  f[17]... fragility: reordering columns can no longer silently corrupt results.
"""

from __future__ import annotations

from multiprocessing import Pool

import numpy as np
import pandas as pd

from anchovy import barcodes
from anchovy.config import ExtractConfig
from anchovy.io import read_sam, read_whitelist
from anchovy.schema import SamColumns, AnchovyColumns


# --------------------------------------------------------------------------- #
# Multiprocessing workers (module-level so they pickle cleanly)
# --------------------------------------------------------------------------- #
def _match_worker(seq_query: tuple[str, str]) -> tuple[int, int, str]:
    """Worker for the signature-search pass. (was: blockDist)"""
    seq, query = seq_query
    return barcodes.best_query_match(seq, query)


# The whitelist and its derived query blocks are READ-ONLY and identical for
# every read, so they are handed to each worker process once at startup instead
# of travelling inside every payload.
#
# They used to be two of the ten fields in each payload tuple. In-process that
# is just a reference, but Pool.map PICKLES the payloads, and the big objects
# are re-serialized once per chunk. At 10X scale that is not a micro-
# optimization: a 737,280-barcode whitelist makes a single payload 185 MB, and
# 402k reads over 16 processes is ~64 chunks -- roughly 12 GB serialized each
# way, and the parent buffers those chunks as it fills the task queue. The
# worker needed all of it to do one positional lookup, whitelist.iat[pos, 0].
#
# With the fork start method the children inherit initargs through the fork
# itself, so this costs nothing at all; under spawn it is paid once per worker
# rather than once per chunk.
_ASSIGN_STATE: dict = {}


def _assign_init(whitelist, barcode_blocks, umi_start, umi_end_trim) -> None:
    """Pool initializer: publish the shared read-only state to this worker."""
    _ASSIGN_STATE["whitelist"] = whitelist
    _ASSIGN_STATE["barcode_blocks"] = barcode_blocks
    _ASSIGN_STATE["umi_start"] = umi_start
    _ASSIGN_STATE["umi_end_trim"] = umi_end_trim


def _assign_worker(payload: tuple) -> tuple:
    """Worker for the barcode-assignment pass. (was: cellMatch)

    payload = (read_id, full_seq, match_pos, matchseq, read_len, offset)
    The whitelist, barcode blocks and UMI offsets come from _ASSIGN_STATE,
    which _assign_init populated once when this process started.

    Returns the 10 fields in schema.AnchovyColumns.ORDER order.
    """
    read_id, full_seq, match_pos, matchseq, read_len, offset = payload

    whitelist = _ASSIGN_STATE["whitelist"]
    barcode_blocks = _ASSIGN_STATE["barcode_blocks"]
    umi_start = _ASSIGN_STATE["umi_start"]
    umi_end_trim = _ASSIGN_STATE["umi_end_trim"]

    read_seq = full_seq[offset:read_len]
    barcode, min_d, min_pos, matchblock = barcodes.assign_barcode(
        matchseq, barcode_blocks, whitelist
    )
    umi = barcodes.extract_umi(matchseq, umi_start, umi_end_trim)

    # Order must match AnchovyColumns.ORDER:
    # CBC, minD, BC_ID, readPos, clipLength, reconQuery, matchseq, read, mappedSeq, UMI
    return (barcode, min_d, min_pos, match_pos, offset,
            matchblock, matchseq, read_id, read_seq, umi)


# --------------------------------------------------------------------------- #
# Pipeline passes
# --------------------------------------------------------------------------- #
def find_signature_positions(sam: pd.DataFrame, query: str,
                             config: ExtractConfig) -> pd.DataFrame:
    """Pass 1: locate the best signature match per read. (was: poolBlocks)

    Searches a window ending at each read's mapped offset. Adds minD, minPos,
    matchseq columns to the frame.
    """
    query_upper = query.upper()
    window = config.upstream_window

    search_inputs = [
        (row.seq[max(0, row.offset - window):row.offset].upper(), query_upper)
        for row in sam.itertuples()
    ]

    with Pool(config.nthreads) as pool:
        print("\n1. Computing minimum distance hit position for {} reads."
              .format(len(sam)))
        results = pool.map(_match_worker, search_inputs)

    # minD/minPos/matchseq are intermediate columns (not part of the SAM schema),
    # named literally to match what pass 2 reads by attribute.
    sam = sam.copy()
    sam["minD"], sam["minPos"], sam["matchseq"] = zip(*results)
    return sam


def assign_barcodes(sam: pd.DataFrame, query: str, whitelist: pd.DataFrame,
                    config: ExtractConfig) -> pd.DataFrame:
    """Pass 2: assign each read to its closest barcode + extract UMI.

    (was: cellIDPool) Filters to reads whose signature match passed the distance
    cutoff, builds the barcode-augmented query blocks once, then matches each
    read in parallel.
    """
    barcode_blocks = barcodes.build_barcode_query_blocks(query, whitelist.CBC)

    kept = sam[sam.minD < config.min_distance_cutoff]

    payloads = [
        (row.read, row.seq, row.minPos, row.matchseq, row.readLen, row.offset)
        for row in kept.itertuples()
    ]

    with Pool(config.nthreads, initializer=_assign_init,
              initargs=(whitelist, barcode_blocks,
                        config.umi_start_offset, config.umi_end_trim)) as pool:
        print("\n2. Identifying Cell Barcodes...")
        results = pool.map(_assign_worker, payloads)

    out = pd.DataFrame(results, columns=AnchovyColumns.ORDER)
    return out


# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #
def run(sam: str, whitelist: str, signature: str | None = None,
        config: ExtractConfig | None = None) -> pd.DataFrame:
    """Run the full extract stage and return the anchovy table.

    Args:
        sam: path to the mapped SAM/BAM.
        whitelist: path to the 10X barcode whitelist.
        signature: 10X signature to search for. Defaults to config.signature.
            Must have the standard layout (22 nt handle, N-run covering the
            16 nt barcode plus the UMI, 10 nt handle) -- validated up front.
            Changing chemistry is exactly this plus a matching whitelist: the
            UMI width is derived from the signature's length, so the v2 (26 N)
            and v3 (28 N) signatures need no other configuration.
        config: an ExtractConfig; defaults to ExtractConfig() (original values).

    Returns:
        DataFrame in schema.AnchovyColumns.ORDER order.

    Raises:
        ValueError: if the signature does not have the expected layout.
    """
    config = config or ExtractConfig()
    query = (signature or config.signature).upper()

    # Fail here rather than three stages later. The slice points that pull the
    # barcode and UMI out of a matched block assume a specific signature layout;
    # a signature that does not have it produces a wrong-but-plausible UMI that
    # nothing downstream can detect. See barcodes.validate_signature.
    barcodes.validate_signature(query)

    print("Query Length: {}".format(len(query)))
    sam_df = read_sam(sam, config.effective_min_read_length())
    sam_df = sam_df[sam_df[SamColumns.TEMPLATE] != "*"]
    wl_df = read_whitelist(whitelist)

    sam_df = find_signature_positions(sam_df, query, config)
    print("\nMapped hits in {} reads."
          .format(int(np.sum([int(i) >= 0 for i in sam_df.minPos]))))

    sam_df = sam_df[sam_df.matchseq != ""]
    return assign_barcodes(sam_df, query, wl_df, config)
