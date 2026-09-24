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
from anchovy.schema import (SamColumns, AnchovyColumns,
                            SIGNATURE_NON_UMI_LEN,
                            SIGNATURE_PREFIX_LEN, SIGNATURE_UMI_START)


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


def _assign_init(whitelist, barcode_blocks, barcode_index,
                 umi_start, umi_end_trim, max_barcode_errors=None) -> None:
    """Pool initializer: publish the shared read-only state to this worker."""
    _ASSIGN_STATE["whitelist"] = whitelist
    _ASSIGN_STATE["barcode_blocks"] = barcode_blocks
    _ASSIGN_STATE["barcode_index"] = barcode_index
    _ASSIGN_STATE["umi_start"] = umi_start
    _ASSIGN_STATE["umi_end_trim"] = umi_end_trim
    _ASSIGN_STATE["max_barcode_errors"] = max_barcode_errors


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
    barcode_index = _ASSIGN_STATE["barcode_index"]
    umi_start = _ASSIGN_STATE["umi_start"]
    umi_end_trim = _ASSIGN_STATE["umi_end_trim"]

    read_seq = full_seq[offset:read_len]
    assigned = barcodes.assign_barcode(
        matchseq, barcode_blocks, whitelist, barcode_index=barcode_index,
        max_barcode_errors=_ASSIGN_STATE["max_barcode_errors"],
    )
    if assigned is None:
        return None          # barcode beyond the limit; dropped by the caller
    barcode, min_d, min_pos, matchblock = assigned
    umi = barcodes.extract_umi(matchseq, umi_start, umi_end_trim)

    # Order must match AnchovyColumns.ORDER:
    # CBC, minD, BC_ID, readPos, clipLength, reconQuery, matchseq, read, mappedSeq, UMI
    return (barcode, min_d, min_pos, match_pos, offset,
            matchblock, matchseq, read_id, read_seq, umi)


# --------------------------------------------------------------------------- #
# Pipeline passes
# --------------------------------------------------------------------------- #
def _chunksize(n_items: int, n_workers: int) -> int:
    """How many items a worker takes at a time.

    Pool.map's own rule is ceil(n / (4 * workers)), which on a multi-million
    read run hands each worker a chunk of hundreds of thousands of payloads --
    every one of which carries a full read sequence, and all of which are
    pickled before the first result comes back. Capping the chunk bounds what
    the queue holds at any moment; the floor keeps the per-chunk overhead off
    the critical path on small inputs, where the whole thing fits anyway.
    """
    return max(1, min(1000, -(-n_items // max(1, n_workers * 4))))


def find_signature_positions(sam: pd.DataFrame, query: str,
                             config: ExtractConfig) -> pd.DataFrame:
    """Pass 1: locate the best signature match per read. (was: poolBlocks)

    Searches a window ending at each read's mapped offset. Adds minD, minPos,
    matchseq columns to the frame.
    """
    query_upper = query.upper()
    window = config.upstream_window

    # A GENERATOR, NOT A LIST, and imap rather than map. Pool.map materializes
    # the whole input, then pickles all of it into the task queue before any
    # worker starts, so the parent holds the windows twice over -- once as
    # objects and once as bytes. imap over a generator cuts the windows as the
    # queue drains, and the parent never holds more than a few chunks.
    search_inputs = (
        (row.seq[max(0, row.offset - window):row.offset].upper(), query_upper)
        for row in sam.itertuples()
    )

    with Pool(config.nthreads) as pool:
        print("\n1. Computing minimum distance hit position for {} reads."
              .format(len(sam)))
        results = list(pool.imap(_match_worker, search_inputs,
                                 chunksize=_chunksize(len(sam), config.nthreads)))

    # minD/minPos/matchseq are intermediate columns (not part of the SAM schema),
    # named literally to match what pass 2 reads by attribute.
    #
    # ASSIGNED IN PLACE. This used to copy the frame first, so that the caller's
    # frame was left untouched -- but the only caller reassigns its variable to
    # the return value, so the original became garbage the moment this returned
    # and the copy bought nothing. It cost a full duplicate of the frame,
    # sequences and all, live at the same time as the original: on a PacBio run
    # that was the single largest allocation in the stage.
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

    # Exact-match index over the whitelist. Every block differs from every
    # other in only the 16 barcode characters, so a read whose barcode region
    # is a whitelist entry can be resolved by lookup instead of by scanning all
    # of them. See barcodes.assign_barcode for why that is the same answer.
    barcode_index = barcodes.build_barcode_index(whitelist.CBC)

    kept = sam[sam.minD < config.min_distance_cutoff]

    payloads = [
        (row.read, row.seq, row.minPos, row.matchseq, row.readLen, row.offset)
        for row in kept.itertuples()
    ]

    # Report the split, because it is what decides this stage's runtime: a
    # lookup is ~46,000x cheaper than the scan, so the reads that miss are
    # essentially the whole cost.
    lo, hi = SIGNATURE_PREFIX_LEN, SIGNATURE_UMI_START
    exact = sum(1 for p in payloads if p[3][lo:hi] in barcode_index)
    total = len(payloads)
    if total:
        print("\n2. Identifying Cell Barcodes...")
        rest = ("the rest are resolved within {} error(s) or dropped"
                .format(config.max_barcode_errors)
                if config.max_barcode_errors is not None
                else "the rest fall back to the full scan")
        print("   {:,}/{:,} reads ({:.1%}) match a whitelist barcode exactly; {}."
              .format(exact, total, exact / total, rest))

    with Pool(config.nthreads, initializer=_assign_init,
              initargs=(whitelist, barcode_blocks, barcode_index,
                        config.umi_start_offset, config.umi_end_trim,
                        config.max_barcode_errors)) as pool:
        # imap for the same reason as pass 1, and it matters more here: a
        # payload carries the read's full sequence, so map's up-front pickle of
        # every one of them is a second copy of every base in the run.
        results = list(pool.imap(_assign_worker, payloads,
                                 chunksize=_chunksize(len(payloads),
                                                      config.nthreads)))

    # Reads whose barcode exceeded the limit come back as None. They were
    # never scanned -- that is the point of bounding the search rather than
    # filtering after it.
    if config.max_barcode_errors is not None:
        kept_results = [r for r in results if r is not None]
        dropped = len(results) - len(kept_results)
        print("   dropped {:,} read(s) ({:.1%}) whose barcode carried more than "
              "{} error(s); {:,} remain."
              .format(dropped, dropped / max(len(results), 1),
                      config.max_barcode_errors, len(kept_results)))
        results = kept_results

    out = pd.DataFrame(results, columns=AnchovyColumns.ORDER)

    # A perfect barcode still scores the UMI width, because every block is
    # padded with one N per UMI base and an N never matches a real base. So
    # that width is the floor, and anything above it is error in the barcode.
    _report_barcode_distances(out[AnchovyColumns.MIN_DISTANCE],
                              len(query) - SIGNATURE_NON_UMI_LEN)
    return out


def _report_barcode_distances(distances, baseline: int) -> None:
    """Print how well the assigned barcodes actually matched.

    Worth printing unprompted: the search returns the NEAREST whitelist entry
    whatever the distance, so without this a read whose barcode was unreadable
    is indistinguishable downstream from one that matched perfectly. The shape
    of this distribution is what says whether the cell assignments mean
    anything.
    """
    if len(distances) == 0:
        return
    exact = int((distances == baseline).sum())
    one = int((distances == baseline + 1).sum())
    two = int((distances == baseline + 2).sum())
    worse = int((distances > baseline + 2).sum())
    total = len(distances)

    print("   barcode match quality (a perfect barcode scores {}):"
          .format(baseline))
    for label, count in (("exact", exact), ("1 error", one),
                         ("2 errors", two), ("3+ errors", worse)):
        print("     {:<10} {:>9,} ({:5.1%})".format(label, count, count / total))


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

    # SAY WHAT THIS IS ABOUT TO COST, BEFORE SPENDING IT.
    #
    # The stage holds the surviving reads and forks a pool over them, so its
    # footprint is set here -- and when it exceeds what the machine has, the
    # kernel sends SIGKILL, which cannot be caught. The run then dies with
    # "died with <Signals.SIGKILL: 9>" and nothing else: no traceback, no stage
    # name, nothing to distinguish it from a crash. A colleague's run was
    # debugged from the snakemake log alone for want of this line.
    #
    # The constant is measured, not derived: 4.6 KB per surviving read per
    # worker, on a PacBio-shaped SAM, across the whole process tree.
    projected = 0.3 + config.nthreads * 4.6e-6 * len(sam_df)
    print("\n{:,} reads x {} worker(s) -- this stage needs roughly {:.1f} GB. "
          "If it is killed with no message, that is the kernel: lower "
          "extract_threads."
          .format(len(sam_df), config.nthreads, projected))

    sam_df = find_signature_positions(sam_df, query, config)
    print("\nMapped hits in {} reads."
          .format(int(np.sum([int(i) >= 0 for i in sam_df.minPos]))))

    sam_df = sam_df[sam_df.matchseq != ""]
    return assign_barcodes(sam_df, query, wl_df, config)
