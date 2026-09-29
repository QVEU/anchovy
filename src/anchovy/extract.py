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

from collections import Counter
from multiprocessing import Pool

import numpy as np
import pandas as pd

from anchovy import barcodes
from anchovy.config import ExtractConfig
from anchovy.io import AnchovyCsvSink, iter_sam_chunks, read_whitelist
from anchovy.schema import (SamColumns, AnchovyColumns,
                            SIGNATURE_NON_UMI_LEN,
                            SIGNATURE_PREFIX_LEN, SIGNATURE_UMI_START)


# --------------------------------------------------------------------------- #
# What the stage will cost
# --------------------------------------------------------------------------- #
# ONE DEFINITION, BECAUSE TWO DRIFTED. The workflow sizes its SLURM reservation
# from this and the stage prints it on the way in; when they were separate
# expressions the printed one kept a formula the reservation had already
# outgrown, and the run that most needed the warning was told 7.7 GB where 54
# was right.
def whitelist_worker_gb(n_barcodes: int, bounded: bool) -> float:
    """Per-worker cost of the whitelist lookup tables, in GB.

    A line through two measured points per mode -- 737,280 and 6,794,880
    barcodes, the real v2 and v3 lists:

        bounded    0.19 GB and 1.12 GB
        unbounded  0.36 GB and 2.75 GB

    PAID ONCE PER WORKER, not once. The tables are read-only and identical in
    every worker, but fork's copy-on-write does not save them: CPython's
    refcounts touch every object, so each worker ends up holding its own copy.

    Bounded is cheaper because the per-barcode template array is only built
    when the exhaustive scan will read it -- see _BarcodeState. The ~0.07 GB
    intercept is what a worker pays whatever the list holds.

    Args:
        n_barcodes: how many barcodes the whitelist actually holds.
        bounded: whether max_barcode_errors is set.
    """
    if bounded:
        return 0.077 + 1.535e-7 * n_barcodes
    return 0.069 + 3.945e-7 * n_barcodes


def projected_memory_gb(n_barcodes: int, threads: int, chunk_size: int,
                        bounded: bool) -> float:
    """What extract will hold at once, in GB.

    The tables once in the parent and once per worker, plus the chunk each
    worker is holding. Bounded by the chunk, NOT by the size of the run: this
    is the same number whether the SAM carries forty thousand reads or eleven
    million.

    On a full whitelist the table term dwarfs the chunk term -- 2.75 GB against
    0.46 GB at a 100,000-read chunk -- so `threads` is effectively the only
    memory knob that matters, and lowering chunk_size buys much less than it
    does on a small run-specific list.
    """
    w = whitelist_worker_gb(n_barcodes, bounded)
    return w + threads * (w + 4.6e-6 * chunk_size)


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
                 umi_start, umi_end_trim, max_barcode_errors=None,
                 query=None) -> None:
    """Pool initializer: publish the shared read-only state to this worker."""
    _ASSIGN_STATE["whitelist"] = whitelist
    # None when max_barcode_errors is set: see _BarcodeState. `query` then
    # stands in for it, one rebuilt template at a time.
    _ASSIGN_STATE["barcode_blocks"] = barcode_blocks
    _ASSIGN_STATE["barcode_index"] = barcode_index
    _ASSIGN_STATE["umi_start"] = umi_start
    _ASSIGN_STATE["umi_end_trim"] = umi_end_trim
    _ASSIGN_STATE["max_barcode_errors"] = max_barcode_errors
    _ASSIGN_STATE["query"] = query


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
        query=_ASSIGN_STATE["query"],
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
                             config: ExtractConfig, pool=None) -> pd.DataFrame:
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

    # `pool` is passed in by run(), which keeps one alive across every chunk:
    # forking a pool per chunk would pay the fork cost a hundred times over on a
    # large run. Left None, one is made here for the single-frame callers.
    owned = pool is None
    pool = pool or Pool(config.nthreads)
    try:
        results = list(pool.imap(_match_worker, search_inputs,
                                 chunksize=_chunksize(len(sam), config.nthreads)))
    finally:
        if owned:
            pool.terminate()

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


class _BarcodeState:
    """The pass-2 lookup tables, built once and reused for every chunk.

    WHY THIS IS A CLASS. build_barcode_query_blocks makes one search template
    per whitelist entry -- 6.8 million of them for a v3 list -- and the pool
    publishes them to its workers through an initializer. Both used to happen
    inside assign_barcodes, which was fine while the stage ran once over the
    whole file and ruinous the moment it runs once per chunk: a hundred chunks
    would rebuild and re-ship the whole thing a hundred times. Built here, the
    cost is paid once however the reads arrive.
    """

    def __init__(self, query, whitelist, config):
        self.whitelist = whitelist
        self.config = config
        self.query = query
        # THE BIG ARRAY IS ONLY BUILT WHEN SOMETHING WILL SCAN IT.
        #
        # One 60-character template per whitelist entry, numpy UCS-4, rebuilt in
        # every worker because fork's copy-on-write does not survive CPython's
        # refcounting. On a v3 list that is 1.6 GB per worker -- 59% of what a
        # worker holds -- and with max_barcode_errors set nothing ever reads it:
        # an index hit and a bounded search each resolve to one index, and
        # barcodes.barcode_query_block rebuilds that single template from the
        # barcode. Unset, the exhaustive scan needs all of them and this is the
        # cost of asking for it.
        self.blocks = (None if config.max_barcode_errors is not None
                       else barcodes.build_barcode_query_blocks(query,
                                                                whitelist.CBC))
        # Exact-match index over the whitelist. Every block differs from every
        # other in only the 16 barcode characters, so a read whose barcode
        # region is a whitelist entry can be resolved by lookup instead of by
        # scanning all of them. See barcodes.assign_barcode for why that is the
        # same answer.
        self.index = barcodes.build_barcode_index(whitelist.CBC)

    def pool(self):
        return Pool(self.config.nthreads, initializer=_assign_init,
                    initargs=(self.whitelist, self.blocks, self.index,
                              self.config.umi_start_offset,
                              self.config.umi_end_trim,
                              self.config.max_barcode_errors,
                              self.query))


class _Tally:
    """Counts accumulated across chunks, reported once at the end.

    Chunking would otherwise turn one summary into one per chunk, which on a
    hundred chunks is noise rather than reporting. Everything here is a scalar
    or a small Counter, so carrying it across the whole run costs nothing --
    the distances in particular are counted, not collected, which is the
    difference between a few hundred bytes and one integer per read.
    """

    def __init__(self):
        self.mapped_hits = 0
        self.payloads = 0
        self.exact = 0
        self.dropped = 0
        self.assigned = 0
        self.distances = Counter()

    def report(self, config, baseline):
        if self.payloads:
            rest = ("the rest are resolved within {} error(s) or dropped"
                    .format(config.max_barcode_errors)
                    if config.max_barcode_errors is not None
                    else "the rest fall back to the full scan")
            print("\n2. Identifying Cell Barcodes...")
            print("   {:,}/{:,} reads ({:.1%}) match a whitelist barcode "
                  "exactly; {}.".format(self.exact, self.payloads,
                                        self.exact / self.payloads, rest))
        if config.max_barcode_errors is not None:
            total = self.assigned + self.dropped
            print("   dropped {:,} read(s) ({:.1%}) whose barcode carried more "
                  "than {} error(s); {:,} remain."
                  .format(self.dropped, self.dropped / max(total, 1),
                          config.max_barcode_errors, self.assigned))
        _report_barcode_distances(self.distances, baseline)


def _assign_chunk(sam, state, pool, tally):
    """Pass 2 over one chunk. Returns the chunk's output rows."""
    config = state.config
    kept = sam[sam.minD < config.min_distance_cutoff]

    payloads = [
        (row.read, row.seq, row.minPos, row.matchseq, row.readLen, row.offset)
        for row in kept.itertuples()
    ]

    # Count the split, because it is what decides this stage's runtime: a
    # lookup is ~46,000x cheaper than the scan, so the reads that miss are
    # essentially the whole cost.
    lo, hi = SIGNATURE_PREFIX_LEN, SIGNATURE_UMI_START
    tally.exact += sum(1 for p in payloads if p[3][lo:hi] in state.index)
    tally.payloads += len(payloads)

    # imap over map for the same reason as pass 1, and it matters more here: a
    # payload carries the read's full sequence, so map's up-front pickle of
    # every one of them is a second copy of every base in the chunk.
    results = list(pool.imap(_assign_worker, payloads,
                             chunksize=_chunksize(len(payloads),
                                                  config.nthreads)))

    # Reads whose barcode exceeded the limit come back as None. They were
    # never scanned -- that is the point of bounding the search rather than
    # filtering after it.
    if config.max_barcode_errors is not None:
        kept_results = [r for r in results if r is not None]
        tally.dropped += len(results) - len(kept_results)
        results = kept_results

    tally.assigned += len(results)
    tally.distances.update(r[1] for r in results)
    return results


def assign_barcodes(sam: pd.DataFrame, query: str, whitelist: pd.DataFrame,
                    config: ExtractConfig) -> pd.DataFrame:
    """Pass 2 over a whole frame at once. (was: cellIDPool)

    The single-chunk path, kept for callers that already hold the entire frame.
    run() drives _assign_chunk directly so it can share one pool and one
    _BarcodeState across every chunk.
    """
    state = _BarcodeState(query, whitelist, config)
    tally = _Tally()
    with state.pool() as pool:
        results = _assign_chunk(sam, state, pool, tally)
    # A perfect barcode still scores the UMI width, because every block is
    # padded with one N per UMI base and an N never matches a real base. So
    # that width is the floor, and anything above it is error in the barcode.
    tally.report(config, len(query) - SIGNATURE_NON_UMI_LEN)
    return pd.DataFrame(results, columns=AnchovyColumns.ORDER)


def _report_barcode_distances(distances, baseline: int) -> None:
    """Print how well the assigned barcodes actually matched.

    Worth printing unprompted: the search returns the NEAREST whitelist entry
    whatever the distance, so without this a read whose barcode was unreadable
    is indistinguishable downstream from one that matched perfectly. The shape
    of this distribution is what says whether the cell assignments mean
    anything.
    """
    # A COUNTER, NOT A SERIES. Chunking means these arrive a chunk at a time and
    # the histogram is the only thing that needs all of them -- counted, that is
    # a handful of keys however many reads there were; collected, it was one
    # integer per read held until the very end.
    total = sum(distances.values())
    if total == 0:
        return
    exact = distances.get(baseline, 0)
    one = distances.get(baseline + 1, 0)
    two = distances.get(baseline + 2, 0)
    worse = sum(n for d, n in distances.items() if d > baseline + 2)

    print("   barcode match quality (a perfect barcode scores {}):"
          .format(baseline))
    for label, count in (("exact", exact), ("1 error", one),
                         ("2 errors", two), ("3+ errors", worse)):
        print("     {:<10} {:>9,} ({:5.1%})".format(label, count, count / total))


# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #
def run(sam: str, whitelist: str, signature: str | None = None,
        config: ExtractConfig | None = None, sink=None):
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
    wl_df = read_whitelist(whitelist)
    state = _BarcodeState(query, wl_df, config)
    tally = _Tally()

    # MEMORY IS BOUNDED BY THE CHUNK, NOT BY THE FILE.
    #
    # This used to read the whole SAM into one frame, run both passes over it,
    # and build the entire output table before returning any of it -- so the
    # stage cost scaled with the run. At eleven million reads the sequences
    # alone are ~22 GB, and the kernel killed it with SIGKILL, which cannot be
    # caught: the run died reporting only "died with <Signals.SIGKILL: 9>".
    #
    # Chunked, each of those is bounded by chunk_size instead, and the lookup
    # tables and the worker pools are built once outside the loop rather than
    # once per chunk. The output is byte-identical either way: chunks are
    # processed in file order and concatenated in that order.
    emit = sink if sink is not None else []

    # ONE POOL FOR BOTH PASSES, AND THE REASON IS FORK SAFETY, not tidiness.
    # A live Pool runs three management threads IN THE PARENT
    # (_handle_workers, _handle_tasks, _handle_results), so creating a second
    # pool forks a multi-threaded process -- which CPython 3.12 warns about and
    # which can genuinely deadlock: the child inherits a lock whose holding
    # thread does not exist in it. Two pools here made CI emit 68 of those
    # warnings where it had one.
    #
    # Sharing is free. _match_worker reads none of the state _assign_init
    # publishes, so a pool built for pass 2 serves pass 1 unchanged, and the
    # single fork happens before any pool thread exists.
    with state.pool() as pool:
        for n_chunk, chunk in enumerate(
                iter_sam_chunks(sam, config.effective_min_read_length(),
                                config.chunk_size), start=1):
            chunk = chunk[chunk[SamColumns.TEMPLATE] != "*"]
            if chunk.empty:
                continue
            if n_chunk == 1:
                print("\nProcessing in chunks of up to {:,} reads, {} worker(s)"
                      " against {:,} barcodes -- roughly {:.0f} GB at a time. "
                      "A run killed with no message is the kernel: lower "
                      "extract_threads or extract_chunk_size."
                      .format(config.chunk_size, config.nthreads, len(wl_df),
                              projected_memory_gb(
                                  len(wl_df), config.nthreads,
                                  config.chunk_size,
                                  config.max_barcode_errors is not None)))

            chunk = find_signature_positions(chunk, query, config, pool=pool)
            tally.mapped_hits += int((chunk.minPos >= 0).sum())
            chunk = chunk[chunk.matchseq != ""]

            rows = _assign_chunk(chunk, state, pool, tally)
            frame = pd.DataFrame(rows, columns=AnchovyColumns.ORDER)
            if sink is not None:
                sink.write(frame)
            else:
                emit.append(frame)
            print("   chunk {}: {:,} reads in, {:,} assigned ({:,} so far)"
                  .format(n_chunk, len(chunk), len(frame), tally.assigned))

    print("\nMapped hits in {} reads.".format(tally.mapped_hits))
    # A perfect barcode still scores the UMI width, because every block is
    # padded with one N per UMI base and an N never matches a real base. So
    # that width is the floor, and anything above it is error in the barcode.
    tally.report(config, len(query) - SIGNATURE_NON_UMI_LEN)

    if sink is not None:
        return tally.assigned
    if not emit:
        return pd.DataFrame(columns=AnchovyColumns.ORDER)
    return pd.concat(emit, ignore_index=True)


def run_to_csv(sam: str, whitelist: str, out_path: str,
               signature: str | None = None,
               config: ExtractConfig | None = None) -> int:
    """Stream the extract stage straight to its CSV. Returns rows written.

    The difference from run() is only where the output goes. run() concatenates
    every chunk and hands back one frame, which means the whole table is in
    memory at the end however small each chunk was -- fine for a test, and the
    thing this stage was dying of on a real run. Here each chunk is appended as
    it is produced and never held.
    """
    with AnchovyCsvSink(out_path) as sink:
        return run(sam=sam, whitelist=whitelist, signature=signature,
                   config=config, sink=sink)
