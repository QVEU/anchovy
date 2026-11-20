"""
cell_lists.py -- the run's matched-barcode list, and custom cell lists against it.

THE SPLIT THIS MODULE MAKES
---------------------------
Extract always matches reads against the FULL 10x whitelist for the chemistry.
Which of those barcodes are real cells is a separate question, answered by
Cell Ranger or Seurat, and the answer changes whenever the single-cell analysis
is re-run with different parameters. So it is kept out of extract entirely:

  1. barcode_table()   -- every barcode extract assigned reads to, with its
                          read and UMI counts (and, for an On-Chip
                          Multiplexing channel, the sample its overhang names).
                          This is the run's MASTER list: `<sample>_barcodes.csv`.
  2. match_cell_list() -- the rows of a master list whose barcode is in a custom
                          cell list. Cheap, so a new list from a new clustering
                          is a seconds-long `anchovy match-cells`, not a re-run
                          of extract.

WHY NOT MATCH AGAINST THE CELL LIST DIRECTLY
--------------------------------------------
That is what anchovy used to recommend, and it had three costs. A filtered list
covers only some of the library's reads, so the exact-match floor
(config.min_exact_rate) tripped on data that was fine -- most often one
sample's cells from a multiplexed channel. The answer was tied to one version
of the cell calls, so a new clustering meant re-running the most expensive
stage. And with mismatches allowed, a read from a cell NOT in the list was
pulled onto its nearest listed neighbour instead of being left out.
"""

from __future__ import annotations

import gzip
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from anchovy.schema import AnchovyColumns, SIGNATURE_BARCODE_LEN

# Columns of the master list.
CBC = "CBC"
READS = "reads"
UMIS = "umis"
OCM_SAMPLE = "ocm_sample"

_BARCODE = re.compile(r"[ACGT]{%d}" % SIGNATURE_BARCODE_LEN)


# --------------------------------------------------------------------------- #
# The master list
# --------------------------------------------------------------------------- #
def barcode_table(anchovy_csv: str, ocm: dict[str, set[str]] | None = None,
                  chunksize: int = 200_000) -> pd.DataFrame:
    """Every assigned barcode in an anchovy CSV, with reads and distinct UMIs.

    Streams the CSV, so memory is one chunk plus one entry per barcode. Sorted
    by reads, most first. With `ocm` ({sample: overhangs}), adds the sample
    each barcode's overhang belongs to (blank for an overhang no sample owns).
    """
    reads: Counter = Counter()
    umis: dict[str, set] = defaultdict(set)
    for chunk in pd.read_csv(anchovy_csv, chunksize=chunksize,
                             usecols=[AnchovyColumns.CBC, AnchovyColumns.UMI],
                             dtype=str, keep_default_na=False):
        reads.update(chunk[AnchovyColumns.CBC])
        for cbc, umi in zip(chunk[AnchovyColumns.CBC],
                            chunk[AnchovyColumns.UMI]):
            umis[cbc].add(umi)

    table = pd.DataFrame({
        CBC: list(reads),
        READS: list(reads.values()),
        UMIS: [len(umis[c]) for c in reads],
    }, columns=[CBC, READS, UMIS])
    if ocm:
        from anchovy.ocm import overhang
        owner = {seq: sample for sample, seqs in ocm.items() for seq in seqs}
        table[OCM_SAMPLE] = [owner.get(overhang(c), "") for c in table[CBC]]
    return (table.sort_values([READS, CBC], ascending=[False, True])
            .reset_index(drop=True))


def _write_atomic(frame: pd.DataFrame, out_path: str) -> None:
    target = Path(out_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_name(target.name + ".partial")
    try:
        frame.to_csv(tmp, index=False)
        os.replace(tmp, target)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def write_barcode_table(anchovy_csv: str, out_path: str,
                        ocm: dict[str, set[str]] | None = None) -> pd.DataFrame:
    table = barcode_table(anchovy_csv, ocm=ocm)
    _write_atomic(table, out_path)
    total = int(table[READS].sum())
    print(f"{len(table):,} barcodes matched from {total:,} reads -> {out_path}")
    if ocm:
        for sample, n in table.groupby(OCM_SAMPLE)[READS].sum().items():
            print(f"   {sample or '(no sample)'}: {n:,} reads "
                  f"({n / max(total, 1):.1%})")
    return table


def read_barcode_table(path: str) -> pd.DataFrame:
    table = pd.read_csv(path, dtype={CBC: str}, keep_default_na=False)
    if CBC not in table.columns:
        raise ValueError(
            f"{path} has no {CBC!r} column -- expected a `_barcodes.csv` master "
            f"list written by `anchovy barcodes` or the workflow.")
    return table


# --------------------------------------------------------------------------- #
# Custom cell lists
# --------------------------------------------------------------------------- #
def read_cell_list(path: str) -> set[str]:
    """The 16 nt barcodes in a cell list, however it was exported.

    Accepts what Cell Ranger and Seurat actually write, so nobody has to clean
    a list by hand first: barcodes.tsv(.gz) with the `-1` GEM suffix, a CSV
    with a header, Seurat column names with a sample prefix (`P10_5_AAAC...-1`)
    or a `_1` suffix. Each line contributes the first 16-base A/C/G/T run in
    its first field. A line with none -- a header -- is skipped; a list where
    most lines have none is refused, because that is the wrong file.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    found, skipped, lines = set(), [], 0
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            lines += 1
            field = re.split(r"[,\t ]", line, maxsplit=1)[0].strip('"\'')
            hit = _BARCODE.search(field.upper())
            if hit:
                found.add(hit.group(0))
            else:
                skipped.append(field)
    if not found or len(skipped) > max(1, lines // 10):
        example = f" (e.g. {skipped[0]!r})" if skipped else ""
        raise ValueError(
            f"{path}: {len(skipped):,} of {lines:,} line(s) carry no 16 nt "
            f"barcode{example}. Expected one cell barcode per line, as in "
            f"Cell Ranger's barcodes.tsv.gz.")
    return found


def match_cell_list(table: pd.DataFrame, cells: set[str]) -> pd.DataFrame:
    """Rows of a master list whose barcode is in `cells`."""
    return table[table[CBC].isin(cells)].reset_index(drop=True)


def run_match(barcodes_csv: str, cell_list: str, out_path: str,
              label: str | None = None) -> dict:
    """`anchovy match-cells`: write the master-list rows in a cell list.

    Prints how much of the list was found and what share of the run's reads it
    holds -- the two numbers that say whether the list belongs to this run.
    """
    table = read_barcode_table(barcodes_csv)
    cells = read_cell_list(cell_list)
    matched = match_cell_list(table, cells)
    _write_atomic(matched, out_path)

    total = int(table[READS].sum()) if READS in table else 0
    in_list = int(matched[READS].sum()) if READS in matched else 0
    stats = {"cells": len(cells), "found": len(matched),
             "reads": in_list, "total_reads": total}
    name = label or Path(cell_list).name
    print(f"{name}: {len(matched):,} of {len(cells):,} listed cells "
          f"({len(matched) / len(cells):.1%}) have reads in this run, holding "
          f"{in_list:,} of {total:,} reads ({in_list / max(total, 1):.1%}).")
    if OCM_SAMPLE in matched and len(matched):
        by_sample = matched.groupby(OCM_SAMPLE).size()
        print("   by OCM sample: " + ", ".join(
            f"{s or '(none)'} {n:,}" for s, n in by_sample.items()))
    if cells and len(matched) / len(cells) < 0.5:
        print(f"   NOTE: under half of {name} was found. A list from another "
              f"run, another sample, or the other (3'/5') kit looks like this.")
    print(f"Wrote {out_path}")
    return stats


def select_cells(anchovy_csv: str, cells: set[str], out_path: str,
                 chunksize: int = 200_000) -> tuple[int, int]:
    """Write the rows of an anchovy CSV whose barcode is in `cells`.

    The per-cell stages then see only the listed cells. Streams in chunks and
    always writes a header, so an empty selection still gives the next stage a
    readable table. Returns (kept, total).
    """
    target = Path(out_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_name(target.name + ".partial")
    kept = total = 0
    try:
        with tmp.open("w", newline="") as handle:
            wrote_header = False
            for chunk in pd.read_csv(anchovy_csv, chunksize=chunksize,
                                     index_col=0):
                total += len(chunk)
                sub = chunk[chunk[AnchovyColumns.CBC].astype(str).isin(cells)]
                sub = sub[AnchovyColumns.ORDER]
                sub.index = range(kept, kept + len(sub))
                kept += len(sub)
                if len(sub) or not wrote_header:
                    sub.to_csv(handle, header=not wrote_header)
                    wrote_header = True
            if not wrote_header:
                pd.DataFrame(columns=AnchovyColumns.ORDER).to_csv(handle)
        os.replace(tmp, target)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    return kept, total


def run_select(anchovy_csv: str, cell_list: str, out_path: str) -> int:
    """`anchovy select-cells`: an anchovy CSV restricted to a cell list."""
    cells = read_cell_list(cell_list)
    kept, total = select_cells(anchovy_csv, cells, out_path)
    print(f"{Path(cell_list).name}: kept {kept:,} of {total:,} reads "
          f"({kept / max(total, 1):.1%}) from the {len(cells):,} listed cells "
          f"-> {out_path}")
    return kept


def label_for(path: str) -> str:
    """A filename-safe label for a cell list, from its file name."""
    name = Path(path).name
    for ext in (".gz", ".tsv", ".csv", ".txt"):
        if name.endswith(ext):
            name = name[:-len(ext)]
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name) or "cells"
