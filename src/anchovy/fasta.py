"""
fasta.py -- the fasta stage: anchovy CSV -> one FASTA per cell barcode.

Migrated from CBCtoFasta.py. The original did everything in one function that
looped over barcodes and wrote files as a side effect, reading columns by
positional index (i[9] = read ID, i[10] = mapped sequence). That made it both
untestable (you had to write files and read them back to check anything) and
fragile (positional access).

THE SPLIT
---------
We separate the WHAT from the WHERE:

  build_cell_fastas(df, config) -> dict[barcode, fasta_text]
      Pure. Groups reads by barcode, applies the min-reads filter, formats each
      surviving cell's reads as FASTA text. No I/O. Trivially unit-testable:
      pass a small DataFrame, assert on the returned dict.

  write_cell_fastas(fastas, out_dir) -> list[path]
      Thin. Takes the dict and writes one .fa file per barcode. The only part
      that touches disk.

  run(csv, out_dir, config) -> list[path]
      Convenience wrapper: read CSV, build, write. What the CLI calls.

This mirrors the extract.py pattern (pure logic in barcodes.py, orchestration in
extract.py) and makes the golden test simple: compare the returned dict against
frozen expected FASTA text, no file plumbing in the assertion.

BEHAVIOR PRESERVATION
---------------------
Reproduces CBCtoFasta exactly:
  - only barcodes with >= min_reads_per_cbc reads are emitted,
  - each record is ">{read_id}\n{mapped_seq}",
  - records within a cell are joined by "\n", files are named "{barcode}.fa"
    with the barcode stripped of whitespace.
Uses schema column names instead of positional indices; the min-reads threshold
comes from config instead of the literal 5.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from anchovy.config import FastaConfig
from anchovy.schema import AnchovyColumns


def build_cell_fastas(df: pd.DataFrame,
                      config: FastaConfig | None = None) -> dict[str, str]:
    """Group reads by barcode and format per-cell FASTA text (pure).

    Args:
        df: anchovy table with at least the CBC, read, and mappedSeq columns
            (see schema.AnchovyColumns).
        config: FastaConfig; defaults to FastaConfig() (min_reads_per_cbc=5).

    Returns:
        dict mapping each surviving barcode (stripped) to its FASTA text.

    Only barcodes supported by at least config.min_reads_per_cbc reads are
    included, matching the original `if len(...) >= 5` filter.
    """
    config = config or FastaConfig()

    cbc_col = AnchovyColumns.CBC
    read_col = AnchovyColumns.READ
    seq_col = AnchovyColumns.MAPPED_SEQ

    fastas: dict[str, str] = {}
    # np.unique to iterate barcodes in a deterministic (sorted) order, matching
    # the original's np.unique(anchout.CBC) traversal.
    for barcode in np.unique(df[cbc_col]):
        cell_reads = df[df[cbc_col] == barcode]

        if len(cell_reads) >= config.min_reads_per_cbc:
            records = [
                ">{}\n{}".format(row[read_col], row[seq_col])
                for _, row in cell_reads.iterrows()
            ]
            fastas[str(barcode).strip()] = "\n".join(records)

    return fastas


def write_cell_fastas(fastas: dict[str, str], out_dir: str | Path) -> list[Path]:
    """Write each barcode's FASTA text to {out_dir}/{barcode}.fa (thin I/O).

    Returns the list of paths written, sorted for determinism.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    written: list[Path] = []
    for barcode, text in fastas.items():
        path = out_dir / f"{barcode}.fa"
        path.write_text(text)
        written.append(path)
    return sorted(written)


def run(csv: str, out_dir: str | Path,
        config: FastaConfig | None = None) -> list[Path]:
    """Read an anchovy CSV, build per-cell FASTAs, and write them.

    Args:
        csv: path to an anchovy output CSV (from the extract stage).
        out_dir: directory to write the per-cell .fa files into.
        config: FastaConfig; defaults to FastaConfig().

    Returns:
        list of written FASTA paths.
    """
    df = pd.read_csv(csv)
    fastas = build_cell_fastas(df, config)
    return write_cell_fastas(fastas, out_dir)
