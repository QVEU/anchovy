"""
ocm.py -- split one channel's anchovy CSV into its On-Chip Multiplexing samples.

WHAT OCM IS, AND WHERE THE SAMPLE LIVES
---------------------------------------
10x's On-Chip Multiplexing (GEM-X Universal 3' v4 / 5' v3 "4-plex") loads up to
four samples into the four inlets of ONE channel. Each inlet gets gel beads
from a different set (OB1-OB4), and the set is written into the cell barcode
itself: Cell Ranger's chemistry definitions for every OCM chemistry
(SC5P-R2-OCM-v3, SC3Pv4-polyA-OCM, ...) declare a 2 nt "overhang" barcode at
offset 7 of the 16 nt gel-bead barcode, checked against an `overhang` list that
translates to OB1-OB4. So the sample of every read is readable from its
barcode alone; no tag reads and no Cell Ranger output are needed to split.

The overhang -> OB table ships inside Cell Ranger, not in its public source, so
it is not hard-coded here. table_from_cellranger() reads it back out of a
`cellranger multi` run's per_sample_outs, which is also the check that the
offset holds for the chemistry in hand: every sample should show its own small
set of dinucleotides, and no dinucleotide should appear in two samples.

WHY SPLIT AFTER EXTRACT
-----------------------
Extract is the expensive stage and the samples share one library, so it runs
once per channel. Everything after it -- cells, consensus, genotypes, networks,
allele frequencies -- is only meaningful within one sample, so the split comes
before the per-cell stages and each sample then runs as an ordinary one.
"""

from __future__ import annotations

import gzip
import os
from collections import Counter
from pathlib import Path

import pandas as pd

from anchovy.schema import AnchovyColumns, SIGNATURE_BARCODE_LEN

# Cell Ranger chemistry_defs.json, every *-OCM chemistry: kind "overhang",
# offset 7, length 2, read R1 -- inside the 16 nt gel-bead barcode at offset 0.
OVERHANG_OFFSET = 7
OVERHANG_LEN = 2

_BASES = set("ACGT")


def overhang(barcode: str) -> str:
    """The OCM overhang (sample) bases of a 16 nt cell barcode."""
    return barcode[OVERHANG_OFFSET:OVERHANG_OFFSET + OVERHANG_LEN]


def validate_mapping(mapping: dict[str, list[str] | str]) -> dict[str, set[str]]:
    """Check a {sample: overhangs} mapping and return it as sets (raises).

    Overhangs may be given as a list or as one comma-separated string. Each must
    be OVERHANG_LEN bases of ACGT, and no overhang may belong to two samples --
    a read carrying it could not be assigned, and silently giving it to
    whichever sample came first would mix two samples' cells.
    """
    if not mapping:
        raise ValueError("OCM mapping is empty: name at least one sample.")
    out: dict[str, set[str]] = {}
    owner: dict[str, str] = {}
    for sample, overhangs in mapping.items():
        if isinstance(overhangs, str):
            overhangs = [o for o in overhangs.split(",") if o.strip()]
        seqs = {str(o).strip().upper() for o in overhangs}
        if not seqs:
            raise ValueError(f"OCM sample {sample!r} has no overhangs.")
        for seq in seqs:
            if len(seq) != OVERHANG_LEN or set(seq) - _BASES:
                raise ValueError(
                    f"OCM sample {sample!r}: {seq!r} is not an overhang -- "
                    f"expected {OVERHANG_LEN} bases of A/C/G/T (barcode "
                    f"positions {OVERHANG_OFFSET + 1}-"
                    f"{OVERHANG_OFFSET + OVERHANG_LEN}).")
            if seq in owner:
                raise ValueError(
                    f"OCM overhang {seq} is given to both {owner[seq]!r} and "
                    f"{sample!r}; a read carrying it could not be assigned.")
            owner[seq] = sample
        out[sample] = seqs
    return out


def split_csv(csv: str, overhangs: set[str], out_path: str,
              chunksize: int = 200_000) -> tuple[int, int]:
    """Write the rows of `csv` whose barcode carries one of `overhangs`.

    Streams the channel CSV in chunks, so memory is bounded by the chunk and
    not by the run. Always writes a header, so a sample with no reads still
    gives the next stage a readable, empty table. Returns (kept, total).
    """
    target = Path(out_path)
    tmp = target.with_name(target.name + ".partial")
    kept = total = 0
    index = 0
    with tmp.open("w", newline="") as handle:
        wrote_header = False
        for chunk in pd.read_csv(csv, chunksize=chunksize, index_col=0):
            total += len(chunk)
            mask = chunk[AnchovyColumns.CBC].astype(str).str[
                OVERHANG_OFFSET:OVERHANG_OFFSET + OVERHANG_LEN].isin(overhangs)
            sub = chunk[mask][AnchovyColumns.ORDER]
            sub.index = range(index, index + len(sub))
            index += len(sub)
            kept += len(sub)
            if len(sub) or not wrote_header:
                sub.to_csv(handle, header=not wrote_header)
                wrote_header = True
        if not wrote_header:
            pd.DataFrame(columns=AnchovyColumns.ORDER).to_csv(handle)
    os.replace(tmp, target)
    return kept, total


def run_split(csv: str, overhangs: list[str] | str, out_path: str,
              sample: str = "sample") -> int:
    seqs = validate_mapping({sample: overhangs})[sample]
    kept, total = split_csv(csv, seqs, out_path)
    share = kept / total if total else 0.0
    print(f"{sample}: kept {kept:,} of {total:,} reads ({share:.1%}) whose "
          f"barcode carries overhang {', '.join(sorted(seqs))}.")
    return kept


# --------------------------------------------------------------------------- #
# Reading the overhang table back out of Cell Ranger's own sample calls
# --------------------------------------------------------------------------- #
def _barcodes_file(sample_dir: Path) -> Path:
    for rel in ("count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz",
                "count/sample_filtered_feature_bc_matrix/barcodes.tsv"):
        path = sample_dir / rel
        if path.exists():
            return path
    raise FileNotFoundError(
        f"{sample_dir}: no count/sample_filtered_feature_bc_matrix/"
        f"barcodes.tsv(.gz) -- is this a cellranger multi per_sample_outs "
        f"sample folder?")


def overhang_counts(per_sample_outs: str) -> dict[str, Counter]:
    """{sample: Counter(overhang -> called cells)} from a per_sample_outs dir."""
    root = Path(per_sample_outs)
    counts: dict[str, Counter] = {}
    for sample_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        path = _barcodes_file(sample_dir)
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt") as handle:
            barcodes = [line.strip().split("-")[0] for line in handle
                        if line.strip()]
        bad = next((b for b in barcodes if len(b) != SIGNATURE_BARCODE_LEN),
                   None)
        if bad is not None:
            raise ValueError(f"{path}: {bad!r} is not a 16 nt cell barcode.")
        counts[sample_dir.name] = Counter(overhang(b) for b in barcodes)
    if not counts:
        raise ValueError(f"{root}: no sample folders found.")
    return counts


def table_from_cellranger(per_sample_outs: str,
                          min_share: float = 0.01) -> dict[str, list[str]]:
    """Infer {sample: overhangs} from Cell Ranger's per-sample cell calls.

    An overhang belongs to a sample when it carries at least `min_share` of
    that sample's cells. Raises when one overhang qualifies for two samples,
    which would mean the offset is wrong for this chemistry, or the run was not
    OCM.
    """
    counts = overhang_counts(per_sample_outs)
    table = {}
    for sample, c in counts.items():
        n = sum(c.values())
        table[sample] = sorted(o for o, k in c.items() if n and k / n >= min_share)
    validate_mapping(table)
    return table


def format_table(per_sample_outs: str, min_share: float = 0.01) -> str:
    """The table as YAML for the workflow's `ocm` key, with the evidence."""
    counts = overhang_counts(per_sample_outs)
    table = table_from_cellranger(per_sample_outs, min_share)
    lines = ["# Overhang (barcode bases 8-9) share of each sample's called "
             "cells:"]
    for sample, c in counts.items():
        n = sum(c.values())
        top = ", ".join(f"{o} {k / n:.1%}" for o, k in c.most_common(6))
        lines.append(f"#   {sample} ({n:,} cells): {top}")
    lines.append("ocm:")
    lines.append("  <channel FASTQ name>:")
    for sample, seqs in table.items():
        lines.append(f"    {sample}: [{', '.join(seqs)}]")
    return "\n".join(lines)
