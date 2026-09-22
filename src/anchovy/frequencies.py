"""
frequencies.py -- allele frequencies over ALL mapped cells, with per-position
denominators.

WHY THIS IS A SEPARATE STAGE FROM THE GENOTYPE NETWORK
------------------------------------------------------
The network and the frequency table want opposite things from the same data.

A genotype network keys on genotype STRINGS, so it needs cells whose strings
mean the same thing: a cell covering half the coding sequence emits a short
token list that is indistinguishable from a fully covered cell which happens to
be clean in the missing half. Admitting cells of differing breadth makes
genotype identity partly an artifact of what got sequenced, so the network is
built from a strict, complete-coverage subset.

A frequency does not need that, and filtering for it actively hurts. What a
frequency needs is a denominator that is honest about who could have been
counted -- so it is taken PER POSITION. A cell contributes at the positions it
covers and is absent elsewhere. Nothing is discarded for being partial, and no
cell inflates a denominator at a position it cannot speak to.

Both tables therefore come out of the same run, and the population table carries
the subset's own counts alongside the full population's, so the subset can be
checked for bias rather than assumed representative.

TWO LEVELS, WHICH ARE NOT INTERCHANGEABLE
-----------------------------------------
Per-cell (CellAlleleColumns): reads carrying an allele over depth at that
position, within one cell. A within-cell quasispecies frequency.

Population (AlleleFrequencyColumns): cells calling an allele over cells with a
confident call there. One cell, one vote.

The read-summed columns in the population table are provided because they are
asked for, but they weight each cell by its sequencing depth, and reads within a
cell are RT/PCR copies of a few templates rather than independent observations.
A 200-read cell outvotes a 20-read cell tenfold in those columns and not at all
in the cell columns. For population claims, use the cell columns.
"""

from __future__ import annotations

import csv
import warnings
from collections import Counter
from pathlib import Path

from anchovy.consensus import CALLED_BASES, parse_consensus_fasta
from anchovy.schema import AlleleFrequencyColumns, CellAlleleColumns

# The bases a pileup column is tallied over. 'N' and 'gap' are read from the
# counts file and count toward DEPTH -- a read with a deletion is still a read
# at that position -- but are never reported as alleles.
PILEUP_BASES = ("A", "C", "G", "T")
PILEUP_ALL = PILEUP_BASES + ("N", "gap")


def read_counts(path: str | Path) -> dict[int, dict[str, int]]:
    """Read one cell's per-position base counts (sam2consensus --counts).

    Returns {1-based position: {A, C, G, T, N, gap}}. POSITIONS ABSENT FROM THE
    FILE HAD NO COVERAGE -- sam2consensus omits them rather than writing a
    reference-length block of zeroes per cell, so a missing key means depth 0
    and not missing data.
    """
    out: dict[int, dict[str, int]] = {}
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            out[int(row["position"])] = {b: int(row[b]) for b in PILEUP_ALL}
    return out


def counts_path_for(cbc_id: str, counts_dir: str | Path,
                    reference_name: str) -> Path | None:
    """Locate a cell's counts file, or None.

    sam2consensus names its outputs "{reference}__{prefix}[...]", where prefix
    is the per-cell SAM minus '.sam', and writes FASTA headers as
    "{prefix}|c{threshold}". So the cell's counts file is recoverable from the
    consensus record's id by dropping the threshold suffix.
    """
    prefix = cbc_id.split("|")[0]
    path = Path(counts_dir) / f"{reference_name}__{prefix}_counts.tsv"
    return path if path.exists() else None


def cell_allele_rows(cbc_id: str, counts: dict[int, dict[str, int]],
                     reference: str, min_alt_reads: int = 2,
                     min_alt_freq: float = 0.0) -> list[dict]:
    """Per-cell pileup rows for one cell's non-reference alleles.

    THE THRESHOLDS ARE NOT OPTIONAL IN PRACTICE. Every covered position carries
    some minority read count from sequencing and PCR error, so reporting every
    allele with at least one read yields a row per covered position per cell --
    hundreds of thousands of rows that are overwhelmingly noise. min_alt_reads
    of 2 drops singletons, which is the bulk of it.
    """
    rows: list[dict] = []
    for pos in sorted(counts):
        if pos > len(reference):
            continue                      # insertion columns have no ref base
        col = counts[pos]
        depth = sum(col[b] for b in PILEUP_ALL)
        if depth == 0:
            continue
        ref_base = reference[pos - 1].upper()
        for base in PILEUP_BASES:
            if base == ref_base:
                continue
            reads = col[base]
            freq = reads / depth
            if reads < min_alt_reads or freq < min_alt_freq:
                continue
            rows.append({
                CellAlleleColumns.CBC_ID: cbc_id,
                CellAlleleColumns.POSITION: pos,
                CellAlleleColumns.REF_BASE: ref_base,
                CellAlleleColumns.ALLELE: base,
                CellAlleleColumns.READS: reads,
                CellAlleleColumns.DEPTH: depth,
                CellAlleleColumns.FREQ: round(freq, 6),
            })
    return rows


def population_rows(records, reference: str, subset_ids: set[str] | None = None,
                    read_totals: dict[int, dict[str, int]] | None = None,
                    regions=None, min_cells: int = 1) -> list[dict]:
    """Per (position, allele) rows over every cell's consensus.

    A row exists for an allele called in at least `min_cells` cells. Read
    support is reported for those same alleles rather than for every allele with
    a read behind it, which keeps the table to observed mutations instead of the
    error spectrum.

    DENOMINATORS ARE PER POSITION AND EXCLUDE UNCALLED BASES. A cell whose
    consensus holds a gap (no coverage) or an IUPAC code (reads disagreed) at a
    position is counted in neither numerator nor denominator there: it has no
    confident call, so including it would understate every frequency at that
    position by the number of cells that simply did not resolve it.
    """
    subset_ids = subset_ids if subset_ids is not None else set()
    ncols = min(len(reference), min((len(r.seq) for r in records), default=0))
    rows: list[dict] = []

    for i in range(ncols):
        ref_base = reference[i].upper()
        all_counts: Counter = Counter()
        sub_counts: Counter = Counter()
        for r in records:
            base = r.seq[i].upper()
            if base not in CALLED_BASES:
                continue
            all_counts[base] += 1
            if r.cbc_id in subset_ids:
                sub_counts[base] += 1

        cells_called = sum(all_counts.values())
        cells_called_sub = sum(sub_counts.values())
        if cells_called == 0:
            continue

        for base in sorted(b for b in all_counts if b != ref_base):
            cells_alt = all_counts[base]
            if cells_alt < min_cells:
                continue
            col = (read_totals or {}).get(i + 1)
            reads_alt = col[base] if col else ""
            reads_depth = sum(col[b] for b in PILEUP_ALL) if col else ""
            rows.append({
                AlleleFrequencyColumns.POSITION: i + 1,
                AlleleFrequencyColumns.REGION: _region_names(regions, i + 1),
                AlleleFrequencyColumns.REF_BASE: ref_base,
                AlleleFrequencyColumns.ALLELE: base,
                AlleleFrequencyColumns.CELLS_ALT: cells_alt,
                AlleleFrequencyColumns.CELLS_CALLED: cells_called,
                AlleleFrequencyColumns.FREQ_CELLS: round(cells_alt / cells_called, 6),
                AlleleFrequencyColumns.CELLS_ALT_SUBSET: sub_counts[base],
                AlleleFrequencyColumns.CELLS_CALLED_SUBSET: cells_called_sub,
                AlleleFrequencyColumns.FREQ_CELLS_SUBSET: (
                    round(sub_counts[base] / cells_called_sub, 6)
                    if cells_called_sub else ""),
                AlleleFrequencyColumns.READS_ALT: reads_alt,
                AlleleFrequencyColumns.READS_DEPTH: reads_depth,
                AlleleFrequencyColumns.FREQ_READS: (
                    round(reads_alt / reads_depth, 6)
                    if reads_depth else ""),
            })
    return rows


def _region_names(regions, genome_pos: int) -> str:
    """Names of the GFF regions containing a 1-based position, ';'-joined."""
    if not regions:
        return ""
    from anchovy.regions import classify
    return ";".join(r.name for r in classify(regions, genome_pos))


def _write(path: Path, rows: list[dict], order: list[str]) -> Path:
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=order)
        writer.writeheader()
        writer.writerows(rows)
    return path


def run(fasta: str, reference: str, out_prefix: str,
        counts_dir: str | None = None, reference_name: str = "",
        subset_csv: str | None = None, gff: str | None = None,
        min_alt_reads: int = 2, min_alt_freq: float = 0.0,
        min_cells: int = 1) -> dict:
    """Full frequency stage: population table, and per-cell table when counts exist.

    `counts_dir` is optional on purpose. Without it the population table is
    still produced from the consensuses alone -- the cell columns need nothing
    else -- and only the read columns are left blank. That keeps the stage
    runnable against results produced before sam2consensus --counts existed.
    """
    records = parse_consensus_fasta(fasta)

    # A REFERENCE OF THE WRONG LENGTH IS THE SILENT FAILURE HERE. Positions are
    # compared by index, so a reference that is shifted or truncated relative to
    # the consensuses does not raise -- it reports nearly every position as
    # mutant at a frequency near 1.0, which reads as a spectacular biological
    # result rather than as mismatched inputs. Length is the cheap tell, and it
    # catches the common cause: passing a trimmed _consensus_reference.txt from
    # an earlier run instead of the genome the reads were mapped to.
    if records:
        modal = Counter(len(r.seq) for r in records).most_common(1)[0][0]
        if len(reference) != modal:
            warnings.warn(
                f"reference is {len(reference)} nt but the consensuses are "
                f"{modal} nt.\n"
                f"  Positions are compared by index, so a mismatched reference "
                f"reports nearly every\n"
                f"  position as mutant instead of failing. Pass the genome the "
                f"reads were mapped to,\n"
                f"  not a trimmed _consensus_reference.txt from an earlier run.",
                stacklevel=2)

    subset_ids: set[str] = set()
    if subset_csv:
        with open(subset_csv) as handle:
            subset_ids = {row["CBC_ID"] for row in csv.DictReader(handle)}

    regions = None
    if gff:
        from anchovy.regions import parse_gff3
        regions = parse_gff3(gff)

    # Per-cell rows, and the summed pileup the population table's read columns
    # use. Both come from one pass so the counts files are read only once.
    cell_rows: list[dict] = []
    read_totals: dict[int, dict[str, int]] = {}
    cells_with_counts = 0
    if counts_dir:
        for record in records:
            path = counts_path_for(record.cbc_id, counts_dir, reference_name)
            if path is None:
                continue
            cells_with_counts += 1
            counts = read_counts(path)
            cell_rows.extend(cell_allele_rows(
                record.cbc_id, counts, reference,
                min_alt_reads=min_alt_reads, min_alt_freq=min_alt_freq))
            for pos, col in counts.items():
                total = read_totals.setdefault(
                    pos, {b: 0 for b in PILEUP_ALL})
                for base in PILEUP_ALL:
                    total[base] += col[base]

    pop_rows = population_rows(
        records, reference, subset_ids=subset_ids,
        read_totals=read_totals or None, regions=regions, min_cells=min_cells)

    written = {"population": _write(
        Path(f"{out_prefix}_alleleFrequencies.csv"), pop_rows,
        AlleleFrequencyColumns.ORDER)}
    if counts_dir:
        written["per_cell"] = _write(
            Path(f"{out_prefix}_cellAlleleFreq.csv"), cell_rows,
            CellAlleleColumns.ORDER)

    stats = {
        "cells": len(records),
        "cells_in_subset": len(subset_ids),
        "cells_with_counts": cells_with_counts,
        "population_rows": len(pop_rows),
        "cell_rows": len(cell_rows),
    }
    return {"population": pop_rows, "per_cell": cell_rows,
            "written": written, "stats": stats}
