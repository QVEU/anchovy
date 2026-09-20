"""
consensus.py -- the consensus stage: filter aligned per-cell consensus sequences,
compute (or accept) a reference, and call per-cell genotypes.

Migrated from ConsensusTool.py. Two behavior-relevant decisions are baked in:

1. THE if/else UNIFICATION (a deliberate, chosen behavior change).
   The original genotypeSummary had two branches: an if-branch that computed its
   own consensus (the only path main() ever used) and an else-branch for a
   supplied reference that was DEAD and BROKEN -- it referenced an undefined
   `consensusSeq` and would raise NameError if ever reached.
   Per the design decision, we replace both with ONE path where the reference is
   a parameter: compute the consensus when none is given, otherwise use the one
   supplied. This removes the broken dead code and makes the reference path
   actually work.

   Consequence for testing:
     - reference=None MUST reproduce the original's working behavior exactly.
       This is golden-tested against the frozen (and hand-verified) output.
     - reference=<given> is NEW capability with no prior behavior to match.
       It is unit-tested on small inputs where the correct answer is known.

2. THE SPLIT (pure core + thin I/O), same as extract/fasta.
     parse_consensus_fasta(path)      -> records with seq + parsed coverage/length
     select_sequences(records, ...)   -> pure filter + trim
     genotype_summary(seqs, ref=None) -> (genotypes, consensus)  [PURE, the core]
     run(...)                          -> orchestration + file writing

BEHAVIOR PRESERVATION (reference=None path)
-------------------------------------------
Reproduces genotypeSummary's if-branch exactly:
  - variant sites = columns (of the trimmed alignment) with >1 unique character
  - per-column consensus = most frequent character (argmax of unique counts)
  - each sequence's genotype = "_".join((i+1)+base) over variant sites where the
    sequence differs from the per-column consensus, 1-based within the trimmed region
Filtering reproduces selectSeqs: keep depth > depth_min AND (< max_gaps '-' in the
trimmed region), then trim to [start:end].
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from anchovy.config import ConsensusConfig
from anchovy.io import parse_description
from anchovy.schema import ConsensusColumns


@dataclass
class ConsensusRecord:
    """One consensus sequence plus its parsed metadata.

    A small typed container so downstream code uses .cbc_id / .seq / .coverage
    instead of positional tuple access -- the same anti-fragility principle as
    schema.py, applied to in-memory records.
    """
    cbc_id: str
    seq: str
    description: str
    coverage: float | None
    length: float | None


def parse_consensus_fasta(path: str | Path) -> list[ConsensusRecord]:
    """Read an *_allConsensus.fasta into typed records.

    Uses io.parse_description (the shared coverage/length parser) instead of the
    original's positional description.split(" ")[2] approach, so a reordered or
    differently-spaced header no longer breaks parsing.

    Minimal FASTA reader (no biopython dependency needed for this simple format):
    a header line starting with '>' followed by one sequence line.
    """
    records: list[ConsensusRecord] = []
    header: str | None = None
    seq_parts: list[str] = []

    def flush():
        if header is not None:
            desc = header[1:].strip()          # drop '>'
            cbc_id = desc.split()[0] if desc else ""
            meta = parse_description(desc)
            records.append(ConsensusRecord(
                cbc_id=cbc_id,
                seq="".join(seq_parts),
                description=desc,
                coverage=meta["coverage"],
                length=meta["length"],
            ))

    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            flush()
            header = line
            seq_parts = []
        elif line.strip():
            seq_parts.append(line.strip())
    flush()
    return records


def select_sequences(records: list[ConsensusRecord],
                     start: int | None = None, end: int | None = None,
                     config: ConsensusConfig | None = None) -> list[ConsensusRecord]:
    """Filter records by coverage and gap count (pure). (was: selectSeqs)

    Keeps a record if coverage > depth_min AND it has fewer than max_gaps '-'
    characters within the analysis window. Sequences are NOT trimmed -- they are
    kept full-length so variant positions stay genome-relative (required for
    region annotation). The window only bounds the GAP CHECK (and later, which
    positions get called), never renumbers coordinates.

    Args:
        start, end: optional 0-based analysis window [start, end). If both None,
            the whole sequence is the window (whole-reference analysis).
    """
    config = config or ConsensusConfig()
    kept: list[ConsensusRecord] = []
    for r in records:
        depth = r.coverage if r.coverage is not None else 0
        lo = 0 if start is None else start
        hi = len(r.seq) if end is None else end
        gaps_in_region = sum(1 for i in range(lo, min(hi, len(r.seq)))
                             if r.seq[i] == "-")
        if depth > config.depth_min and gaps_in_region < config.max_gaps_in_region:
            # keep full-length; no trimming
            kept.append(r)
    return kept


def _column_consensus(sequences: list[str]) -> str:
    """Per-column most-frequent character across equal-length sequences.

    Reproduces the maxChars computation: for each column, the character with the
    highest count (np.unique + argmax), matching the original's tie-breaking
    (argmax returns the first maximal index over the unique() ordering).
    """
    columns = np.transpose([list(s) for s in sequences])
    out = []
    for col in columns:
        chars, counts = np.unique(col, return_counts=True)
        out.append(chars[np.argmax(counts)])
    return "".join(out)


def genotype_summary(sequences: list[str], reference: str | None = None,
                     window: tuple[int, int] | None = None
                     ) -> tuple[list[str], str]:
    """Call per-sequence genotypes against a reference (PURE, the core).

    Sequences are full-length (untrimmed); variant positions are reported as
    1-based GENOME coordinates, so downstream region annotation can place each
    variant in its region(s).

    Args:
        sequences: equal-length aligned full-length sequences.
        reference: optional reference of the same length. If None, the per-column
            consensus is computed and used.
        window: optional (start, end) 0-based half-open range restricting WHICH
            columns are eligible to be called (to exclude ragged/low-coverage
            flanks). Positions are still numbered genome-relative regardless.
            If None, the whole sequence is eligible.

    Returns:
        (genotypes, reference_used). genotypes[i] is the "_"-joined token string
        for sequences[i]; tokens are "{genomePos}{base}" with genomePos 1-based.
    """
    if not sequences:
        return [], reference or ""

    columns = np.transpose([list(s) for s in sequences])
    n_cols = len(columns)

    lo = 0 if window is None else max(0, window[0])
    hi = n_cols if window is None else min(n_cols, window[1])

    # variant sites: columns (within the eligible window) where >1 char appears
    variant_sites = [i for i in range(lo, hi)
                     if len(np.unique(columns[i])) > 1]

    ref = reference if reference is not None else _column_consensus(sequences)

    genotypes = []
    for seq in sequences:
        # position i (0-based column) -> 1-based genome coordinate (i+1)
        tokens = [f"{i + 1}{seq[i]}" for i in variant_sites if seq[i] != ref[i]]
        genotypes.append("_".join(tokens))
    return genotypes, ref


def run(fasta: str, start: int | None = None, end: int | None = None,
        reference: str | None = None,
        config: ConsensusConfig | None = None,
        out_prefix: str | None = None) -> dict:
    """Full consensus stage: read, filter, genotype (genome-relative), write.

    Args:
        fasta: path to *_allConsensus.fasta.
        start, end: OPTIONAL analysis window (1-based inclusive, as a user would
            specify). Restricts which positions are called (to exclude ragged
            flanks) but does NOT trim sequences or renumber coordinates -- variant
            positions are always 1-based genome coordinates. If both None, the
            whole reference is analyzed (so non-coding variants surface).
        reference: optional reference; None computes the consensus.
        config: ConsensusConfig; defaults to ConsensusConfig().
        out_prefix: base path for outputs.

    Returns dict with 'reference', 'records', and written paths.
    """
    config = config or ConsensusConfig()
    records = parse_consensus_fasta(fasta)

    # Convert optional 1-based inclusive window to 0-based half-open for internals.
    win0 = None
    if start is not None and end is not None:
        win0 = (start - 1, end)

    slice_start = None if win0 is None else win0[0]
    slice_end = None if win0 is None else win0[1]
    kept = select_sequences(records, slice_start, slice_end, config)

    genotypes, ref_used = genotype_summary(
        [r.seq for r in kept], reference, window=win0)

    rows = [{
        ConsensusColumns.CBC_ID: r.cbc_id,
        ConsensusColumns.GENOTYPE: g,
        ConsensusColumns.SEQUENCE: r.seq,
        ConsensusColumns.DESCRIPTION: r.description,
    } for r, g in zip(kept, genotypes)]

    written = {}
    if out_prefix is None:
        out_prefix = str(fasta).replace("_allConsensus.fasta", "")

    ref_path = Path(f"{out_prefix}_consensus_reference.txt")
    ref_path.write_text(ref_used)
    written["reference"] = ref_path

    import csv as _csv
    csv_path = Path(f"{out_prefix}_filtConsensus.csv")
    with open(csv_path, "w", newline="") as fh:
        writer = _csv.DictWriter(fh, fieldnames=ConsensusColumns.ORDER)
        writer.writeheader()
        writer.writerows(rows)
    written["csv"] = csv_path

    return {"reference": ref_used, "records": rows, "written": written}
