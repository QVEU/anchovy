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


def select_sequences(records: list[ConsensusRecord], start: int, end: int,
                     config: ConsensusConfig | None = None) -> list[ConsensusRecord]:
    """Filter and trim records to the target region (pure). (was: selectSeqs)

    Keeps a record if coverage > depth_min AND it has fewer than max_gaps '-'
    characters within [start, end); surviving records are trimmed to [start:end].
    Reproduces selectSeqs, with thresholds from config instead of literals.
    """
    config = config or ConsensusConfig()
    kept: list[ConsensusRecord] = []
    for r in records:
        depth = r.coverage if r.coverage is not None else 0
        gaps_in_region = sum(1 for i in range(start, end) if r.seq[i] == "-")
        if depth > config.depth_min and gaps_in_region < config.max_gaps_in_region:
            trimmed = ConsensusRecord(
                cbc_id=r.cbc_id,
                seq=r.seq[start:end],
                description=r.description,
                coverage=r.coverage,
                length=r.length,
            )
            kept.append(trimmed)
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


def genotype_summary(sequences: list[str], reference: str | None = None
                     ) -> tuple[list[str], str]:
    """Call per-sequence genotypes against a reference (PURE, the core).

    This is the unified replacement for the original's if/else. When `reference`
    is None, the reference is the computed per-column consensus (the original's
    working if-branch); otherwise the supplied reference is used (the new,
    previously-broken else-branch, now functional).

    Args:
        sequences: equal-length aligned sequences (already trimmed to the region).
        reference: optional reference sequence of the same length. If None, the
            per-column consensus of `sequences` is computed and used.

    Returns:
        (genotypes, reference_used) where genotypes[i] is the "_"-joined mutation
        token string for sequences[i], and reference_used is the reference the
        genotypes were called against.

    Variant sites are columns where more than one character appears across
    `sequences`; a token is emitted for a sequence at a variant site only where
    it differs from the reference. Positions are 1-based within the region.
    """
    if not sequences:
        return [], reference or ""

    columns = np.transpose([list(s) for s in sequences])
    variant_sites = [i for i in range(len(columns))
                     if len(np.unique(columns[i])) > 1]

    ref = reference if reference is not None else _column_consensus(sequences)

    genotypes = []
    for seq in sequences:
        tokens = [f"{i + 1}{seq[i]}" for i in variant_sites if seq[i] != ref[i]]
        genotypes.append("_".join(tokens))
    return genotypes, ref


def run(fasta: str, start: int, end: int, reference: str | None = None,
        config: ConsensusConfig | None = None,
        out_prefix: str | None = None) -> dict:
    """Full consensus stage: read, filter/trim, genotype, and write outputs.

    Args:
        fasta: path to *_allConsensus.fasta.
        start, end: region of interest (nt).
        reference: optional reference; None computes the consensus (original path).
        config: ConsensusConfig; defaults to ConsensusConfig().
        out_prefix: base path for outputs; defaults to the fasta path with
            '_allConsensus.fasta' stripped (matching the original's naming).

    Returns:
        dict with 'reference', 'records' (list of dicts with CBC_ID/genotype/
        sequence/description), and the paths written.
    """
    config = config or ConsensusConfig()
    records = parse_consensus_fasta(fasta)
    kept = select_sequences(records, start, end, config)

    genotypes, ref_used = genotype_summary([r.seq for r in kept], reference)

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
