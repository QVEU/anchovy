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
  - variant sites = columns where any sequence differs from the reference. The
    original asked whether the column varied across the cells, which is the SAME
    SET when the reference is the computed consensus (the consensus is always one
    of the characters present), but differs for a SUPPLIED reference -- see
    genotype_summary().
  - per-column consensus = most frequent character (argmax of unique counts)
  - each sequence's genotype = "_".join((i+1)+base) over variant sites where the
    sequence differs from the per-column consensus, 1-based within the trimmed region
Filtering reproduces selectSeqs: keep depth > depth_min AND (< max_gaps '-' in the
trimmed region), then trim to [start:end].

3. WHOLE-REFERENCE MODE (run(trim=False)) -- added for region-aware annotation.
   The legacy path trims to the ORF and renumbers positions from 1 within it,
   which throws away the genome coordinates the GFF3-driven annotate stage needs.
   With trim=False the full reference length is kept, so genotype token positions
   ARE genome coordinates, and start/end degrade from a cut to an optional
   analysis window. Legacy remains the default, so every frozen golden still
   passes unchanged.
"""

from __future__ import annotations

import warnings
from collections import Counter
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


def select_sequences(records: list[ConsensusRecord], start: int | None = None,
                     end: int | None = None,
                     config: ConsensusConfig | None = None,
                     trim: bool = True,
                     stats: dict | None = None) -> list[ConsensusRecord]:
    """Filter (and in legacy mode trim) records to the target region (pure).

    (was: selectSeqs)

    Two modes, selected by `trim`:

    trim=True (DEFAULT, the legacy path)
        Reproduces selectSeqs exactly: keep a record if coverage > depth_min AND
        it has fewer than max_gaps '-' characters within [start, end), then cut
        the sequence down to [start:end]. Downstream positions are therefore
        numbered from 1 WITHIN the trimmed region.

    trim=False (whole-reference mode)
        The sequence is kept at full reference length, so downstream positions
        are genome coordinates that never renumber. `start`/`end` become an
        optional ANALYSIS WINDOW rather than a cut.

    THE GAP FILTER IS WINDOW-CONDITIONAL, and that is load-bearing. sam2consensus
    emits a full-reference-length consensus and gap-fills every position that had
    no read coverage, so a whole-genome gap count is dominated by the ragged
    uncovered flanks. Applying max_gaps to that count filters out every cell. The
    filter is therefore applied only when a window is actually given -- which is
    the point of having a window: it names the covered core to judge cells on.

    BREADTH AND DEPTH-WHERE-CALLED (config.min_breadth / min_depth_called)
    ---------------------------------------------------------------------
    Optional, off unless set, and applied in addition to depth_min. They exist
    because depth_min tests a number that conflates depth with breadth -- see
    the long note on ConsensusConfig for why that silently discards deep cells
    that happen to span less of the genome.

    Both are derived from what is already on hand, so nothing needs recomputing
    upstream. With `ncols = len(r.seq)` (the full alignment length, which is
    what sam2consensus divided by):

        n_called     = positions holding a real base, i.e. not '-'
        breadth      = n_called / ncols       (over the WINDOW when given)
        depth_called = coverage * ncols / n_called

    The second identity holds because `coverage * ncols` recovers sam2consensus'
    `sumcov` -- the summed depth it accumulated. NOTE IT IS A SLIGHT
    OVERESTIMATE: sumcov accumulates at every position with AT LEAST ONE read,
    while n_called counts only positions that reached cons_min_depth, so depth
    from the 1..cons_min_depth-1 positions lands in the numerator while those
    positions are excluded from the denominator. The bias is upward and small
    (those positions are shallow by construction), and it is stated here rather
    than papered over because the alternative -- the true covered-position
    count -- does not survive into the FASTA and cannot be recovered from it.

    Breadth is restricted to [start, end) when a window is given, matching the
    gap filter above -- see the comment at the computation. depth_called never
    is: `coverage` is a whole-reference number and the per-position depths that
    would let it be windowed are not recoverable from the FASTA, so windowing
    only its denominator would compare incomparable scales.
    """
    config = config or ConsensusConfig()
    has_window = start is not None and end is not None
    check_shape = (config.min_breadth is not None
                   or config.min_depth_called is not None)
    n_low_depth = n_narrow = n_shallow = n_gappy = 0
    kept: list[ConsensusRecord] = []
    for r in records:
        depth = r.coverage if r.coverage is not None else 0
        if depth <= config.depth_min:
            n_low_depth += 1
            continue
        if check_shape:
            ncols = len(r.seq)
            # BREADTH IS WINDOW-CONDITIONAL, for the same reason the gap filter
            # below is: a whole-genome count is dominated by the ragged
            # uncovered flanks that every amplicon run has, so "fully covered"
            # over the full reference is a bar essentially no cell clears.
            # Given a window it means fully covered ACROSS THE CORE, which is
            # the question worth asking and the point of naming a window.
            lo, hi = (start, end) if has_window else (0, ncols)
            span = hi - lo
            n_in_window = sum(1 for i in range(lo, hi) if r.seq[i] != "-")
            breadth = (n_in_window / span) if span else 0.0
            # DEPTH-CALLED IS NOT, and cannot be: `coverage` is sumcov over the
            # FULL length, and the per-position depths needed to restrict it to
            # a window do not survive into the FASTA. So it stays a
            # whole-sequence number even when breadth is windowed.
            # coverage == sumcov / ncols, so sumcov == depth * ncols.
            n_called = sum(1 for ch in r.seq if ch != "-")
            depth_called = (depth * ncols / n_called) if n_called else 0.0
            if config.min_breadth is not None and breadth < config.min_breadth:
                n_narrow += 1
                continue
            if (config.min_depth_called is not None
                    and depth_called < config.min_depth_called):
                n_shallow += 1
                continue
        if has_window:
            gaps_in_region = sum(1 for i in range(start, end) if r.seq[i] == "-")
            if gaps_in_region >= config.max_gaps_in_region:
                n_gappy += 1
                continue
        kept.append(ConsensusRecord(
            cbc_id=r.cbc_id,
            seq=r.seq[start:end] if trim else r.seq,
            description=r.description,
            coverage=r.coverage,
            length=r.length,
        ))
    if stats is not None:
        stats["dropped_low_depth"] = n_low_depth
        stats["dropped_narrow"] = n_narrow
        stats["dropped_shallow_called"] = n_shallow
        stats["dropped_gappy"] = n_gappy
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


# A variant call requires a definite base on both sides. Everything else that
# sam2consensus can emit at a position -- an IUPAC ambiguity code (R, Y, S...),
# its lowercase form (a gap or N was among the observed bases), or N itself --
# records UNCERTAINTY, not a difference.
CALLED_BASES = frozenset("ACGT")


def genotype_summary(sequences: list[str], reference: str | None = None,
                     window: tuple[int, int] | None = None,
                     skip_gaps: bool = False,
                     skip_ambiguous: bool = False,
                     stats: dict | None = None) -> tuple[list[str], str]:
    """Call per-sequence genotypes against a reference (PURE, the core).

    This is the unified replacement for the original's if/else. When `reference`
    is None, the reference is the computed per-column consensus (the original's
    working if-branch); otherwise the supplied reference is used (the new,
    previously-broken else-branch, now functional).

    Args:
        sequences: equal-length aligned sequences. In the legacy path these are
            already trimmed to the region; in whole-reference mode they are full
            reference length, so token positions ARE genome coordinates.
        reference: optional reference sequence of the same length. If None, the
            per-column consensus of `sequences` is computed and used.
        window: optional (lo, hi) 0-based half-open ANALYSIS WINDOW. Only columns
            inside it may become variant sites. It restricts WHICH columns are
            called, never how they are NUMBERED -- positions stay 1-based over
            the sequences as given, so a window never renumbers anything.
        skip_gaps: when True, never emit a token at a column where either the
            sequence or the reference holds a '-'. See the note below.
        skip_ambiguous: when True, never emit a token where either side is not
            a definite A/C/G/T. See the second note below.
        stats: optional dict, populated with counts of what was skipped, so a
            caller can report it rather than dropping data in silence.

    Returns:
        (genotypes, reference_used) where genotypes[i] is the "_"-joined mutation
        token string for sequences[i], and reference_used is the reference the
        genotypes were called against.

    A variant site is any column where at least one sequence differs from the
    reference, so a mutation FIXED across every cell is still reported. Positions
    are 1-based over the input columns.

    WHY skip_ambiguous EXISTS
    -------------------------
    sam2consensus does not only emit A/C/G/T. Where a position's reads disagree
    it emits an IUPAC ambiguity code -- R for A-or-G, Y for C-or-T, and so on --
    and a lowercase form of it when a gap or an N was among the observed bases.
    Those record that the cell's reads DISAGREED, not that the cell carries a
    mutation.

    Treated as variants they are actively misleading, because the genotype
    machinery works on token identity: two cells both reading R at position 3185
    are grouped as sharing a mutation, when what they share is uncertainty. On a
    real run that mattered -- 6 of 11 distinct tokens were ambiguity codes, so
    over half the genotype network was built from positions nobody had called.

    They are also self-inflicted at low coverage: at 3-5x a single discordant
    read can stop any base reaching the consensus threshold. More reads, not a
    different threshold, is what resolves them -- raising the threshold emits
    MORE ambiguity, since the algorithm accumulates bases until their combined
    coverage reaches it.

    Dropping them is the conservative reading and the default in
    whole-reference mode. It does discard any genuine within-cell mixed
    population, which for a viral quasispecies could be real signal, so the
    count is reported rather than silently absorbed, and it can be turned off.

    WHY skip_gaps EXISTS (whole-reference mode only, hence default False)
    --------------------------------------------------------------------
    sam2consensus gap-fills positions with no read coverage. Once the reference
    is kept whole, those '-' columns vary across cells purely because cells were
    sequenced over different spans. Calling "201-" (cell has a gap) or "201A"
    against a '-' reference (reference has the gap) reports a coverage
    difference as a substitution, which it is not. A substitution call requires
    a real base on both sides, so both directions are skipped. The legacy
    trimmed path keeps its original behavior untouched.
    """
    if not sequences:
        return [], reference or ""

    ncols = len(sequences[0])
    lo, hi = (0, ncols) if window is None else (max(0, window[0]), min(ncols, window[1]))

    if reference is not None and len(reference) != ncols:
        raise ValueError(
            f"reference length {len(reference)} does not match the alignment "
            f"length {ncols}. Every position is compared by index, so a "
            f"mismatched reference would shift or truncate the calls. In "
            f"whole-reference mode pass the full genome; in the default trimmed "
            f"mode pass a reference already cut to the same region.")

    ref = reference if reference is not None else _column_consensus(sequences)

    # A site is worth calling wherever ANY sequence differs from the reference.
    #
    # This used to ask a different question -- whether the column varied ACROSS
    # THE CELLS -- which quietly discarded FIXED DIFFERENCES: a position where
    # every cell agrees with the others but none matches the supplied reference
    # produced no call at all, so the whole population read as wild-type there.
    # On a passaged or lab-adapted stock that can be most of the real variants.
    #
    # For the computed-reference path the two questions are provably the same,
    # because the column consensus is always one of the characters present, so a
    # column varies if and only if some sequence differs from it. Verified over
    # 20000 randomized alignments (including gap-bearing alphabets): zero
    # disagreements. That is why widening this cannot move any reference=None
    # output, and why the frozen goldens still pass unchanged.
    variant_sites = [i for i in range(lo, hi)
                     if any(seq[i] != ref[i] for seq in sequences)]

    n_gap = n_ambiguous = 0
    genotypes = []
    for seq in sequences:
        tokens = []
        for i in variant_sites:
            if seq[i] == ref[i]:
                continue
            if skip_gaps and (seq[i] == "-" or ref[i] == "-"):
                n_gap += 1
                continue
            if skip_ambiguous and not (seq[i] in CALLED_BASES
                                       and ref[i] in CALLED_BASES):
                n_ambiguous += 1
                continue
            tokens.append(f"{i + 1}{seq[i]}")
        genotypes.append("_".join(tokens))

    if stats is not None:
        stats["gap_calls_skipped"] = n_gap
        stats["ambiguous_calls_skipped"] = n_ambiguous

    return genotypes, ref


def run(fasta: str, start: int | None = None, end: int | None = None,
        reference: str | None = None,
        config: ConsensusConfig | None = None,
        out_prefix: str | None = None,
        trim: bool = True,
        keep_ambiguous: bool = False) -> dict:
    """Full consensus stage: read, filter (trim), genotype, and write outputs.

    Args:
        fasta: path to *_allConsensus.fasta.
        start, end: region of interest (nt, 0-based half-open). Their meaning
            depends on `trim` -- see below.
        reference: optional reference; None computes the consensus (original path).
        config: ConsensusConfig; defaults to ConsensusConfig().
        out_prefix: base path for outputs; defaults to the fasta path with
            '_allConsensus.fasta' stripped (matching the original's naming).
        keep_ambiguous: whole-reference mode only. By default a position is
            called only where both the cell and the reference hold a definite
            A/C/G/T, so IUPAC ambiguity codes do not become mutations. Set True
            to keep them -- relevant if you are after genuine within-cell mixed
            populations rather than clean per-cell genotypes. How many calls
            this dropped is reported in the return value and printed.
        trim: True (default) keeps the legacy behavior exactly -- cut every
            sequence to [start:end] and number genotype tokens from 1 within that
            region. False selects WHOLE-REFERENCE MODE: sequences stay full
            length, token positions are genome coordinates that never renumber,
            and start/end (if given) act only as an analysis window that excludes
            ragged flanks from variant calling.

    Returns:
        dict with 'reference', 'records' (list of dicts with CBC_ID/genotype/
        sequence/description), and the paths written.

    WHOLE-REFERENCE MODE is what makes region-aware annotation possible: the
    annotate stage needs genome coordinates to look mutations up in a GFF3, and
    the legacy trimmed path renumbers positions relative to the ORF, destroying
    exactly that information.
    """
    config = config or ConsensusConfig()
    records = parse_consensus_fasta(fasta)
    stats: dict = {}
    kept = select_sequences(records, start, end, config, trim=trim, stats=stats)
    stats["input_records"] = len(records)

    # EQUAL LENGTH IS ASSUMED EVERYWHERE BELOW AND GUARANTEED BY NOTHING.
    # sam2consensus appends insertion columns to the consensus it emits, so a
    # cell carrying a called insertion comes out LONGER than the reference.
    # Genotypes are called by index against ref[i], so every position after
    # that insertion is shifted, and the shift reads as a run of substitutions
    # covering the whole rest of the genome.
    #
    # Nothing caught it. The only length check compares the reference against
    # sequences[0], so whether a run raised ValueError or silently mis-called
    # depended on the length of whichever cell happened to come first -- and
    # annotate's max_mutations_per_cell cap then discarded the shifted cells as
    # "hypermutated", which is precisely what a frame shift looks like. The
    # failure was therefore invisible from the outputs: the cells did not
    # appear, and nothing said why.
    #
    # Insertion columns cannot be identified from the FASTA alone (sam2consensus
    # does not mark them), so a shifted cell cannot be repaired here -- only
    # recognised. Coordinates that cannot be trusted must not produce
    # genotypes, so those cells are dropped, counted, and reported.
    #
    # When every length already agrees this is a no-op, which is why the
    # goldens are untouched.
    if kept:
        expected = (len(reference) if reference is not None
                    else Counter(len(r.seq) for r in kept).most_common(1)[0][0])
        conformant = [r for r in kept if len(r.seq) == expected]
        stats["expected_length"] = expected
        stats["dropped_length_mismatch"] = len(kept) - len(conformant)
        kept = conformant
    stats["kept"] = len(kept)
    if trim:
        genotypes, ref_used = genotype_summary([r.seq for r in kept], reference)
    else:
        window = (start, end) if start is not None and end is not None else None
        genotypes, ref_used = genotype_summary(
            [r.seq for r in kept], reference, window=window, skip_gaps=True,
            skip_ambiguous=not keep_ambiguous, stats=stats)

    # Loud, because the old behaviour here was to lose these cells without a
    # word -- as a ValueError with an unrelated message, or as silent
    # hypermutants two stages downstream.
    if stats.get("dropped_length_mismatch"):
        print(f"Dropped {stats['dropped_length_mismatch']} cell(s) whose consensus "
              f"is not {stats['expected_length']} nt. sam2consensus emits called "
              f"insertions as extra columns, which shift every position after "
              f"them out of genome coordinates, so these cells cannot be "
              f"genotyped by position. Raise cons_min_depth if the insertions "
              f"are low-depth artifacts.")

    # Same reasoning as the ambiguity report below: a filter that silently
    # halves the cell count is the kind of thing that gets mistaken for a
    # biological result, so each cutoff reports its own toll separately rather
    # than leaving one aggregate "kept N" to be reverse-engineered.
    if stats.get("dropped_narrow") or stats.get("dropped_shallow_called"):
        print(f"Filtered {stats['input_records']} consensus sequences to "
              f"{stats['kept']}: {stats['dropped_low_depth']} below depth_min, "
              f"{stats['dropped_narrow']} below min_breadth, "
              f"{stats['dropped_shallow_called']} below min_depth_called.")

    # Say what was dropped. Filtering ambiguity codes is the right default, but
    # a variant quietly missing from a table is the kind of thing that costs
    # someone an afternoon, so the number is never left implicit.
    if stats.get("ambiguous_calls_skipped"):
        print(f"Skipped {stats['ambiguous_calls_skipped']} call(s) at positions "
              "with an ambiguous base (IUPAC code, N, or lowercase). These are "
              "positions where the reads disagreed, not called mutations; pass "
              "keep_ambiguous to retain them.")

    # A COMPUTED reference needs at least two cells to mean anything. With one,
    # the per-column consensus IS that cell, so no position can differ from it
    # and the genotype is empty by construction -- not because the cell matches
    # the virus, but because there was nothing to compare it against. Say so,
    # because an empty genotype column otherwise reads as a biological result.
    if reference is None and len(kept) < 2:
        warnings.warn(
            f"only {len(kept)} cell(s) passed filtering, and no reference was "
            f"supplied.\n"
            f"  The reference is computed as the consensus ACROSS cells, so with "
            f"fewer than two\n"
            f"  there is nothing to compare against and every genotype comes out "
            f"empty.\n"
            f"  Either lower depth_min to keep more cells, or pass an explicit "
            f"--reference\n"
            f"  (the genome FASTA) so each cell is called against that instead.",
            stacklevel=2)

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

    return {"reference": ref_used, "records": rows, "written": written,
            "stats": stats}
