"""
schema.py -- the single source of truth for data that flows between stages.

WHY THIS FILE EXISTS
--------------------
In the original pipeline, each stage knew the internal layout of the others:
anchovy.py wrote a CSV, CBCtoFasta.py read it by positional index (i[9], i[10]),
ConsensusTool.py packed metadata into FASTA description strings, and the R script
parsed those strings back apart by whitespace offset. Every one of those is an
implicit, unenforced contract. Change one column's order and everything downstream
silently corrupts.

This module makes the contract explicit and central. Every reader and writer in
the package imports these names. To add or rename a column, you edit it HERE, once,
and every stage stays consistent. The positional-index problem doesn't get patched --
it stops being possible, because nothing addresses data by position anymore.
"""

from __future__ import annotations


class AnchovyColumns:
    """Columns in the anchovy extract output (SAM -> per-read CBC/UMI table).

    Produced by extract.py, consumed by fasta.py. These names replace the old
    positional itertuples() indices that coupled the two scripts.
    """
    CBC = "CBC"                 # assigned cell barcode
    MIN_DISTANCE = "minD"       # Levenshtein distance of the best barcode match
    BARCODE_ID = "BC_ID"        # index of the matched barcode in the whitelist
    READ_POSITION = "readPos"   # position of the query match within the read
    CLIP_LENGTH = "clipLength"  # clipped read length
    RECON_QUERY = "reconQuery"  # reconstructed barcode-augmented query block
    MATCH_SEQ = "matchseq"      # the matched query block
    READ = "read"               # read ID  (was positional index 9 in CBCtoFasta)
    MAPPED_SEQ = "mappedSeq"    # the mapped sequence (was positional index 10)
    UMI = "UMI"                 # extracted UMI

    # Canonical output order. Writing through this guarantees the column order
    # is stable, so even code that hasn't migrated to named access stays correct.
    ORDER = [CBC, MIN_DISTANCE, BARCODE_ID, READ_POSITION, CLIP_LENGTH,
             RECON_QUERY, MATCH_SEQ, READ, MAPPED_SEQ, UMI]


class SamColumns:
    """The 11 standard SAM fields captured by the manual text parse in io.py."""
    READ = "read"
    FLAG = "FLAG"
    TEMPLATE = "template"
    POS = "pos"
    MAPQ = "mapq"
    CIGAR = "cigar"
    RNEXT = "Rnext"
    PNEXT = "Pnext"
    TLEN = "Tlen"
    SEQ = "seq"
    QSCORE = "Qscore"

    ORDER = [READ, FLAG, TEMPLATE, POS, MAPQ, CIGAR,
             RNEXT, PNEXT, TLEN, SEQ, QSCORE]

    # Columns added downstream by the CIGAR parse and the matching stage.
    READ_LEN = "readLen"
    CLIP_READ_LEN = "clipReadLen"
    OFFSET = "offset"
    LENGTH = "length"


class ConsensusColumns:
    """Columns in the filtered-consensus table (consensus.py -> R annotation).

    The R script's CBC_ID/genotype/description contract lives here so both the
    Python writer and the R reader agree on it. `description` still carries the
    'coverage:NN length:MM' convention; see io.parse_description() for the parser
    that both sides should use instead of ad-hoc string splitting.
    """
    CBC_ID = "CBC_ID"           # cell barcode / cell identifier
    GENOTYPE = "genotype"       # mutation tokens joined by "_", e.g. "12A_45T"
    SEQUENCE = "sequence"       # the consensus nucleotide sequence
    DESCRIPTION = "description" # free-text metadata (coverage/length)

    ORDER = [CBC_ID, GENOTYPE, SEQUENCE, DESCRIPTION]


# The 10X signature and description-field conventions were previously hardcoded
# as string literals in multiple files. Centralize the structural constants here;
# the *tunable* ones (cutoffs, thresholds) live in config.py instead.
DESCRIPTION_COVERAGE_KEY = "coverage:"
DESCRIPTION_LENGTH_KEY = "length:"
