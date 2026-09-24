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

    # What read_sam actually retains. QSCORE is the per-base quality string --
    # the same length as SEQ, so together they are the bulk of the frame -- and
    # nothing in the pipeline reads it. ORDER still describes the SAM line; this
    # is the subset kept in memory. Both are needed: the reader slices a parsed
    # line by len(KEPT), and ORDER says what it is slicing.
    KEPT = [READ, FLAG, TEMPLATE, POS, MAPQ, CIGAR,
            RNEXT, PNEXT, TLEN, SEQ]

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


class CellAlleleColumns:
    """Columns in the per-cell pileup table (frequencies.py).

    One row per (cell, position, non-reference allele) with read support. This
    is a WITHIN-CELL frequency: reads carrying the allele over read depth at
    that position, in that cell. It answers whether a cell is mixed at a site.
    """
    CBC_ID = "CBC_ID"
    POSITION = "position"        # 1-based genome coordinate
    REF_BASE = "ref_base"
    ALLELE = "allele"
    READS = "reads"              # reads carrying ALLELE in this cell
    DEPTH = "depth"              # reads covering POSITION in this cell
    FREQ = "freq"                # READS / DEPTH

    ORDER = [CBC_ID, POSITION, REF_BASE, ALLELE, READS, DEPTH, FREQ]


class AlleleFrequencyColumns:
    """Columns in the population allele-frequency table (frequencies.py).

    One row per (position, non-reference allele) observed in any cell's
    consensus. THE DENOMINATORS ARE PER POSITION, which is the point: a cell
    that covers part of the CDS counts where it has data and is simply absent
    where it does not, so partial cells contribute without being discarded and
    without inflating any denominator they cannot speak to.

    Three views of the same allele, deliberately side by side:

      *_cells    one cell, one vote, over every mapped cell. The population
                 frequency -- what fraction of cells carry this mutation.
      *_subset   the same over the cells that passed filtering (the ones the
                 genotype network is built from), so the subset can be checked
                 for bias against the full population rather than assumed
                 representative.
      *_reads    reads carrying the allele over reads covering the position,
                 summed across cells. NOTE this weights each cell by its depth,
                 and within a cell reads are amplification copies of a few
                 templates -- so it is a sequencing-level summary, not an
                 independent-observation frequency. Use *_cells for population
                 claims.
    """
    POSITION = "position"        # 1-based genome coordinate
    REGION = "region"            # GFF region name(s) containing POSITION
    REF_BASE = "ref_base"
    ALLELE = "allele"

    CELLS_ALT = "cells_alt"
    CELLS_CALLED = "cells_called"
    FREQ_CELLS = "freq_cells"

    CELLS_ALT_SUBSET = "cells_alt_subset"
    CELLS_CALLED_SUBSET = "cells_called_subset"
    FREQ_CELLS_SUBSET = "freq_cells_subset"

    READS_ALT = "reads_alt"
    READS_DEPTH = "reads_depth"
    FREQ_READS = "freq_reads"

    ORDER = [POSITION, REGION, REF_BASE, ALLELE,
             CELLS_ALT, CELLS_CALLED, FREQ_CELLS,
             CELLS_ALT_SUBSET, CELLS_CALLED_SUBSET, FREQ_CELLS_SUBSET,
             READS_ALT, READS_DEPTH, FREQ_READS]


# The 10X signature and description-field conventions were previously hardcoded
# as string literals in multiple files. Centralize the structural constants here;
# the *tunable* ones (cutoffs, thresholds) live in config.py instead.
DESCRIPTION_COVERAGE_KEY = "coverage:"
DESCRIPTION_LENGTH_KEY = "length:"


# --------------------------------------------------------------------------- #
# 10X signature layout
# --------------------------------------------------------------------------- #
# The signature itself is a per-run CHOICE and lives in config.py. Its SHAPE is
# structural, and belongs here:
#
#     <---- 22 ----><------- barcode 16 -------><-- UMI n --><--- 10 --->
#     CTACACGACG...  NNNNNNNNNNNNNNNN            NNNNNNNNNN   TTTCTTATAT
#      constant 5'         barcode                  UMI       constant 3'
#
# barcodes.py encoded these as bare numbers (22, 48, 10, and 38 = 22 + 16), and
# they are load-bearing: the UMI width is derived as len(signature) - 48, which
# is what lets one signature string carry the whole chemistry. Swap the 26-N v2
# signature for the 28-N v3 one and the UMI follows automatically, because the
# prefix, barcode and suffix lengths are identical between the two.
#
# That only holds while a signature actually HAS this shape. Naming the parts
# lets validate_signature() check it instead of assuming it.
SIGNATURE_PREFIX_LEN = 22    # constant 5' handle
SIGNATURE_BARCODE_LEN = 16   # cell barcode; 16 in both 10X v2 and v3
SIGNATURE_SUFFIX_LEN = 10    # constant 3' handle

# Offset of the UMI within a matched block, i.e. past the prefix and barcode.
SIGNATURE_UMI_START = SIGNATURE_PREFIX_LEN + SIGNATURE_BARCODE_LEN          # 38
# Everything that is not UMI: len(signature) minus this is the UMI width.
SIGNATURE_NON_UMI_LEN = SIGNATURE_UMI_START + SIGNATURE_SUFFIX_LEN          # 48
