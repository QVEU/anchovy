"""
config.py -- every tunable parameter for the pipeline, in one place.

WHY THIS FILE EXISTS
--------------------
In the original scripts the same kind of number was hardcoded in many spots:
the distance cutoff (42) and UMI offsets (38/10) lived inside anchovy.py, the
minimum reads per cell (5) inside CBCtoFasta.py, the depth minimum (10) inside
ConsensusTool.py, and so on. To run on a different 10X chemistry, a noisier
Nanopore dataset, or a different ORF, a researcher had to grep through source
in three languages.

Here every one of those becomes a named field with (a) a default that exactly
reproduces the current behavior and (b) a docstring explaining what it does and
when you'd change it. Nothing about behavior changes by introducing this file --
the defaults ARE the old literals. What changes is that the values are now
discoverable, documented, and overridable from the CLI or the workflow config.

DESIGN NOTE
-----------
We use frozen dataclasses. "Frozen" makes a config object immutable once built,
so a parameter can't be accidentally mutated halfway through a run -- a config
should be decided once, up front, then only read. Grouping by pipeline stage
(Extract / Fasta / Consensus) keeps each command's knobs together and mirrors
the module layout.

STRUCTURAL vs TUNABLE
---------------------
Structural constants that define the data format (column names, the
'coverage:'/'length:' description keys) live in schema.py, not here. This file
holds only values a user might legitimately want to change per run. The 10X
signature sits on the line between the two; it lives here because different
chemistries genuinely use different signatures, so it's a per-run choice.
"""

from __future__ import annotations

from dataclasses import dataclass, field


# The default 10X signature (v2/v3 3' chemistry). The N-run is the CBC+UMI region.
# Different chemistries use different signatures, so this is a per-run parameter,
# not a fixed structural constant.
DEFAULT_TENX_SIGNATURE = (
    "CTACACGACGCTCTTCCGATCT"
    "NNNNNNNNNNNNNNNNNNNNNNNNNN"
    "TTTCTTATAT"
)


@dataclass(frozen=True)
class ExtractConfig:
    """Parameters for `anchovy extract` (SAM -> per-read CBC/UMI table).

    Migrated from anchovy.py, where these were bare literals inside the
    functions noted in each comment.
    """

    # The 10X barcode signature searched for in each read. (was: tenXsignature)
    signature: str = DEFAULT_TENX_SIGNATURE

    # Number of worker processes for the multiprocessing pools.
    # (was: nthreads=16 default on poolBlocks / cellIDPool)
    nthreads: int = 16

    # How far upstream of the mapped offset to search for the signature, in nt.
    # (was: the literal 200 in poolBlocks' slice c[10][max(0,(c[14]-200)):c[14]])
    upstream_window: int = 200

    # Keep only candidate reads whose best signature match is closer than this
    # Levenshtein distance before attempting barcode assignment.
    # (was: pdSam[pdSam.minD<42] in cellIDPool)
    min_distance_cutoff: int = 42

    # Offsets used to slice the UMI out of a matched signature block:
    # matchseq[umi_start : len(matchseq) - umi_end_trim].
    # (was: matchseq[38:(len(matchseq)-10)] in cellMatch)
    umi_start_offset: int = 38
    umi_end_trim: int = 10

    # Maximum number of errors tolerated INSIDE the 16 nt barcode before a read
    # is dropped rather than assigned to a cell.
    #
    # None (the default) keeps every read, which is the original behavior: the
    # barcode search returns the NEAREST whitelist entry with no floor, so a
    # read whose barcode region is noise is still assigned to some cell, and
    # goes on to become per-cell reads, a consensus, and a genotype. On a real
    # run 43.7% of reads had no exact whitelist barcode, and all of them were
    # assigned anyway.
    #
    # Counted in ERRORS, not in raw Levenshtein distance, because the distance
    # a perfect match scores is the width of the UMI: every block is padded
    # with one N per UMI base, and an N never equals a real base. That is 10
    # for v2 and 12 for v3, so a raw threshold would silently mean something
    # different per chemistry. 0 admits only exact barcodes, 1 allows a single
    # substitution (the usual 10X correction), and so on.
    max_barcode_errors: int | None = None

    # Minimum read length filter, applied as length > this value. The original
    # used the query/signature length itself (minL = quL) as the threshold.
    # None means "use len(signature)", preserving the original behavior exactly;
    # set an explicit int to override.
    min_read_length: int | None = None

    # Reads held in memory at once. The stage reads the SAM, runs both passes
    # and builds its output a chunk at a time, so this -- not the size of the
    # run -- is one of the two terms in its footprint:
    #
    #     GB = W + nthreads x ( W + 4.6 KB x chunk_size )
    #
    # where W is the whitelist's lookup tables, paid once per worker. See
    # extract.projected_memory_gb, which is what both the stage and the
    # workflow's SLURM reservation compute this from.
    #
    # CHUNK SIZE IS THE SMALLER KNOB ON A FULL WHITELIST. At v3's 6.8 million
    # barcodes W is 2.75 GB against the chunk term's 0.46 GB, so nthreads is
    # what to lower on a tight node; chunking only dominates on a small
    # run-specific list. There is little to gain from raising this either way,
    # because the pool is already saturated well below it.
    chunk_size: int = 100_000

    def effective_min_read_length(self) -> int:
        """Resolve the read-length threshold, defaulting to signature length.

        This reproduces the original `minL = quL` behavior while still letting a
        user pin an explicit value. Kept as a method so the fallback logic lives
        with the config rather than being re-derived at each call site.
        """
        if self.min_read_length is None:
            return len(self.signature)
        return self.min_read_length


@dataclass(frozen=True)
class FastaConfig:
    """Parameters for `anchovy fasta` (per-read table -> per-cell FASTAs).

    Migrated from CBCtoFasta.py.
    """

    # Only emit a FASTA for barcodes supported by at least this many reads.
    # (was: the literal 5 in `if len(anchout.CBC[anchout.CBC==j])>=5`)
    min_reads_per_cbc: int = 5


@dataclass(frozen=True)
class ConsensusConfig:
    """Parameters for `anchovy consensus` (filter + genotype summary).

    Migrated from ConsensusTool.py. Note: start/end (the ORF/region bounds) are
    NOT here -- they're per-invocation positional arguments, not settings with a
    sensible default, so they stay as CLI arguments.
    """

    # Minimum per-sequence coverage/depth to keep a consensus, applied as
    # depth > this value. (was: depthMin=10 passed to selectSeqs in main)
    depth_min: int = 10

    # Minimum sequence length to keep. (was: lengthMin=1)
    length_min: int = 1

    # Maximum number of gap ('-') characters tolerated within the target region
    # before a sequence is filtered out. (was: the literal 3 in the selectSeqs
    # comprehension: sum(... == "-" ...) < 3)
    max_gaps_in_region: int = 3

    # --- breadth and depth-where-called (the unconflated filters) ----------- #
    # WHY THESE EXIST. `depth_min` is applied to the `coverage:` field that
    # sam2consensus writes into each consensus header, and that number is
    #
    #     (sum of depth over positions WITH READS) / (FULL reference length)
    #
    # -- the numerator skips uncovered positions, the denominator does not. So
    # it is not a depth at all but a PRODUCT of two different things: how deep
    # the reads were, and how much of the genome they spanned. A cell with 40x
    # over 40% of the genome scores 16 and is dropped; a cell with 21x across
    # the whole genome scores 21 and is kept -- despite having worse calls at
    # every position it calls.
    #
    # That is not a hypothetical. On the EV-A71 passage run, of 382 cells cut
    # by depth_min: 20, a full 122 had BOTH >20x where they called a base AND
    # >=50% of the genome covered. They were discarded for an averaging
    # artifact, not for quality -- a 27% loss of the usable population.
    #
    # Nor does depth_min protect anything the pipeline needs protecting from:
    # cons_min_depth already gap-fills any position under its threshold, so
    # every CALLED base is backed by real reads, and genotype_summary runs with
    # skip_gaps=True in whole-reference mode, so uncovered positions cannot
    # manufacture a variant in either direction.
    #
    # Splitting the metric lets each half be set for what it actually protects:
    #
    #   min_depth_called -- call quality WHERE a base was called.
    #   min_breadth      -- how much of the genome the cell saw. This is the
    #                       half that matters for the genotype and epistatic
    #                       NETWORKS: a cell covering 40% emits a short token
    #                       list that is indistinguishable from a fully covered
    #                       cell which happens to be clean in the missing
    #                       regions, so admitting cells of wildly different
    #                       breadth makes genotype identity partly an artifact
    #                       of what got sequenced. Keep this high when the
    #                       networks are the output you care about; lower it
    #                       when you are hunting variants at particular sites
    #                       and a deep partial cell is still evidence.
    #
    # Both default to None (off), so every existing config behaves exactly as
    # it did. Turning them on usually means setting `depth_min: 0`, since
    # depth_min is roughly their product and would otherwise re-impose the
    # conflated cutoff underneath the two clean ones.
    min_breadth: float | None = None
    min_depth_called: float | None = None


@dataclass(frozen=True)
class AnnotationConfig:
    """Parameters for the annotation stage (anchovy.annotate).

    Previously described an R stage and two settings that controlled nothing:
    `plot_haplotypes` (the port dropped plotting) and `build_network` (which
    duplicated annotate.run's own `network` argument). Both removed.
    """

    # Cells with this many called mutations or more are dropped as likely
    # artifacts. (was: mutantTable[BCMutCount < 200] in Consensus_Annotation.R,
    # and note the bound is exclusive -- a cell AT this count is dropped.)
    max_mutations_per_cell: int = 200


@dataclass(frozen=True)
class PipelineConfig:
    """Top-level config aggregating every stage.

    A single object to pass around. Each stage config is built with its own
    defaults unless overridden.
    """

    extract: ExtractConfig = field(default_factory=ExtractConfig)
    fasta: FastaConfig = field(default_factory=FastaConfig)
    consensus: ConsensusConfig = field(default_factory=ConsensusConfig)
    annotation: AnnotationConfig = field(default_factory=AnnotationConfig)
