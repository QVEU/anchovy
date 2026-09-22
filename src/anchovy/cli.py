"""
cli.py -- command-line entry point for the anchovy pipeline.

Wired to pyproject.toml's [project.scripts] as `anchovy = "anchovy.cli:main"`,
so after `pip install -e .` the shell command `anchovy` dispatches here.

Each subcommand is a thin adapter: parse args -> build the relevant config from
any overrides -> call the stage's run() -> report what was produced. All the real
logic lives in extract.py / fasta.py / consensus.py; this file only translates
between the command line and those functions.

DESIGN
------
- Every tunable exposed as an optional flag defaults to the config's default, so
  running with no flags reproduces the original behavior exactly. Flags only
  override when explicitly given.
- Commands return an int exit code (0 = success) so the shell and any wrapping
  workflow (Snakemake, later) can detect failure.
"""

from __future__ import annotations

import argparse
import sys

from anchovy import __version__
from anchovy.config import ExtractConfig, FastaConfig, ConsensusConfig


# --------------------------------------------------------------------------- #
# extract
# --------------------------------------------------------------------------- #
def _cmd_extract(args: argparse.Namespace) -> int:
    from anchovy import extract

    # Build config from overrides; unspecified flags keep config defaults.
    defaults = ExtractConfig()
    config = ExtractConfig(
        signature=args.signature or defaults.signature,
        nthreads=args.threads if args.threads is not None else defaults.nthreads,
        min_distance_cutoff=(args.max_distance if args.max_distance is not None
                             else defaults.min_distance_cutoff),
        max_barcode_errors=args.max_barcode_errors,
    )

    df = extract.run(sam=args.sam, whitelist=args.whitelist,
                     signature=config.signature, config=config)

    from anchovy.io import write_anchovy_csv
    out = args.out or args.sam.replace(".sam", "_anchovy.csv").replace(".bam", "_anchovy.csv")
    write_anchovy_csv(df, out)
    print(f"Wrote {out} ({len(df)} reads assigned).")
    return 0


# --------------------------------------------------------------------------- #
# fasta
# --------------------------------------------------------------------------- #
def _cmd_fasta(args: argparse.Namespace) -> int:
    from anchovy import fasta

    defaults = FastaConfig()
    config = FastaConfig(
        min_reads_per_cbc=(args.min_reads if args.min_reads is not None
                           else defaults.min_reads_per_cbc),
    )

    written = fasta.run(csv=args.csv, out_dir=args.outdir, config=config)
    print(f"Wrote {len(written)} per-cell FASTA(s) to {args.outdir}.")
    return 0


# --------------------------------------------------------------------------- #
# consensus
# --------------------------------------------------------------------------- #
def _cmd_consensus(args: argparse.Namespace) -> int:
    from anchovy import consensus

    defaults = ConsensusConfig()
    config = ConsensusConfig(
        depth_min=args.depth_min if args.depth_min is not None else defaults.depth_min,
        max_gaps_in_region=(args.max_gaps if args.max_gaps is not None
                            else defaults.max_gaps_in_region),
        # Left as None when unset, which is what switches these filters off --
        # so no default is substituted here the way it is above.
        min_breadth=args.min_breadth,
        min_depth_called=args.min_depth_called,
    )

    reference = None
    if args.reference:
        # Via the shared reader, so a FASTA reference works. A bare read() would
        # splice the '>' header into the sequence and shift every coordinate.
        from anchovy.io import read_reference_sequence
        reference = read_reference_sequence(args.reference)

    trim = not args.whole_reference
    if trim and (args.start is None or args.end is None):
        print("error: start and end are required unless --whole-reference is given.",
              file=sys.stderr)
        return 2

    result = consensus.run(
        fasta=args.fasta, start=args.start, end=args.end,
        reference=reference, config=config, out_prefix=args.out_prefix,
        trim=trim, keep_ambiguous=args.keep_ambiguous,
    )
    written = result["written"]
    print(f"Kept {len(result['records'])} sequences.")
    print(f"Wrote {written['reference']}")
    print(f"Wrote {written['csv']}")
    return 0


def _cmd_annotate(args: argparse.Namespace) -> int:
    from anchovy import annotate
    from anchovy.config import AnnotationConfig

    defaults = AnnotationConfig()
    config = AnnotationConfig(
        max_mutations_per_cell=(args.max_mutations if args.max_mutations
                                is not None else defaults.max_mutations_per_cell),
    )
    result = annotate.run(
        filt_consensus_csv=args.csv,
        reference_file=args.reference,
        out_prefix=args.out_prefix,
        network=not args.no_network,
        gff=args.gff,
        self_edges=args.self_edges,
        config=config,
    )
    for label, path in result["written"].items():
        print(f"Wrote {path}")
    return 0


# --------------------------------------------------------------------------- #
# parser
# --------------------------------------------------------------------------- #
def _cmd_frequencies(args: argparse.Namespace) -> int:
    from anchovy import frequencies
    from anchovy.io import read_reference_sequence

    reference = read_reference_sequence(args.reference)
    result = frequencies.run(
        fasta=args.fasta, reference=reference, out_prefix=args.out_prefix,
        counts_dir=args.counts_dir, reference_name=args.reference_name,
        subset_csv=args.subset, gff=args.gff,
        min_alt_reads=args.min_alt_reads, min_alt_freq=args.min_alt_freq,
        min_cells=args.min_cells,
    )
    stats = result["stats"]
    print(f"{stats['cells']} cells, {stats['cells_in_subset']} in the subset, "
          f"{stats['cells_with_counts']} with pileup counts.")
    if args.counts_dir and stats["cells_with_counts"] == 0:
        print("warning: --counts-dir was given but no counts file matched any "
              "cell. Check --reference-name matches the sam2consensus output "
              "names; the read columns will be empty.", file=sys.stderr)
    for key, path in result["written"].items():
        print(f"Wrote {path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="anchovy",
        description="Single-cell viral consensus and genotype-network pipeline.",
    )
    parser.add_argument("--version", action="version",
                        version=f"%(prog)s {__version__}")

    sub = parser.add_subparsers(dest="command", required=True,
                                metavar="{extract,fasta,consensus,annotate,frequencies}")

    # --- extract ---
    p_extract = sub.add_parser(
        "extract", help="Extract cell barcodes/UMIs from a mapped SAM/BAM.")
    p_extract.add_argument("sam", help="Path to input mapped SAM/BAM file.")
    p_extract.add_argument("whitelist", help="Path to 10X barcode whitelist.")
    p_extract.add_argument("-o", "--out", help="Output CSV path (default: alongside input).")
    p_extract.add_argument("--signature", help="10X signature (default: v2/v3 3').")
    p_extract.add_argument("--threads", type=int, help="Worker processes (default: 16).")
    p_extract.add_argument("--max-distance", type=int, dest="max_distance",
                           help="Max Levenshtein distance for the SIGNATURE match "
                                "to keep a read (default: 42).")
    p_extract.add_argument("--max-barcode-errors", type=int,
                           dest="max_barcode_errors",
                           help="Drop reads whose assigned barcode carries more "
                                "than this many errors. Unset, every read is "
                                "assigned to its nearest whitelist barcode with "
                                "no floor, however poorly it matched. 0 admits "
                                "only exact barcodes, 1 allows one substitution.")
    p_extract.set_defaults(func=_cmd_extract)

    # --- fasta ---
    p_fasta = sub.add_parser(
        "fasta", help="Write one FASTA per cell barcode from an anchovy CSV.")
    p_fasta.add_argument("csv", help="anchovy output CSV (from `anchovy extract`).")
    p_fasta.add_argument("outdir", help="Directory to write per-cell FASTAs into.")
    p_fasta.add_argument("--min-reads", type=int, dest="min_reads",
                         help="Minimum reads per cell to emit a FASTA (default: 5).")
    p_fasta.set_defaults(func=_cmd_fasta)

    # --- consensus ---
    p_cons = sub.add_parser(
        "consensus", help="Filter consensus sequences and summarize genotypes.")
    p_cons.add_argument("fasta", help="Path to <NAME>_allConsensus.fasta.")
    p_cons.add_argument("start", type=int, nargs="?", default=None,
                        help="ORF/region start (nt). With --whole-reference this "
                             "is the analysis window start instead, and optional.")
    p_cons.add_argument("end", type=int, nargs="?", default=None,
                        help="ORF/region end (nt). With --whole-reference this is "
                             "the analysis window end instead, and optional.")
    p_cons.add_argument("--whole-reference", dest="whole_reference",
                        action="store_true",
                        help="Keep the full reference instead of trimming to "
                             "start/end, so genotype positions are genome "
                             "coordinates. Required for region-aware annotation "
                             "(`anchovy annotate --gff`).")
    p_cons.add_argument("--keep-ambiguous", dest="keep_ambiguous",
                        action="store_true",
                        help="Keep positions whose base is an IUPAC ambiguity "
                             "code (R, Y, S...), N, or lowercase. By default "
                             "these are skipped: they mark positions where the "
                             "reads disagreed, not called mutations, and "
                             "treating them as variants groups cells by shared "
                             "uncertainty. Use this if you are after genuine "
                             "within-cell mixed populations.")
    p_cons.add_argument("--reference",
                        help="Reference sequence, as FASTA or a raw sequence file "
                             "(default: compute the consensus of the input).")
    p_cons.add_argument("--out-prefix", dest="out_prefix",
                        help="Output path prefix (default: derived from input).")
    p_cons.add_argument("--depth-min", type=int, dest="depth_min",
                        help="Minimum coverage to keep a sequence (default: 10).")
    p_cons.add_argument("--max-gaps", type=int, dest="max_gaps",
                        help="Max gaps allowed in region (default: 3).")
    p_cons.add_argument("--min-breadth", type=float, dest="min_breadth",
                        help="Minimum fraction of the reference a cell must "
                             "actually call (0-1) to be kept. Off by default. "
                             "Unlike --depth-min this measures only how much of "
                             "the genome the cell saw, which is what keeps "
                             "genotype strings comparable between cells.")
    p_cons.add_argument("--min-depth-called", type=float, dest="min_depth_called",
                        help="Minimum mean depth at the positions a cell "
                             "actually called. Off by default. This is the "
                             "quality half of --depth-min, with the breadth "
                             "half factored out, so a deep cell spanning less "
                             "of the genome is no longer penalised for its "
                             "uncovered flanks.")
    p_cons.set_defaults(func=_cmd_consensus)

    # --- annotate ---
    p_annot = sub.add_parser(
        "annotate", help="Annotate mutations and build genotype networks.")
    p_annot.add_argument("csv", help="filtConsensus.csv (from `anchovy consensus`).")
    p_annot.add_argument("reference", help="Reference sequence file.")
    p_annot.add_argument("out_prefix", help="Output path prefix.")
    p_annot.add_argument("--no-network", action="store_true",
                         help="Skip generating the network CSVs.")
    p_annot.add_argument("--self-edges", dest="self_edges", action="store_true",
                         help="Keep each genotype's edge to itself in the "
                              "epistatic network. Off by default: Cytoscape "
                              "draws them as a loop on every node and no node "
                              "is lost by removing them. Use this to reproduce "
                              "the original R output exactly.")
    p_annot.add_argument("--gff",
                         help="GFF3 of genome regions. Enables region-aware "
                              "annotation: mutation positions are treated as "
                              "genome coordinates and annotated against every "
                              "containing region, written to "
                              "<out_prefix>_regionAnnotations.csv. Requires the "
                              "consensus stage to have run --whole-reference.")
    p_annot.add_argument("--max-mutations", type=int, dest="max_mutations",
                         help="Drop cells carrying this many called mutations "
                              "or more, as likely artifacts (default: 200).")
    p_annot.set_defaults(func=_cmd_annotate)

    # --- frequencies ---
    p_freq = sub.add_parser(
        "frequencies",
        help="Allele frequencies over all mapped cells, with per-position "
             "denominators.")
    p_freq.add_argument("fasta", help="Path to <NAME>_allConsensus.fasta.")
    p_freq.add_argument("--reference", required=True,
                        help="Reference genome, as FASTA or a raw sequence file.")
    p_freq.add_argument("--out-prefix", dest="out_prefix", required=True,
                        help="Output path prefix.")
    p_freq.add_argument("--counts-dir", dest="counts_dir",
                        help="Directory of per-cell <ref>__<cell>_counts.tsv "
                             "files from `sam2consensus --counts`. Without it "
                             "the cell-vote columns are still produced and only "
                             "the read columns are left empty.")
    p_freq.add_argument("--reference-name", dest="reference_name", default="",
                        help="Reference name as it appears in the counts "
                             "filenames (the SAM's reference). Required to "
                             "locate them when --counts-dir is given.")
    p_freq.add_argument("--subset", dest="subset",
                        help="filtConsensus.csv naming the cells the genotype "
                             "network is built from. Their counts are reported "
                             "alongside the full population's, so the subset "
                             "can be checked for bias rather than assumed "
                             "representative.")
    p_freq.add_argument("--gff", dest="gff",
                        help="GFF3 for region names on each position.")
    p_freq.add_argument("--min-alt-reads", type=int, dest="min_alt_reads",
                        default=2,
                        help="Per-cell table: minimum reads carrying an allele "
                             "to report it (default: 2). Every covered position "
                             "carries singleton error reads, so 1 reports the "
                             "error spectrum rather than the variants.")
    p_freq.add_argument("--min-alt-freq", type=float, dest="min_alt_freq",
                        default=0.0,
                        help="Per-cell table: minimum within-cell allele "
                             "fraction to report (default: 0.0).")
    p_freq.add_argument("--min-cells", type=int, dest="min_cells", default=1,
                        help="Population table: minimum cells calling an allele "
                             "for it to get a row (default: 1).")
    p_freq.set_defaults(func=_cmd_frequencies)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
