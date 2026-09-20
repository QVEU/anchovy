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
    )

    reference = None
    if args.reference:
        reference = open(args.reference).read().strip()

    result = consensus.run(
        fasta=args.fasta, start=args.start, end=args.end,
        reference=reference, config=config, out_prefix=args.out_prefix,
    )
    written = result["written"]
    print(f"Kept {len(result['records'])} sequences.")
    print(f"Wrote {written['reference']}")
    print(f"Wrote {written['csv']}")
    return 0


def _cmd_annotate(args: argparse.Namespace) -> int:
    from anchovy import annotate

    result = annotate.run(
        filt_consensus_csv=args.csv,
        reference_file=args.reference,
        out_prefix=args.out_prefix,
        network=not args.no_network,
        gff=args.gff,          # None -> legacy frame-1; set -> region-aware
    )
    for label, path in result["written"].items():
        print(f"Wrote {path}")
    return 0


# --------------------------------------------------------------------------- #
# parser
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="anchovy",
        description="Single-cell viral consensus and genotype-network pipeline.",
    )
    parser.add_argument("--version", action="version",
                        version=f"%(prog)s {__version__}")

    sub = parser.add_subparsers(dest="command", required=True,
                                metavar="{extract,fasta,consensus}")

    # --- extract ---
    p_extract = sub.add_parser(
        "extract", help="Extract cell barcodes/UMIs from a mapped SAM/BAM.")
    p_extract.add_argument("sam", help="Path to input mapped SAM/BAM file.")
    p_extract.add_argument("whitelist", help="Path to 10X barcode whitelist.")
    p_extract.add_argument("-o", "--out", help="Output CSV path (default: alongside input).")
    p_extract.add_argument("--signature", help="10X signature (default: v2/v3 3').")
    p_extract.add_argument("--threads", type=int, help="Worker processes (default: 16).")
    p_extract.add_argument("--max-distance", type=int, dest="max_distance",
                           help="Max Levenshtein distance to keep a read (default: 42).")
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
                        help="Optional analysis window start (1-based, nt). "
                             "Omit to analyze the whole reference.")
    p_cons.add_argument("end", type=int, nargs="?", default=None,
                        help="Optional analysis window end (1-based, nt).")
    p_cons.add_argument("--reference", help="Reference sequence file (default: compute consensus).")
    p_cons.add_argument("--out-prefix", dest="out_prefix",
                        help="Output path prefix (default: derived from input).")
    p_cons.add_argument("--depth-min", type=int, dest="depth_min",
                        help="Minimum coverage to keep a sequence (default: 10).")
    p_cons.add_argument("--max-gaps", type=int, dest="max_gaps",
                        help="Max gaps allowed in window (default: 3).")
    p_cons.set_defaults(func=_cmd_consensus)

    # --- annotate ---
    p_annot = sub.add_parser(
        "annotate", help="Annotate mutations and build genotype networks.")
    p_annot.add_argument("csv", help="filtConsensus.csv (from `anchovy consensus`).")
    p_annot.add_argument("reference", help="Reference sequence file.")
    p_annot.add_argument("out_prefix", help="Output path prefix.")
    p_annot.add_argument("--gff", default=None,
                         help="GFF3 file of regions. If given, annotation is "
                              "region-aware (coding + non-coding, both strands) "
                              "and a _regionAnnotations.csv is written. If omitted, "
                              "legacy single-reference frame-1 annotation is used.")
    p_annot.add_argument("--no-network", action="store_true",
                         help="Skip generating the network CSVs.")
    p_annot.set_defaults(func=_cmd_annotate)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
