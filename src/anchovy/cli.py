"""
cli.py -- command-line entry point for the anchovy pipeline.

This is wired to pyproject.toml's [project.scripts] as `anchovy = "anchovy.cli:main"`,
so after `pip install -e .` the shell command `anchovy` calls main() here.

Right now the subcommands are stubs that just report they were reached. That's
deliberate: we stand up an installable, runnable skeleton FIRST, confirm the
plumbing works end to end, then replace each stub body with real logic (importing
from extract.py, fasta.py, consensus.py). Building the skeleton before the flesh
means every later change is verified against something that already runs.
"""

from __future__ import annotations

import argparse
import sys

from anchovy import __version__


def _cmd_extract(args: argparse.Namespace) -> int:
    # Will call extract.py (the former anchovy.py core).
    print(f"[anchovy extract] sam={args.sam} whitelist={args.whitelist}")
    print("  (not implemented yet -- skeleton stub)")
    return 0


def _cmd_fasta(args: argparse.Namespace) -> int:
    # Will call fasta.py (the former CBCtoFasta.py).
    print(f"[anchovy fasta] indir={args.indir} csv={args.csv}")
    print("  (not implemented yet -- skeleton stub)")
    return 0


def _cmd_consensus(args: argparse.Namespace) -> int:
    # Will call consensus.py (the former ConsensusTool.py).
    print(f"[anchovy consensus] fasta={args.fasta} start={args.start} end={args.end}")
    print("  (not implemented yet -- skeleton stub)")
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Construct the top-level parser and its subcommands.

    Keeping parser construction in its own function (rather than inline in main)
    makes it testable: a test can build the parser and assert on how it parses
    example argument lists, without actually running any command.
    """
    parser = argparse.ArgumentParser(
        prog="anchovy",
        description="Single-cell viral consensus and genotype-network pipeline.",
    )
    parser.add_argument("--version", action="version",
                        version=f"%(prog)s {__version__}")

    sub = parser.add_subparsers(dest="command", required=True,
                                metavar="{extract,fasta,consensus}")

    # anchovy extract: SAM -> per-read CBC/UMI table
    p_extract = sub.add_parser("extract",
                               help="Extract cell barcodes/UMIs from a mapped SAM.")
    p_extract.add_argument("sam", help="Path to input mapped SAM file.")
    p_extract.add_argument("whitelist", help="Path to 10X barcode whitelist.")
    p_extract.set_defaults(func=_cmd_extract)

    # anchovy fasta: per-read table -> per-cell FASTAs
    p_fasta = sub.add_parser("fasta",
                             help="Write one FASTA per cell barcode from anchovy CSV.")
    p_fasta.add_argument("indir", help="Directory containing the anchovy CSV / for output.")
    p_fasta.add_argument("csv", help="anchovy output CSV filename.")
    p_fasta.set_defaults(func=_cmd_fasta)

    # anchovy consensus: filter + genotype summary
    p_cons = sub.add_parser("consensus",
                            help="Filter consensus sequences and summarize genotypes.")
    p_cons.add_argument("fasta", help="Path to <NAME>_allConsensus.fasta.")
    p_cons.add_argument("start", type=int, help="ORF/region start (nt).")
    p_cons.add_argument("end", type=int, help="ORF/region end (nt).")
    p_cons.set_defaults(func=_cmd_consensus)

    return parser


def main(argv: list[str] | None = None) -> int:
    """Entry point. Returns an exit code so it's testable and shell-friendly.

    argv defaults to None so argparse reads sys.argv; passing a list lets tests
    drive it directly, e.g. main(["extract", "in.sam", "wl.txt"]).
    """
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
