#!/usr/bin/env python
"""
genbank_to_gff3.py -- turn a GenBank record into an anchovy region GFF3.

WHY THIS EXISTS, rather than a hand-written GFF3 shipped in the repo: the
coordinates would be copied by hand, and a coordinate typo in a region file does
not crash anything -- it silently renumbers every amino acid downstream. Deriving
them from the GenBank record the pipeline already downloads removes that failure
mode entirely, and keeps the example correct for whatever reference you point it
at rather than only the one it was written against.

It also gets you the mature peptides for free. Picornavirus records annotate the
polyprotein CDS plus a mat_peptide for each cleavage product (VP1-VP4, 2A-2C,
3A-3D), which is exactly the overlapping-region case anchovy annotates in long
format: a mutation in VP1 is numbered both within VP1 and within the polyprotein.

    python genbank_to_gff3.py record.gb regions.gff3

Features mapped to anchovy's vocabulary (see src/anchovy/regions.py):
    CDS          -> CDS                    (coding)
    mat_peptide  -> mature_protein_region  (coding)
    5'UTR        -> five_prime_UTR         (non-coding)
    3'UTR        -> three_prime_UTR        (non-coding)
    stem_loop    -> stem_loop              (non-coding)
Everything else in the record is skipped, which anchovy would ignore anyway.
"""

from __future__ import annotations

import sys

from Bio import SeqIO

# GenBank feature key -> the GFF3 type anchovy recognizes.
FEATURE_MAP = {
    "CDS": "CDS",
    "mat_peptide": "mature_protein_region",
    "5'UTR": "five_prime_UTR",
    "3'UTR": "three_prime_UTR",
    "stem_loop": "stem_loop",
}

# Qualifiers to try, in order, for a human-readable region name.
NAME_QUALIFIERS = ("gene", "product", "note", "locus_tag")


def _name_for(feature, fallback: str) -> str:
    for key in NAME_QUALIFIERS:
        values = feature.qualifiers.get(key)
        if values and values[0].strip():
            # GFF3 attributes are ';'-separated and '='-delimited, so a name
            # containing either would corrupt column 9.
            return values[0].strip().replace(";", ",").replace("=", "-")
    return fallback


def _phase_for(feature) -> str:
    """GFF3 phase from the GenBank codon_start qualifier.

    codon_start is 1-based ("start translating at the Nth base"); GFF3 phase is
    the number of bases to SKIP. They differ by one, which is an easy off-by-one
    to inherit silently.
    """
    codon_start = feature.qualifiers.get("codon_start", ["1"])[0]
    try:
        return str(int(codon_start) - 1)
    except ValueError:
        return "0"


def convert(record) -> list[str]:
    """GenBank record -> GFF3 lines (header included)."""
    seqid = record.id
    lines = ["##gff-version 3",
             f"##sequence-region {seqid} 1 {len(record.seq)}"]

    for feature in record.features:
        gff_type = FEATURE_MAP.get(feature.type)
        if gff_type is None:
            continue

        # Biopython locations are 0-based half-open; GFF3 is 1-based inclusive.
        start = int(feature.location.start) + 1
        end = int(feature.location.end)
        strand = "-" if feature.location.strand == -1 else "+"

        coding = gff_type in ("CDS", "mature_protein_region")
        phase = _phase_for(feature) if coding else "."
        name = _name_for(feature, f"{gff_type}_{start}")

        lines.append("\t".join([
            seqid, "genbank_to_gff3", gff_type, str(start), str(end),
            ".", strand, phase, f"Name={name}",
        ]))
    return lines


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(__doc__.strip(), file=sys.stderr)
        return 2
    gb_path, gff_path = argv[1], argv[2]

    # A failed efetch can return an HTML error page, or a truncated record, with
    # a 200 status. Biopython then raises somewhere deep in its scanner and the
    # user gets a stack trace instead of being told the download is bad.
    try:
        records = list(SeqIO.parse(gb_path, "genbank"))
    except Exception as exc:                      # noqa: BLE001 - reported, not swallowed
        first = ""
        try:
            with open(gb_path) as handle:
                first = handle.readline().strip()[:120]
        except OSError:
            pass
        print(f"error: {gb_path} is not a readable GenBank record ({exc}).\n"
              f"  The download probably failed or returned an error page.\n"
              f"  First line: {first!r}", file=sys.stderr)
        return 1

    if not records:
        print(f"error: no GenBank records in {gb_path} -- the download is empty "
              f"or is not GenBank format.", file=sys.stderr)
        return 1
    if len(records) > 1:
        print(f"error: {gb_path} holds {len(records)} records; anchovy assumes a "
              f"single reference. Split it or pass the one you mean.",
              file=sys.stderr)
        return 1

    lines = convert(records[0])
    feature_rows = [ln for ln in lines if not ln.startswith("#")]
    if not feature_rows:
        print(f"error: {gb_path} has no CDS/mat_peptide/UTR features to convert. "
              f"Region-aware annotation needs at least one.", file=sys.stderr)
        return 1

    with open(gff_path, "w") as handle:
        handle.write("\n".join(lines) + "\n")

    print(f"wrote {gff_path} ({len(feature_rows)} regions)")
    for row in feature_rows:
        cols = row.split("\t")
        print(f"    {cols[8][5:]:<24} {cols[2]:<22} {cols[3]}-{cols[4]} "
              f"strand {cols[6]} phase {cols[7]}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
