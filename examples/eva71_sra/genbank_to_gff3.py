#!/usr/bin/env python
"""
genbank_to_gff3.py -- turn a GenBank record into an anchovy region GFF3.

WHY THIS EXISTS, rather than a hand-written GFF3 shipped in the repo: the
coordinates would be copied by hand, and a coordinate typo in a region file does
not crash anything -- it silently renumbers every amino acid downstream. Deriving
them from the GenBank record the pipeline already downloads removes that failure
mode entirely, and keeps the example correct for whatever reference you point it
at rather than only the one it was written against.

MATURE PEPTIDES ARE NOT GUARANTEED. Many picornavirus records annotate the
polyprotein CDS plus a mat_peptide per cleavage product (VP1-VP4, 2A-2C, 3A-3D),
which is the overlapping-region case anchovy annotates in long format: a VP1
mutation numbered both within VP1 and within the polyprotein. Plenty of records
carry only the CDS -- AF304458 (EV-A71 Tainan/4643/98) is one -- and then every
mutation is reported against the 2194-residue polyprotein and nothing else.

    python genbank_to_gff3.py record.gb regions.gff3
    python genbank_to_gff3.py record.gb regions.gff3 --transfer-from donor.gb

--transfer-from fills that gap WITHOUT hand-typing coordinates, which is the
failure mode this script exists to avoid. It translates both polyproteins,
aligns them, and carries each donor mat_peptide through the alignment onto the
target. Boundaries land where the sequences actually correspond rather than at
an assumed offset, and the alignment identity is reported and enforced, so a
donor that is not the same virus is refused instead of silently renumbering
everything.

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
from Bio.Align import PairwiseAligner, substitution_matrices

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


# Below this, the donor is not the same virus and its cleavage sites cannot be
# trusted to correspond. Two EV-A71 genotypes sit well above 90% at the protein
# level, so this refuses a mistake rather than rejecting legitimate donors.
MIN_TRANSFER_IDENTITY = 0.80


def _polyprotein_cds(record):
    """The record's single CDS, or None. Multiple CDSs are ambiguous here."""
    cds = [f for f in record.features if f.type == "CDS"]
    return cds[0] if len(cds) == 1 else None


def _aa_offsets(feature, cds) -> tuple[int, int]:
    """A mat_peptide's residue span, as offsets from the start of the CDS."""
    start_nt = int(feature.location.start) - int(cds.location.start)
    end_nt = int(feature.location.end) - int(cds.location.start)
    return start_nt // 3, end_nt // 3


def transfer_mat_peptides(donor, target,
                          min_identity: float = MIN_TRANSFER_IDENTITY):
    """Carry the donor's mat_peptides onto the target through a protein alignment.

    Returns (features, report). Each feature is a dict with name/start/end in
    1-based inclusive TARGET genome coordinates. `report` carries the identity
    and any per-peptide notes, so the caller can print what it did rather than
    the user having to trust it.

    Raises ValueError when the transfer cannot be justified -- either record
    lacking a single CDS, the donor carrying no mat_peptides, or the two
    polyproteins aligning below `min_identity`.
    """
    donor_cds, target_cds = _polyprotein_cds(donor), _polyprotein_cds(target)
    if donor_cds is None:
        raise ValueError("donor record does not have exactly one CDS")
    if target_cds is None:
        raise ValueError("target record does not have exactly one CDS")

    peptides = [f for f in donor.features if f.type == "mat_peptide"]
    if not peptides:
        raise ValueError("donor record has no mat_peptide features to transfer")

    donor_aa = str(donor_cds.translate(donor.seq, cds=False)).rstrip("*")
    target_aa = str(target_cds.translate(target.seq, cds=False)).rstrip("*")

    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "global"
    alignment = aligner.align(donor_aa, target_aa)[0]

    # Residue-level map, donor index -> target index, built from the aligned
    # blocks. Positions the alignment does not pair simply do not appear, which
    # is what makes an unmappable boundary detectable rather than silently
    # approximated.
    donor_to_target: dict[int, int] = {}
    matches = 0
    for (d_start, d_end), (t_start, t_end) in zip(*alignment.aligned):
        for offset in range(d_end - d_start):
            d_i, t_i = d_start + offset, t_start + offset
            donor_to_target[d_i] = t_i
            if donor_aa[d_i] == target_aa[t_i]:
                matches += 1

    identity = matches / max(len(donor_aa), len(target_aa))
    if identity < min_identity:
        raise ValueError(
            f"donor and target polyproteins are only {identity:.1%} identical "
            f"(need {min_identity:.0%}). Transferring cleavage sites between "
            f"sequences this different would put every boundary in the wrong "
            f"place. Is the donor the same virus?")

    target_cds_start = int(target_cds.location.start)
    strand = "-" if target_cds.location.strand == -1 else "+"

    features, notes = [], []
    for peptide in peptides:
        name = _name_for(peptide, "mat_peptide")
        d_start_aa, d_end_aa = _aa_offsets(peptide, donor_cds)

        # A peptide's last residue is d_end_aa - 1; map both ends inclusively.
        t_start_aa = donor_to_target.get(d_start_aa)
        t_last_aa = donor_to_target.get(d_end_aa - 1)
        if t_start_aa is None or t_last_aa is None:
            notes.append(f"{name}: boundary fell in an alignment gap, skipped")
            continue

        start = target_cds_start + t_start_aa * 3 + 1          # 1-based
        end = target_cds_start + (t_last_aa + 1) * 3           # inclusive
        features.append({"name": name, "start": start, "end": end,
                         "strand": strand,
                         "shifted": (d_start_aa != t_start_aa)})

    return features, {"identity": identity, "notes": notes,
                      "donor_len": len(donor_aa), "target_len": len(target_aa)}


def convert(record, transferred=None) -> list[str]:
    """GenBank record -> GFF3 lines (header included).

    `transferred` is an optional list of mat_peptide features carried over from
    a donor record by transfer_mat_peptides(); they are emitted after the
    record's own features, marked in the GFF3 so their origin is visible in the
    file rather than only in the terminal that produced it.
    """
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

    for feature in transferred or []:
        lines.append("\t".join([
            seqid, "genbank_to_gff3", "mature_protein_region",
            str(feature["start"]), str(feature["end"]),
            ".", feature["strand"], "0",
            f"Name={feature['name']};Note=transferred_by_alignment",
        ]))
    return lines


def _load_single_record(path: str):
    """Parse a GenBank file that must hold exactly one record."""
    records = list(SeqIO.parse(path, "genbank"))
    if not records:
        raise ValueError(f"no GenBank records in {path}")
    if len(records) > 1:
        raise ValueError(f"{path} holds {len(records)} records; expected one")
    return records[0]


def main(argv: list[str]) -> int:
    args = list(argv[1:])
    donor_path = None
    if "--transfer-from" in args:
        i = args.index("--transfer-from")
        if i + 1 >= len(args):
            print("error: --transfer-from needs a GenBank file", file=sys.stderr)
            return 2
        donor_path = args[i + 1]
        del args[i:i + 2]

    if len(args) != 2:
        print(__doc__.strip(), file=sys.stderr)
        return 2
    gb_path, gff_path = args

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

    record = records[0]

    transferred = None
    if donor_path:
        try:
            donor = _load_single_record(donor_path)
            transferred, report = transfer_mat_peptides(donor, record)
        except (ValueError, OSError) as exc:
            print(f"error: could not transfer mature peptides from "
                  f"{donor_path}: {exc}", file=sys.stderr)
            return 1
        print(f"transferred {len(transferred)} mature peptide(s) from "
              f"{donor.id} ({report['identity']:.1%} identical over "
              f"{report['donor_len']} vs {report['target_len']} residues)")
        for note in report["notes"]:
            print(f"    warning: {note}")

    lines = convert(record, transferred=transferred)
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
        name = cols[8][5:].split(";")[0]
        origin = " (transferred)" if "transferred_by_alignment" in cols[8] else ""
        print(f"    {name:<24} {cols[2]:<22} {cols[3]}-{cols[4]} "
              f"strand {cols[6]} phase {cols[7]}{origin}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
