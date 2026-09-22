"""
regions.py -- region model + GFF3 parsing + region-aware mutation annotation.

This is the foundation for annotating mutations against a genome that has both
coding and non-coding regions, and multiple (possibly overlapping) features. It
replaces the old single-ORF assumption (a bare start/end with implicit frame-1
translation) with a general model driven entirely by a GFF3 file.

KEY DESIGN PRINCIPLES (all settled + verified before this was written):
  - The GFF3 is the source of truth. Strand, phase, coordinates come from the
    file; nothing is inferred from feature nesting.
  - Long-format output: a mutation is annotated against EVERY region containing
    it. A variant inside NS5 (which is inside the polyprotein) yields two rows --
    one per region -- each numbered in that region's own frame of reference.
  - The nucleotide identity is always genome/forward-strand:
    mutation_id = "{regionTag}:{wtBase}{genomePos}{mutBase}".
    Amino-acid columns follow the region's strand/frame (reverse-complemented for
    minus-strand features); the base columns stay in genome/forward terms.

CODON MATH (verified numerically against real Dengue coords and a synthetic
minus-strand gene before coding -- see the project discussion):
  + strand: coding frame begins at (start + phase), residues count upward
  - strand: coding frame begins at (end - phase),   residues count downward,
            bases complemented, feature span reverse-complemented for translation

Coordinates are 1-based inclusive throughout (GFF3 convention).
"""

from __future__ import annotations

from dataclasses import dataclass

from Bio.Seq import Seq


# GFF3 feature types we treat as coding vs non-coding. Others are ignored on
# parse (genes, mRNAs, exons, etc.) -- we annotate against CDS/mature proteins
# and UTR-like non-coding features.
CODING_TYPES = {"CDS", "mature_protein_region", "mature_protein_region_of_CDS"}
NONCODING_TYPES = {
    "five_prime_UTR", "three_prime_UTR", "UTR",
    "five_prime_utr", "three_prime_utr",
    "stem_loop", "ncRNA", "misc_feature", "region",
}

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass
class Region:
    """One annotated region from the GFF3.

    Attributes:
        name: display name / tag (from Name/ID/gene attribute, else generated).
        start, end: 1-based inclusive genome coordinates.
        strand: '+' or '-'.
        phase: 0/1/2 GFF3 phase for coding features (bases to skip to the first
            complete codon, in translation direction). Ignored for non-coding.
        coding: True if this region's mutations get amino-acid annotation.
    """
    name: str
    start: int
    end: int
    strand: str
    phase: int
    coding: bool

    def contains(self, genome_pos: int) -> bool:
        return self.start <= genome_pos <= self.end


# --------------------------------------------------------------------------- #
# GFF3 parsing
# --------------------------------------------------------------------------- #
def _parse_attributes(field: str) -> dict[str, str]:
    """Parse the GFF3 column-9 attribute string into a dict.

    GFF3 attributes are 'key=value;key=value'. We're lenient about whitespace
    and ignore malformed pairs rather than failing the whole parse.
    """
    attrs: dict[str, str] = {}
    for part in field.strip().split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        key, _, value = part.partition("=")
        attrs[key.strip()] = value.strip()
    return attrs


def _region_name(attrs: dict[str, str], feature_type: str, start: int) -> str:
    """Pick a name from attributes, in priority order, else generate one."""
    for key in ("Name", "name", "ID", "gene", "product"):
        if key in attrs and attrs[key]:
            return attrs[key]
    return f"{feature_type}_{start}"


def parse_gff3(path: str) -> list[Region]:
    """Parse a GFF3 file into a list of Region objects (coding + non-coding).

    Supports a pragmatic subset of GFF3: standard 9 tab-separated columns, with
    the feature types in CODING_TYPES / NONCODING_TYPES. Comment lines (starting
    with '#') and blank lines are skipped; unrecognized feature types are ignored.

    Raises ValueError on structurally broken rows (wrong column count, bad
    coordinates) so problems surface loudly rather than silently mis-annotating.
    """
    regions: list[Region] = []
    with open(path) as handle:
        for lineno, raw in enumerate(handle, 1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 8:
                raise ValueError(
                    f"{path}:{lineno}: expected >=8 tab-separated GFF3 columns, "
                    f"got {len(cols)}")

            seqid, source, ftype, start_s, end_s, score, strand, phase_s = cols[:8]
            attrs = _parse_attributes(cols[8]) if len(cols) > 8 else {}

            ftype_norm = ftype.strip()
            is_coding = ftype_norm in CODING_TYPES
            is_noncoding = ftype_norm in NONCODING_TYPES
            if not (is_coding or is_noncoding):
                continue  # ignore feature types we don't annotate against

            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                raise ValueError(
                    f"{path}:{lineno}: non-integer coordinates "
                    f"'{start_s}'/'{end_s}'")
            if start > end:
                raise ValueError(
                    f"{path}:{lineno}: start {start} > end {end}")
            if strand not in ("+", "-"):
                raise ValueError(
                    f"{path}:{lineno}: strand must be '+' or '-', got '{strand}'")

            # Phase only meaningful for coding; default 0 and tolerate '.'.
            if is_coding:
                phase = 0 if phase_s.strip() in (".", "") else int(phase_s)
                if phase not in (0, 1, 2):
                    raise ValueError(
                        f"{path}:{lineno}: phase must be 0/1/2, got '{phase_s}'")
            else:
                phase = 0

            regions.append(Region(
                name=_region_name(attrs, ftype_norm, start),
                start=start, end=end, strand=strand, phase=phase,
                coding=is_coding,
            ))
    return regions


def cds_window(path: str) -> tuple[int, int]:
    """The CDS feature's bounds, as a 0-based half-open analysis window.

    WHY THIS IS NOT LEFT TO THE USER. The window (orf_start/orf_end) is 0-based
    half-open because it indexes the consensus string directly, while GFF3 is
    1-based inclusive. Copying CDS bounds from an annotation into a config by
    hand therefore needs an off-by-one applied to one end and not the other --
    and getting it wrong does not fail, it silently shifts which positions are
    called and which cells are judged complete. That is the same class of error
    genbank_to_gff3.py exists to prevent, so the conversion is done once, here.

    Takes the CDS specifically rather than the union of coding regions: the
    mature peptides parse as coding too, and a record missing some of them
    would quietly yield a shorter window than the coding sequence.

    Raises ValueError if the GFF has no CDS, or more than one -- a multi-CDS
    record has no single coding sequence to take a window from, and guessing
    between them would be a coordinate error of exactly the kind above.
    """
    spans: list[tuple[int, int]] = []
    with open(path) as handle:
        for lineno, raw in enumerate(handle, 1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 5 or cols[2].strip() != "CDS":
                continue
            try:
                spans.append((int(cols[3]), int(cols[4])))
            except ValueError:
                raise ValueError(
                    f"{path}:{lineno}: CDS row has non-integer coordinates "
                    f"({cols[3]!r}, {cols[4]!r})") from None

    if not spans:
        raise ValueError(
            f"{path} has no CDS feature, so there is no coding sequence to take "
            f"an analysis window from. Give explicit orf_start/orf_end instead, "
            f"or use a GFF3 that annotates the CDS.")
    if len(spans) > 1:
        listed = ", ".join(f"{a}-{b}" for a, b in sorted(spans))
        raise ValueError(
            f"{path} has {len(spans)} CDS features ({listed}). A window must "
            f"name one coding sequence; with several there is no unambiguous "
            f"choice, so set orf_start/orf_end explicitly for the one you mean.")

    start, end = spans[0]
    if start < 1 or end < start:
        raise ValueError(
            f"{path}: CDS coordinates {start}-{end} are not a valid 1-based "
            f"inclusive span.")
    # 1-based inclusive -> 0-based half-open.
    return start - 1, end


def validate_regions(regions: list[Region]) -> list[str]:
    """Return a list of warning strings for suspicious region definitions.

    These are the cheap sanity checks that catch common coordinate mistakes
    (e.g. specifying a protein in residue units instead of nucleotides): a coding
    region whose length isn't a multiple of 3 after removing the phase offset is
    almost always a coordinate error.
    """
    warnings: list[str] = []
    for r in regions:
        if r.coding:
            coding_len = (r.end - r.start + 1) - r.phase
            if coding_len % 3 != 0:
                warnings.append(
                    f"coding region '{r.name}' ({r.start}-{r.end}, phase {r.phase}) "
                    f"has coding length {coding_len}, not a multiple of 3 -- "
                    f"check the coordinates.")
    return warnings


# --------------------------------------------------------------------------- #
# Per-region annotation
# --------------------------------------------------------------------------- #
def classify(regions: list[Region], genome_pos: int) -> list[Region]:
    """All regions containing a genome position (order preserved)."""
    return [r for r in regions if r.contains(genome_pos)]


def _coding_frame_pos(region: Region, genome_pos: int) -> int:
    """1-based position within the region's coding frame (translation direction).

    <= 0 means the position falls in the phase-skipped bases before the first
    complete codon (not translatable within this feature).
    """
    if region.strand == "+":
        return genome_pos - (region.start + region.phase) + 1
    else:
        return (region.end - region.phase) - genome_pos + 1


def annotate_in_region(region: Region, genome_pos: int, wt_base: str,
                       mut_base: str, reference: str) -> dict:
    """Annotate one mutation within one region.

    Args:
        region: the containing Region.
        genome_pos: 1-based genome position of the mutation.
        wt_base, mut_base: reference and mutant nucleotide, GENOME/forward strand.
        reference: full reference sequence (1-based externally; str is 0-based).

    Returns a dict of annotation columns. Amino-acid fields are populated only
    for coding regions where the position sits in a complete codon; otherwise
    they are None. The mutation_id is always genome/forward-strand.
    """
    wt_base = wt_base.upper()
    mut_base = mut_base.upper()
    row = {
        "mutation_id": f"{region.name}:{wt_base}{genome_pos}{mut_base}",
        "region": region.name,
        "region_type": "coding" if region.coding else "non-coding",
        "strand": region.strand,
        "genome_pos": genome_pos,
        "wt_base": wt_base,
        "mut_base": mut_base,
        "wt_aa": None, "mut_aa": None, "residue": None,
        "codon_wt": None, "codon_mut": None, "sub_class": None,
    }
    if not region.coding:
        return row

    frame_pos = _coding_frame_pos(region, genome_pos)
    if frame_pos <= 0:
        # In the phase-skipped lead-in; not a complete codon here.
        row["sub_class"] = "non-coding-frame"
        return row

    residue = (frame_pos - 1) // 3 + 1
    pos_in_codon = (frame_pos - 1) % 3
    row["residue"] = residue

    # Extract the WT codon's three genome positions (in translation direction),
    # then build the codon sequence, complementing for minus strand.
    if region.strand == "+":
        codon_genome_start = genome_pos - pos_in_codon          # 1-based
        codon_nts = reference[codon_genome_start - 1: codon_genome_start + 2]
        codon_wt = codon_nts.upper()
        # apply mutation at pos_in_codon within the codon
        cm = list(codon_wt)
        cm[pos_in_codon] = mut_base
        codon_mut = "".join(cm)
    else:
        # On minus strand the codon runs downward in genome coords; residue's
        # first base (translation 5') is at the higher genome coordinate.
        codon_genome_high = genome_pos + pos_in_codon           # 5' base (genome)
        span = reference[codon_genome_high - 3: codon_genome_high]  # 3 genome nts
        codon_wt = _revcomp(span.upper())
        # mutation: complement the genome base; its position in the (revcomp)
        # codon is pos_in_codon from the 5' (high-genome) end.
        cm = list(codon_wt)
        cm[pos_in_codon] = _COMPLEMENT_CHAR(mut_base)
        codon_mut = "".join(cm)

    row["codon_wt"] = codon_wt
    row["codon_mut"] = codon_mut
    # Translate the single codons (guard against non-ACGT / partial).
    row["wt_aa"] = _safe_translate_codon(codon_wt)
    row["mut_aa"] = _safe_translate_codon(codon_mut)
    if row["wt_aa"] and row["mut_aa"] and "X" not in (row["wt_aa"], row["mut_aa"]):
        row["sub_class"] = "Syn" if row["wt_aa"] == row["mut_aa"] else "Non-Syn"
    else:
        row["sub_class"] = "X"
    return row


def _COMPLEMENT_CHAR(base: str) -> str:
    return base.upper().translate(_COMPLEMENT)


def _safe_translate_codon(codon: str) -> str:
    """Translate a 3-nt codon to one amino acid; 'X' if not translatable."""
    if len(codon) != 3 or any(c not in "ACGT" for c in codon):
        return "X"
    return str(Seq(codon).translate())


def annotate_mutation(regions: list[Region], genome_pos: int, wt_base: str,
                      mut_base: str, reference: str) -> list[dict]:
    """Annotate a mutation against every region containing it (long format).

    If the position is in no region, returns a single 'intergenic' row so every
    mutation gets at least one genome-anchored record.
    """
    containing = classify(regions, genome_pos)
    if not containing:
        wt_base, mut_base = wt_base.upper(), mut_base.upper()
        return [{
            "mutation_id": f"intergenic:{wt_base}{genome_pos}{mut_base}",
            "region": "intergenic", "region_type": "non-coding", "strand": ".",
            "genome_pos": genome_pos, "wt_base": wt_base, "mut_base": mut_base,
            "wt_aa": None, "mut_aa": None, "residue": None,
            "codon_wt": None, "codon_mut": None, "sub_class": None,
        }]
    return [annotate_in_region(r, genome_pos, wt_base, mut_base, reference)
            for r in containing]
