"""
make_mapping_fixtures.py -- fixture for the mapping-onward stages (map -> consensus).

Run once:
    python tests/make_mapping_fixtures.py

WHY A SEPARATE FIXTURE
----------------------
The barcode fixture (make_fixtures.py) has reads with a 10X signature spliced
into random sequence -- great for testing barcode extraction, useless for
mapping (random reads don't align to any reference). The mapping/consensus
stages need reads that actually align to the template so minimap2 and
sam2consensus have something to work with.

Rather than contort one fixture to do both jobs (which would force regenerating
the already-frozen extract/fasta goldens), we keep two fixtures testing the two
halves of the pipeline. The DAG wiring that connects them is verified by dry-run.

THE SCENARIO (predicted so we can VERIFY, not just capture)
-----------------------------------------------------------
- A 600 nt deterministic template (tests/data/mapping/template.fasta).
- Per-cell FASTAs whose reads are 500 nt windows [50:550) of the template, so
  they seed and align well under minimap2 map-hifi.
- 6 reads per cell (above sam2consensus min-depth 5).
- Cell "cellA": reads match the template exactly -> consensus == template[50:550].
- Cell "cellB": every read carries ONE planted substitution at absolute template
  position 200 (T -> A) -> consensus == template[50:550] with index 150 changed.
  That position falls inside the CDS, so it exercises the CODING path.
- Cell "cellC": every read carries ONE planted substitution at absolute template
  position 120, which falls inside the 5'UTR, exercising the NON-CODING path --
  a mutation that must be reported with no amino-acid columns rather than
  silently translated as though it were coding.

WHY THREE CELLS, NOT TWO
------------------------
With only two cells the variant column split 50/50, and the computed per-column
consensus broke that tie by np.unique ordering -- alphabetically. The reference
base came out arbitrary, so the cell MATCHING the template was reported as the
mutant and the direction of every call was inverted. A third cell makes the
majority at each variant column the template base, which is both the realistic
case and the one whose output can be read at a glance. The tie-breaking behavior
itself is still covered, at the unit level, in tests/test_consensus.py.

These predictions are recorded in mapping_expected.json so the golden-freeze
step can verify sam2consensus output against them instead of trusting it blindly.
"""

from __future__ import annotations

import json
from pathlib import Path

DATA = Path(__file__).parent / "data" / "mapping"

# Deterministic template via a simple LCG, so the fixture is identical every run.
def _make_template(n: int, seed: int = 42) -> str:
    nt = "ACGT"
    s = seed & 0xFFFFFFFF
    out = []
    for _ in range(n):
        s = (s + 0x6D2B79F5) & 0xFFFFFFFF
        t = (s ^ (s >> 15)) * (1 | s) & 0xFFFFFFFF
        t = (t + ((t ^ (t >> 7)) * (61 | t) & 0xFFFFFFFF)) ^ t & 0xFFFFFFFF
        val = ((t ^ (t >> 14)) & 0xFFFFFFFF) / 4294967296
        out.append(nt[int(val * 4) % 4])
    return "".join(out)


TEMPLATE_LEN = 600
CORE_START = 50
CORE_END = 550
READS_PER_CELL = 6
VARIANT_ABS_POS = 200            # 0-based template index of the CODING variant
UTR_VARIANT_ABS_POS = 120        # 0-based template index of the NON-CODING variant
NT = "ACGT"

# The 5'UTR variant has to satisfy three overlapping constraints at once, which
# leaves a narrow legal band of 0-based [100, 149) -- genome 101..149:
#   * inside the 5'UTR      (genome 1-149, i.e. below CDS_START)
#   * covered by the reads  (genome 51-550, i.e. within [CORE_START, CORE_END))
#   * inside the analysis window used by workflow/config_test.yaml (orf_start 100
#     / orf_end 500, which are 0-based slice indices)
# Index 120 sits clear of all three edges.

# --- Region model for the fixture (regions.gff3) --------------------------- #
# 1-based inclusive genome coordinates, GFF3 convention. Deliberately chosen so
# the fixture exercises the thing the legacy annotator got WRONG: the CDS does
# NOT start at genome position 1, so correct residue numbering is only possible
# if the annotator reads the frame from the GFF instead of assuming frame 1.
#
#   5'UTR   1-149     non-coding
#   CDS   150-500     coding, + strand, phase 0  -> 351 nt = 117 codons
#   3'UTR 501-600     non-coding
#
# The planted variant sits at genome 201 (VARIANT_ABS_POS + 1), which is inside
# the CDS at codon frame position 201-150+1 = 52 -> residue 18, first base of
# its codon. Hand-verified before this fixture was written.
CDS_START = 150
CDS_END = 500
REGION_ROWS = [
    ("five_prime_UTR", 1,         CDS_START - 1, ".", "ID=5UTR;Name=5UTR"),
    ("CDS",            CDS_START, CDS_END,       "0", "ID=poly;Name=polyprotein"),
    ("three_prime_UTR", CDS_END + 1, TEMPLATE_LEN, ".", "ID=3UTR;Name=3UTR"),
]


def build():
    DATA.mkdir(parents=True, exist_ok=True)
    template = _make_template(TEMPLATE_LEN)

    # Write the template FASTA minimap2 maps against.
    (DATA / "template.fasta").write_text(f">testref\n{template}\n")

    core = template[CORE_START:CORE_END]

    # Variant base: something different from the reference base at that position.
    ref_base = template[VARIANT_ABS_POS]
    variant_base = NT[(NT.index(ref_base) + 1) % 4]
    variant_core_idx = VARIANT_ABS_POS - CORE_START

    # cellA: exact copies of the core.
    cellA_reads = [core for _ in range(READS_PER_CELL)]

    # cellB: every read carries the planted CODING variant.
    cb = list(core)
    cb[variant_core_idx] = variant_base
    cellB_reads = ["".join(cb) for _ in range(READS_PER_CELL)]

    # cellC: every read carries the planted NON-CODING (5'UTR) variant.
    utr_ref_base = template[UTR_VARIANT_ABS_POS]
    utr_variant_base = NT[(NT.index(utr_ref_base) + 1) % 4]
    utr_core_idx = UTR_VARIANT_ABS_POS - CORE_START
    cc = list(core)
    cc[utr_core_idx] = utr_variant_base
    cellC_reads = ["".join(cc) for _ in range(READS_PER_CELL)]

    def write_cell(name, reads):
        recs = [f">{name}_read{i}\n{seq}" for i, seq in enumerate(reads)]
        (DATA / f"{name}.fa").write_text("\n".join(recs) + "\n")

    write_cell("cellA", cellA_reads)
    write_cell("cellB", cellB_reads)
    write_cell("cellC", cellC_reads)

    # --- Predicted CONSENSUS output (verified against sam2consensus) --------- #
    # sam2consensus emits a FULL-REFERENCE-LENGTH consensus (TEMPLATE_LEN), with
    # positions lacking read coverage filled with the gap char '-'. Our reads
    # cover only [CORE_START, CORE_END), so positions outside that are gaps.
    # (This was learned empirically: an earlier prediction of just the bare core
    # was wrong about framing -- the tool gap-fills uncovered flanks.)
    gap = "-"
    lead = gap * CORE_START                       # uncovered 5' flank
    tail = gap * (TEMPLATE_LEN - CORE_END)        # uncovered 3' flank

    cellA_consensus = lead + core + tail
    cellB_consensus = lead + "".join(cb) + tail
    cellC_consensus = lead + "".join(cc) + tail

    # --- Region annotations (regions.gff3) ---------------------------------- #
    gff_lines = ["##gff-version 3"]
    for ftype, gstart, gend, phase, attrs in REGION_ROWS:
        gff_lines.append(
            "\t".join(["testref", "anchovy", ftype, str(gstart), str(gend),
                        ".", "+", phase, attrs]))
    (DATA / "regions.gff3").write_text("\n".join(gff_lines) + "\n")

    expected = {
        "template_ref_name": "testref",
        "template_len": TEMPLATE_LEN,
        "core_start": CORE_START,
        "core_end": CORE_END,
        "variant_abs_pos": VARIANT_ABS_POS,
        "variant_ref_base": ref_base,
        "variant_alt_base": variant_base,
        # Full-length gap-filled consensuses as sam2consensus actually produces.
        "consensus": {
            "cellA": cellA_consensus,
            "cellB": cellB_consensus,
            "cellC": cellC_consensus,
        },
        # Region model + the region-aware call for the planted variant, so the
        # annotation tests assert against a prediction rather than whatever the
        # code happens to emit. Hand-verified: genome 201 is CDS frame position
        # 52 -> residue 18, codon TGT -> AGT, C -> S, non-synonymous.
        "cds_start": CDS_START,
        "cds_end": CDS_END,
        "variant_genome_pos": VARIANT_ABS_POS + 1,
        "variant_region": "polyprotein",
        "variant_residue": 18,
        # The 5'UTR variant. Non-coding, so the prediction is the ABSENCE of
        # amino-acid columns -- there is no residue to name.
        "utr_variant_abs_pos": UTR_VARIANT_ABS_POS,
        "utr_variant_genome_pos": UTR_VARIANT_ABS_POS + 1,
        "utr_variant_ref_base": utr_ref_base,
        "utr_variant_alt_base": utr_variant_base,
        "utr_variant_region": "5UTR",
    }
    (DATA / "mapping_expected.json").write_text(json.dumps(expected, indent=2))

    print(f"wrote {DATA}/template.fasta ({TEMPLATE_LEN} nt)")
    print(f"wrote {DATA}/cellA.fa, {DATA}/cellB.fa, {DATA}/cellC.fa "
          f"({READS_PER_CELL} reads each)")
    print(f"wrote {DATA}/regions.gff3 ({len(REGION_ROWS)} regions)")
    print(f"wrote {DATA}/mapping_expected.json (full-length gap-filled predictions)")
    print(f"planted coding variant:     {ref_base}->{variant_base} "
          f"at template pos {VARIANT_ABS_POS} (genome {VARIANT_ABS_POS + 1}, CDS)")
    print(f"planted non-coding variant: {utr_ref_base}->{utr_variant_base} "
          f"at template pos {UTR_VARIANT_ABS_POS} "
          f"(genome {UTR_VARIANT_ABS_POS + 1}, 5'UTR)")
    print("consensus length: {} ({} lead gaps + {} core + {} tail gaps)".format(
        TEMPLATE_LEN, CORE_START, CORE_END - CORE_START, TEMPLATE_LEN - CORE_END))


if __name__ == "__main__":
    build()
