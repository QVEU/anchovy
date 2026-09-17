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
VARIANT_ABS_POS = 200            # absolute template coordinate of planted variant
NT = "ACGT"


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

    # cellB: every read carries the planted variant.
    cb = list(core)
    cb[variant_core_idx] = variant_base
    cellB_reads = ["".join(cb) for _ in range(READS_PER_CELL)]

    def write_cell(name, reads):
        recs = [f">{name}_read{i}\n{seq}" for i, seq in enumerate(reads)]
        (DATA / f"{name}.fa").write_text("\n".join(recs) + "\n")

    write_cell("cellA", cellA_reads)
    write_cell("cellB", cellB_reads)

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
        },
    }
    (DATA / "mapping_expected.json").write_text(json.dumps(expected, indent=2))

    print(f"wrote {DATA}/template.fasta ({TEMPLATE_LEN} nt)")
    print(f"wrote {DATA}/cellA.fa, {DATA}/cellB.fa ({READS_PER_CELL} reads each)")
    print(f"wrote {DATA}/mapping_expected.json (full-length gap-filled predictions)")
    print(f"planted variant: {ref_base}->{variant_base} at template pos {VARIANT_ABS_POS}")
    print("consensus length: {} ({} lead gaps + {} core + {} tail gaps)".format(
        TEMPLATE_LEN, CORE_START, CORE_END - CORE_START, TEMPLATE_LEN - CORE_END))


if __name__ == "__main__":
    build()
