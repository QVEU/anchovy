"""
make_fastq_fixtures.py -- a FASTQ the whole workflow can actually be run on.

Run once (or whenever you intentionally change the scenario):
    python tests/make_fastq_fixtures.py

WHY THIS EXISTS
---------------
workflow/config.yaml is the example the README tells a new user to run first,
and until this fixture it pointed at a tests/data/fastqs that did not exist. It
failed on its first line with "no FASTQs in input_dir" -- so the one command
meant to prove an installation works proved nothing, and a real misconfiguration
produced the identical error.

The two fixtures that already existed could not fill the gap between them:

  * make_fixtures.py       barcodes spliced into RANDOM sequence. Perfect for
                           extraction, but random reads align to nothing, so the
                           pipeline stops at the first mapped stage.
  * make_mapping_fixtures.py  reads that align, but written as PER-CELL FASTAs
                           with no barcodes in them at all -- they enter the
                           workflow at map_cell, past the half this is about.

This fixture is the join: reads that carry a real barcode AND align to a real
template, so `snakemake --configfile workflow/config.yaml` runs the whole thing,
map_reads through report, on data small enough to finish in seconds.

THE SCENARIO
------------
It deliberately reuses tests/data/mapping/template.fasta and its regions.gff3
rather than inventing a third template, so what a mutation MEANS here is already
pinned down and hand-verified: 5'UTR 1-149, CDS 150-500 (+ strand, phase 0),
3'UTR 501-600.

  * 6 cells, one per barcode in the whitelist this script writes NEXT TO the
    reads. That placement is the example's second lesson: a real run usually
    supplies the barcodes called from its own Illumina run by Cell Ranger and
    Seurat, and they belong with the reads they came from, not in the checkout.
    (tests/data/whitelist.txt is not reused: it is frozen against the extract
    golden, and it holds four barcodes differing only in their last three bases
    -- fine for that fixture, needlessly confusable here.)
  * 8 reads each, above the min_reads and cons_min_depth floors in config.yaml.
  * Each read is  [pad] + [signature with barcode + UMI] + [template 50:550).
    Only the last part aligns, so minimap2 emits a LEADING SOFT CLIP -- which is
    not cosmetic: anchovy finds the signature by searching the window upstream of
    offset = readLen - clipReadLen, so a read that aligns end to end has an empty
    search window and no barcode is ever found. Reads shaped like real cDNA
    (adapter and handle hanging off the 5' end) are what make the stage work, and
    this is the fixture that would catch it breaking.
  * Cells 1-3 match the template; cell 4 carries the CDS substitution at genome
    201, cell 5 the 5'UTR substitution at genome 121, and cell 6 carries both.
    So a run produces a coding call, a non-coding call, and a two-mutation
    genotype -- an edge for the networks to draw rather than six isolated nodes.

    THE MAJORITY HAS TO BE CLEAN, which is why three cells carry nothing. Split
    a variant column evenly and the population consensus breaks the tie by
    np.unique ordering, i.e. alphabetically: the reference base comes out
    arbitrary and the cells MATCHING the template get reported as the mutants.
    An earlier draft of this fixture did exactly that. config.yaml also sets
    `reference`, which calls genotypes against the genome instead of the
    population and removes the tie entirely -- belt and braces, because the two
    failures look identical in the output and neither one errors.

Every read is identical within a cell: this fixture exists to prove the wiring,
and a clean consensus makes a wrong one obvious. Within-cell diversity is
covered at the unit level in tests/test_consensus.py.
"""

from __future__ import annotations

import gzip
from pathlib import Path

DATA = Path(__file__).parent / "data"
OUT_DIR = DATA / "fastqs"

# The sample name the workflow derives from this file -- basename minus
# .fastq/.fq/.gz. Outputs land in <input_dir>/results/example/.
FASTQ_NAME = "example.fastq.gz"

# Written beside the reads, and named in config.yaml. See the docstring: this is
# where a run-specific barcode list belongs.
WHITELIST_NAME = "whitelist.txt"

# Six 16-nt barcodes, pairwise Levenshtein distance >= 7, so no read can be
# assigned to the wrong cell however the matcher is tuned and a genuine
# mis-assignment would be unmistakable rather than plausible.
BARCODES = [
    "AAACCCAAGAAACACT",
    "AACGTGATCCTTAGGC",
    "ACTTGATTGCCAGTAC",
    "AGGCTAACGTTCCGAT",
    "CATGCCTAAGGTTACG",
    "CTAGGTCCAATGCAGT",
]

# --- Read layout ------------------------------------------------------------ #
# v3: 16 nt barcode + 12 nt UMI = 28 N. config.yaml sets chemistry: "v3", and the
# two have to agree -- a v2 signature here would leave the UMI 2 nt short and
# every barcode shifted out of its slot.
SIGNATURE = ("CTACACGACGCTCTTCCGATCT" + "N" * 28 + "TTTCTTATAT")
BARCODE_LEN = 16

READS_PER_CELL = 8

# The aligned portion: the same [50:550) window make_mapping_fixtures.py uses, so
# both fixtures cover the same core of the template and the same variants sit
# inside it.
CORE_START, CORE_END = 50, 550

# 0-based template indices, matching make_mapping_fixtures.py. Genome 201 is
# residue 18 of the polyprotein; genome 121 is inside the 5'UTR and has no
# residue at all.
CDS_VARIANT = 200
UTR_VARIANT = 120

NT = "ACGT"


def _read_template() -> str:
    """The mapping fixture's template, as one string."""
    path = DATA / "mapping" / "template.fasta"
    return "".join(line.strip() for line in path.read_text().splitlines()
                   if not line.startswith(">"))


def _write_whitelist() -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / WHITELIST_NAME
    path.write_text("\n".join(BARCODES) + "\n")
    return path


def _pseudo_random(n: int, seed: int) -> str:
    """Deterministic nucleotides, so the fixture is byte-identical every run.

    The same LCG make_mapping_fixtures.py uses, rather than `random`, which would
    make the output depend on the interpreter's seeding details.
    """
    s = seed & 0xFFFFFFFF
    out = []
    for _ in range(n):
        s = (s + 0x6D2B79F5) & 0xFFFFFFFF
        t = (s ^ (s >> 15)) * (1 | s) & 0xFFFFFFFF
        t = (t + ((t ^ (t >> 7)) * (61 | t) & 0xFFFFFFFF)) ^ t & 0xFFFFFFFF
        val = ((t ^ (t >> 14)) & 0xFFFFFFFF) / 4294967296
        out.append(NT[int(val * 4) % 4])
    return "".join(out)


def _substitute(seq: str, index: int) -> str:
    """Change the base at `index` to the next one in ACGT order."""
    base = seq[index]
    return seq[:index] + NT[(NT.index(base) + 1) % 4] + seq[index + 1:]


def _signature_instance(barcode: str, umi: str) -> str:
    n_start = SIGNATURE.index("N")
    n_end = SIGNATURE.rindex("N") + 1
    filled = (barcode + umi)[:n_end - n_start]
    return SIGNATURE[:n_start] + filled + SIGNATURE[n_end:]


# Which planted variants each cell carries, in barcode order. Three clean cells
# keep the template base in the majority at both variant columns; the sixth cell
# carries both variants, which is what gives the networks an edge.
VARIANTS_BY_CELL = [
    [],
    [],
    [],
    [CDS_VARIANT],
    [UTR_VARIANT],
    [CDS_VARIANT, UTR_VARIANT],
]


def build() -> Path:
    template = _read_template()
    barcodes = BARCODES

    if len(barcodes) != len(VARIANTS_BY_CELL):
        raise SystemExit(
            f"{len(barcodes)} barcodes but variants assigned to "
            f"{len(VARIANTS_BY_CELL)} cells. Keep BARCODES and "
            f"VARIANTS_BY_CELL the same length.")
    if sum(1 for v in VARIANTS_BY_CELL if CDS_VARIANT in v) * 2 >= len(barcodes):
        raise SystemExit(
            "half or more of the cells carry the CDS variant. The population "
            "consensus then ties and picks the reference base alphabetically, "
            "inverting every call. Keep the clean cells in the majority.")

    records: list[str] = []
    seed = 7
    for cell_i, (barcode, variants) in enumerate(zip(barcodes, VARIANTS_BY_CELL)):
        cdna = template
        for pos in variants:
            cdna = _substitute(cdna, pos)
        cdna = cdna[CORE_START:CORE_END]

        for read_i in range(READS_PER_CELL):
            seed += 1
            # A distinct UMI per read, because that is what a real library
            # looks like. anchovy records the UMI but does NOT deduplicate on
            # it -- the published method states the sampling is not deep enough
            # per cell to use UMIs for error correction, so the analysis is the
            # cell-level consensus and every read votes.
            umi = _pseudo_random(12, seed)
            # Padding stands in for the rest of the library construct. Its only
            # job is to sit upstream of the aligned part so the signature falls
            # inside the soft clip.
            pad = _pseudo_random(30, seed + 10_000)
            seq = pad + _signature_instance(barcode, umi) + cdna
            name = f"cell{cell_i + 1}_read{read_i + 1}"
            records.append(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / FASTQ_NAME
    # mtime=0 on the GzipFile, not gzip.open: a timestamp in the gzip header
    # would make the fixture a different file on every regeneration and show up
    # as a spurious diff.
    with open(out, "wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=9,
                           mtime=0) as gz:
            gz.write("".join(records).encode())
    return out


if __name__ == "__main__":
    whitelist = _write_whitelist()
    path = build()
    print(f"wrote {whitelist} ({len(BARCODES)} barcodes)")
    print(f"wrote {path} ({len(BARCODES)} cells x {READS_PER_CELL} reads)")
    print("\nRun the workflow on it with:")
    print("    snakemake -s workflow/Snakefile "
          "--configfile workflow/config.yaml --cores 4")
