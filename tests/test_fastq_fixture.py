"""
Guards on tests/data/fastqs -- the fixture workflow/config.yaml runs on.

That config is the example the README points a new user at, so its failure mode
matters more than most: if the fixture drifts, the command meant to prove an
installation works stops proving it, and a real misconfiguration produces the
same output as a broken fixture.

Running the workflow itself needs minimap2 and R, which CI does not have, so
these check the fixture's structure instead -- which is where the interesting
invariants are anyway. The full run is verified by hand; what it produces is
recorded in tests/make_fastq_fixtures.py.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parent.parent
FASTQ_DIR = REPO / "tests" / "data" / "fastqs"
FASTQ = FASTQ_DIR / "example.fastq.gz"
WHITELIST = FASTQ_DIR / "whitelist.txt"
CONFIG = REPO / "workflow" / "config.yaml"

HANDLE_5 = "CTACACGACGCTCTTCCGATCT"      # 22 nt
HANDLE_3 = "TTTCTTATAT"                  # 10 nt
BARCODE_LEN, UMI_LEN = 16, 12            # v3
CORE_START, CORE_END = 50, 550


def _records():
    lines = gzip.open(FASTQ, "rt").read().splitlines()
    return [(lines[i][1:], lines[i + 1]) for i in range(0, len(lines), 4)]


def _template() -> str:
    path = REPO / "tests" / "data" / "mapping" / "template.fasta"
    return "".join(l.strip() for l in path.read_text().splitlines()
                   if not l.startswith(">"))


# --------------------------------------------------------------------------- #
# config.yaml is a runnable example, not a template to fill in
# --------------------------------------------------------------------------- #
# It used to point at a tests/data/fastqs that did not exist, and to omit the
# window keys, so the command the README opens with failed twice over before
# doing any work.
def test_the_example_config_points_at_files_that_exist():
    cfg = yaml.safe_load(CONFIG.read_text())
    for key in ("input_dir", "template", "reference", "gff", "whitelist"):
        assert key in cfg, f"config.yaml lost its `{key}` key"
        assert (REPO / cfg[key]).exists(), (
            f"config.yaml's {key} -> {cfg[key]} does not exist; the example "
            f"command in the README fails before doing any work")


def test_the_example_config_is_complete_enough_to_run():
    cfg = yaml.safe_load(CONFIG.read_text())
    # The Snakefile raises unless one of these says what to call variants over.
    assert cfg.get("whole_reference") or cfg.get("gff"), (
        "without a gff or whole_reference the consensus stage renumbers "
        "positions and the region lookup has nothing to match")
    assert cfg.get("window") or {"orf_start", "orf_end"} <= cfg.keys(), (
        "the Snakefile refuses to build without an analysis window")
    assert cfg["reference"] == cfg["template"], (
        "`reference` must be the genome the reads were mapped to. Pointing it "
        "elsewhere compares consensuses against the wrong coordinates, which "
        "reports nearly every position as mutant instead of failing")


def test_the_example_writes_inside_the_fixture_directory():
    """Outputs are anchored to input_dir, and that directory is gitignored."""
    cfg = yaml.safe_load(CONFIG.read_text())
    assert cfg["input_dir"] == "tests/data/fastqs"
    ignore = (REPO / ".gitignore").read_text().splitlines()
    assert "results/" in ignore, (
        "results/ is no longer ignored, so running the example dirties the "
        "working tree with tests/data/fastqs/results/")


# --------------------------------------------------------------------------- #
# The reads themselves
# --------------------------------------------------------------------------- #
def test_every_read_is_shaped_for_the_barcode_search():
    """A LEADING SOFT CLIP IS THE WHOLE POINT, and it is easy to lose.

    anchovy finds the signature in the window UPSTREAM of a read's alignment
    offset. A read that aligns end to end has an empty window and no barcode is
    ever found -- so the fixture only exercises extraction while the construct
    (pad + handle + barcode + UMI + handle) sits ahead of the part that aligns.
    """
    template = _template()
    core = template[CORE_START:CORE_END]
    for name, seq in _records():
        head, _, tail = seq.partition(HANDLE_5)
        assert head, f"{name} starts with the handle; nothing would soft-clip"
        assert len(tail) > BARCODE_LEN + UMI_LEN + len(HANDLE_3), name
        aligned = seq[seq.index(HANDLE_3, len(head)) + len(HANDLE_3):]
        assert len(aligned) == len(core), (
            f"{name}'s aligned portion is {len(aligned)} nt, not {len(core)}")
        # At most the two planted substitutions; anything else means the
        # fixture and the template have drifted apart.
        assert sum(a != b for a, b in zip(aligned, core)) <= 2, (
            f"{name} differs from the template in more than the two planted "
            f"positions -- regenerate with tests/make_fastq_fixtures.py")


def test_every_read_carries_an_exact_whitelist_barcode():
    """No correction needed, so a mis-assignment downstream is unambiguous."""
    allowed = set(WHITELIST.read_text().split())
    assert len(allowed) == 6
    for name, seq in _records():
        start = seq.index(HANDLE_5) + len(HANDLE_5)
        assert seq[start:start + BARCODE_LEN] in allowed, name


def test_each_cell_has_enough_reads_to_clear_the_configured_floors():
    from collections import Counter
    cfg = yaml.safe_load(CONFIG.read_text())
    per_cell = Counter(
        seq[seq.index(HANDLE_5) + len(HANDLE_5):][:BARCODE_LEN]
        for _, seq in _records())
    assert len(per_cell) == 6
    floor = max(cfg.get("min_reads", 5), cfg.get("cons_min_depth", 5))
    assert min(per_cell.values()) >= floor, (
        f"a cell has {min(per_cell.values())} reads but config.yaml requires "
        f"{floor}; the example would silently produce fewer cells")


def test_the_umis_are_distinct_within_each_cell():
    """A shared UMI is deduplicated, which would starve the consensus."""
    from collections import defaultdict
    by_cell = defaultdict(list)
    for _, seq in _records():
        start = seq.index(HANDLE_5) + len(HANDLE_5)
        by_cell[seq[start:start + BARCODE_LEN]].append(
            seq[start + BARCODE_LEN:start + BARCODE_LEN + UMI_LEN])
    for barcode, umis in by_cell.items():
        assert len(set(umis)) == len(umis), f"{barcode} reuses a UMI"


# --------------------------------------------------------------------------- #
# The scenario the example is meant to demonstrate
# --------------------------------------------------------------------------- #
def test_the_planted_variants_are_one_coding_and_one_non_coding():
    """The example exists to show both, so a CDS-only fixture is a regression."""
    from anchovy.regions import cds_window
    cds_start, cds_end = cds_window(str(REPO / "tests/data/mapping/regions.gff3"))
    template = _template()

    seen = set()
    for _, seq in _records():
        aligned = seq[seq.index(HANDLE_3, 22) + len(HANDLE_3):]
        for i, (a, b) in enumerate(zip(aligned, template[CORE_START:CORE_END])):
            if a != b:
                seen.add(CORE_START + i)
    assert len(seen) == 2, f"expected two planted positions, found {sorted(seen)}"
    assert any(cds_start <= p < cds_end for p in seen), "no coding variant"
    assert any(not (cds_start <= p < cds_end) for p in seen), "no non-coding variant"


def test_the_clean_cells_stay_in_the_majority():
    """An even split makes the consensus pick the reference base alphabetically.

    The cells matching the template are then reported as the mutants, which is
    what an earlier draft of this fixture did. `reference` in config.yaml also
    guards this; both are cheap and the failure is silent.
    """
    from collections import defaultdict
    template = _template()
    carriers = defaultdict(set)
    for _, seq in _records():
        start = seq.index(HANDLE_5) + len(HANDLE_5)
        barcode = seq[start:start + BARCODE_LEN]
        aligned = seq[seq.index(HANDLE_3, 22) + len(HANDLE_3):]
        for i, (a, b) in enumerate(zip(aligned, template[CORE_START:CORE_END])):
            if a != b:
                carriers[CORE_START + i].add(barcode)
    for position, cells in carriers.items():
        assert len(cells) * 2 < 6, (
            f"{len(cells)} of 6 cells carry the variant at {position}; the "
            f"template base is no longer the clear majority there")


def test_the_fixture_is_reproducible_from_its_generator(tmp_path, monkeypatch):
    """Regenerating must be byte-identical, or it shows up as a spurious diff."""
    import make_fastq_fixtures as gen
    monkeypatch.setattr(gen, "OUT_DIR", tmp_path)
    gen._write_whitelist()
    rebuilt = gen.build()
    assert rebuilt.read_bytes() == FASTQ.read_bytes(), (
        "the committed FASTQ is not what the generator produces now")
    assert (tmp_path / "whitelist.txt").read_text() == WHITELIST.read_text()
