"""
test_ambiguous_bases.py -- IUPAC ambiguity codes are uncertainty, not variants.

sam2consensus does not only emit A/C/G/T. Where a position's reads disagree it
emits an IUPAC ambiguity code -- R for A-or-G, Y for C-or-T -- and a lowercase
form when a gap or an N was among the observed bases. Those say the reads
DISAGREED, not that the cell carries a mutation.

Treated as variants they are worse than useless, because the genotype machinery
works on token identity: two cells both reading R at a position are grouped as
sharing a mutation when what they share is uncertainty. On a real run 6 of 11
distinct tokens were ambiguity codes.

They also HIDE structure. A genotype padded with uncalled positions no longer
sits in a subset relationship with its neighbours, and single-step edges are
defined by exactly those relationships -- see the regression test at the bottom,
taken from real output.
"""

from __future__ import annotations

import pandas as pd
import pytest

from anchovy.annotate import hap_network_gen, haplo_analysis
from anchovy.config import ConsensusConfig
from anchovy.consensus import CALLED_BASES, genotype_summary, run


# --------------------------------------------------------------------------- #
# What counts as a called base
# --------------------------------------------------------------------------- #
def test_only_acgt_is_a_called_base():
    assert CALLED_BASES == frozenset("ACGT")


@pytest.mark.parametrize("code", ["R", "Y", "S", "W", "K", "M", "B", "D", "H",
                                  "V", "N"])
def test_uppercase_iupac_codes_are_not_called(code):
    """Every ambiguity code sam2consensus can emit."""
    seqs = ["AA", "A" + code]
    genos, _ = genotype_summary(seqs, reference="AA", skip_ambiguous=True)
    assert genos == ["", ""], f"{code} was treated as a variant"


@pytest.mark.parametrize("code", ["a", "c", "g", "t", "r", "y", "s", "n"])
def test_lowercase_codes_are_not_called(code):
    """Lowercase marks a gap or N among the observed bases -- still uncertain.

    These matter doubly: compared case-sensitively against an uppercase
    reference, a lowercase base looks like a difference at EVERY such position.
    """
    seqs = ["AA", "A" + code]
    genos, _ = genotype_summary(seqs, reference="AA", skip_ambiguous=True)
    assert genos == ["", ""], f"{code!r} was treated as a variant"


def test_real_bases_still_call():
    """The filter must not swallow actual variants."""
    genos, _ = genotype_summary(["AA", "AC"], reference="AA",
                                skip_ambiguous=True)
    assert genos == ["", "2C"]


def test_an_ambiguous_reference_base_also_blocks_the_call():
    """Both sides must be definite: comparing to an unresolved reference base
    says nothing about the cell."""
    genos, _ = genotype_summary(["AC"], reference="AR", skip_ambiguous=True)
    assert genos == [""]


# --------------------------------------------------------------------------- #
# Opting out, and reporting
# --------------------------------------------------------------------------- #
def test_keeping_them_is_possible():
    """Genuine within-cell mixed populations are a real thing to look for."""
    genos, _ = genotype_summary(["AA", "AR"], reference="AA",
                                skip_ambiguous=False)
    assert genos == ["", "2R"]


def test_the_count_is_reported_not_swallowed():
    """Dropped data must be visible, or it costs someone an afternoon."""
    stats: dict = {}
    genotype_summary(["AAA", "ARY"], reference="AAA", skip_ambiguous=True,
                     stats=stats)
    assert stats["ambiguous_calls_skipped"] == 2


def test_run_filters_by_default_in_whole_reference_mode(tmp_path):
    merged = tmp_path / "x_allConsensus.fasta"
    merged.write_text(">c1 ref coverage:9 length:4\nACGT\n"
                      ">c2 ref coverage:9 length:4\nARGT\n")

    filtered = run(fasta=str(merged), trim=False, reference="ACGT",
                   config=ConsensusConfig(depth_min=1),
                   out_prefix=str(tmp_path / "f"))
    assert [r["genotype"] for r in filtered["records"]] == ["", ""]
    assert filtered["stats"]["ambiguous_calls_skipped"] == 1

    kept = run(fasta=str(merged), trim=False, reference="ACGT",
               keep_ambiguous=True, config=ConsensusConfig(depth_min=1),
               out_prefix=str(tmp_path / "k"))
    assert [r["genotype"] for r in kept["records"]] == ["", "2R"]


def test_the_legacy_trimmed_path_is_untouched(tmp_path):
    """Filtering is whole-reference-mode behavior; the frozen path keeps its own."""
    merged = tmp_path / "x_allConsensus.fasta"
    merged.write_text(">c1 ref coverage:9 length:4\nACGT\n"
                      ">c2 ref coverage:9 length:4\nARGT\n")
    legacy = run(fasta=str(merged), start=0, end=4,
                 config=ConsensusConfig(depth_min=1),
                 out_prefix=str(tmp_path / "l"))
    assert any(r["genotype"] for r in legacy["records"])


# --------------------------------------------------------------------------- #
# The regression that motivated this, from real SRR28178313 output
# --------------------------------------------------------------------------- #
def test_filtering_reveals_single_step_edges_hidden_by_ambiguity():
    """Ambiguity codes do not merely add noise -- they DESTROY real structure.

    A single-step edge needs one genotype to be a subset of another, differing
    by exactly one mutation. Pad the larger genotype with uncalled positions and
    that relationship disappears.

    These five genotypes are verbatim from a real run. With the ambiguity codes
    in, the epistatic network has no genotype-to-genotype edges at all. With
    them removed, two appear -- and they are the actual evolutionary
    relationships in the data.
    """
    def edges(genotypes):
        cons = pd.DataFrame([{"CBC_ID": f"c{i}", "genotype": g}
                             for i, g in enumerate(genotypes)])
        ss, _ = hap_network_gen(haplo_analysis(cons))
        return ss[(ss["source"] != ss["target"])
                  & (ss["target"] != "reference")]

    with_codes = ["1294C_2873C_3176S_3185R_3587R_3751C_6058R", "1294C_3751C",
                  "2873C", "3288T_5059C_5431Y", "3288T_6202R"]
    without = ["1294C_2873C_3751C", "1294C_3751C", "2873C",
               "3288T_5059C", "3288T"]

    assert len(edges(with_codes)) == 0
    found = {(e["source"], e["target"]) for _, e in edges(without).iterrows()}
    assert ("1294C_2873C_3751C", "1294C_3751C") in found
    assert ("3288T", "3288T_5059C") in found


def test_workflow_forwards_the_flag():
    from pathlib import Path
    snakefile = (Path(__file__).resolve().parent.parent
                 / "workflow" / "Snakefile").read_text()
    assert "--keep-ambiguous" in snakefile
    assert "KEEP_AMBIGUOUS" in snakefile
