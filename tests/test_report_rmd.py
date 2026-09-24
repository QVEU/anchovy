"""
Guards on visualization/anchovy_report.rmd.

The Rmd cannot be executed here -- there is no R in the dev container and CI
runs pytest only -- so these are structural checks on the document's source.
They are cheap and they pin the two things that went wrong by being invisible:
an unbalanced chunk, and a colour pair nobody measured.
"""

from __future__ import annotations

import re
from pathlib import Path

RMD = (Path(__file__).resolve().parent.parent
       / "visualization" / "anchovy_report.rmd").read_text()


def chunk(label: str, code_only: bool = False) -> str:
    """The body of one named chunk.

    Splitting on the next "```" does not work: a chunk may print a fenced block
    of its own -- the network panel's install hint does -- and the naive split
    then stops inside the chunk, so a guard silently checks the first few lines
    instead of the code it names. The closing fence is one at the START of a
    line.

    code_only drops whole-line R comments, which a guard needs whenever it
    forbids a pattern: the comment explaining why the pattern is gone quotes
    the pattern, so a guard reading the comments fails on the very fix it is
    there to protect.
    """
    after = RMD.split("```{r " + label)[1]
    body = re.split(r"^```\s*$", after, maxsplit=1, flags=re.M)[0]
    assert body.strip(), f"chunk {label!r} came out empty"
    if code_only:
        body = "\n".join(line for line in body.splitlines()
                          if not line.lstrip().startswith("#"))
    return body


def test_code_fences_are_balanced():
    assert len(re.findall(r"^```", RMD, flags=re.M)) % 2 == 0


def test_chunk_labels_are_unique():
    # knitr aborts on a duplicate label, and the message names the label rather
    # than the chunk that collided with it.
    labels = re.findall(r"^```\{r ([A-Za-z0-9_-]+)", RMD, flags=re.M)
    assert len(labels) == len(set(labels)), (
        f"duplicate chunk labels: {[l for l in set(labels) if labels.count(l) > 1]}")


# --------------------------------------------------------------------------- #
# Substitution-class colour
# --------------------------------------------------------------------------- #
# The report used #d1495b/#66a182 for Non-Syn/Syn. Measured, that pair puts a
# deuteranope at a CVD separation of 6.5 (OKLab x100, 8 is the target) on the
# one distinction those panels exist to draw; the green also sat below the
# chroma floor, reading grey, and below 3:1 contrast on white. The replacement
# measures 24.7 CVD and 33.6 normal-vision separation.
def test_the_unmeasured_red_green_pair_is_gone():
    for bad in ("#d1495b", "#66a182"):
        assert bad not in RMD, (
            f"{bad} is back: that pair fails the chroma floor and leaves "
            f"deuteranopes unable to separate Syn from Non-Syn")


def test_subclass_colour_is_defined_once_and_shared():
    assert RMD.count("SUBCLASS_PAL <- c(") == 1, (
        "the substitution-class palette should be defined once, so the panels "
        "cannot drift apart")
    # Every panel that splits on substitution class consumes it.
    assert RMD.count("values = SUBCLASS_PAL") >= 3


def test_per_cell_raster_keeps_cells_that_carry_nothing():
    """A blank row is the result for a clean cell, not a cell to drop."""
    body = chunk("results-cellmap")
    assert "drop = FALSE" in body, (
        "without drop = FALSE the y scale discards cells with no mutations, "
        "and the plot overstates how mutated the population is")
    assert "levels = per_cell$CBC_ID" in body


# --------------------------------------------------------------------------- #
# The network panel must keep the reference node
# --------------------------------------------------------------------------- #
# The chunk filtered `target != "reference"` out of the edge list. reference
# only ever appears as a TARGET in the epistatic network, so the node dropped
# out of the drawing while the node table printed immediately above it still
# listed it -- four genotypes in the table, three in the picture.
#
# The worse case is a population derived from a molecular clone, where every
# genotype sits one step off the modal one and nothing else: every edge then
# targets reference, the filter removes all of them, and the panel reports "no
# single-step edges" for a network that is nothing but single-step edges.
def test_the_network_panel_does_not_filter_out_the_reference():
    body = chunk("results-network", code_only=True)
    assert 'target != "reference"' not in body, (
        "the reference node is being filtered out of the network drawing; on a "
        "hub-and-spoke population that empties the panel entirely")
    assert "source != target" in body, (
        "self-edges should still be dropped from the drawing -- they render as "
        "a loop on every node and the CSV keeps them for a different reason")


def test_the_reference_node_is_identifiable_in_the_drawing():
    """Colour separates it, but nothing says which end of the scale it is."""
    body = chunk("results-network")
    assert 'name == "reference"' in body and "geom_node_text" in body, (
        "nothing labels the reference node, so the anchor of the picture is "
        "just another dot")
