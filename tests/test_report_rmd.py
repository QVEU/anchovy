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
    chunk = RMD.split("```{r results-cellmap")[1].split("```")[0]
    assert "drop = FALSE" in chunk, (
        "without drop = FALSE the y scale discards cells with no mutations, "
        "and the plot overstates how mutated the population is")
    assert "levels = per_cell$CBC_ID" in chunk
