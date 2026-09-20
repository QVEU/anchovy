"""
annotate.py -- the annotation stage: genotype table -> annotated mutations +
genotype/epistatic networks. Python port of Consensus_Annotation.R's ANALYSIS
core (plotting deliberately dropped).

Ported against a golden reference frozen from the original R (see
tests/data/annot/), with every amino-acid call independently hand-verified. The
port reproduces R's behavior faithfully, INCLUDING two quirks noted inline (the
reference self/edge duplication and the network `count` ratio), so its output
matches the R golden. Where the port would otherwise "improve" on the R, it does
not -- behavior changes belong in separate, labeled commits.

STRUCTURE (pure core + thin orchestration, like the other stages):
  codon()            -- codon + amino acid at a nt position (pure)
  annotate_mutation()-- one mutation token -> ref/mut codon, AA, subName, subClass
  haplo_analysis()   -- unroll genotypes -> per-mutation counts/frequencies table
  hap_network_gen()  -- genotype-overlap network edges (pure -> DataFrames)
  run()              -- orchestrate + write the CSVs
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Seq import Seq


# --------------------------------------------------------------------------- #
# Codon / translation
# --------------------------------------------------------------------------- #
def _residue_of(pos: int) -> int:
    """1-based amino-acid residue index containing 1-based nt position `pos`.

    Reproduces R's AApos: (((pos-1)-(pos-1)%%3)/3)+1.
    """
    return ((pos - 1) - ((pos - 1) % 3)) // 3 + 1


def codon(seq: str, pos: int) -> dict:
    """Codon sequence, residue position, and amino acid at 1-based nt `pos`.

    Mirrors the R `codon()`: translate the whole sequence (standard genetic
    code, no init codon), pick the residue containing `pos`, and extract the
    3-nt codon by the same frame math. Returns dict(codon, resPos, AA).
    """
    if pos is None or (isinstance(pos, float) and np.isnan(pos)):
        return {"codon": "WT", "resPos": 0, "AA": "WT"}

    # Translate full sequence with the standard genetic code. Biopython's
    # translate matches R's Biostrings translate for the standard table.
    #
    # Gap characters are mapped to 'N' first. Biopython RAISES on a codon like
    # '--A', so without this a reference carrying gaps crashes the whole stage --
    # which is reachable now that the consensus stage can keep the whole
    # reference, whose uncovered flanks sam2consensus gap-fills. 'N' translates
    # to 'X' ("not translatable"), which is the honest answer for a position with
    # no coverage. On a gap-free reference this substitution does nothing, so
    # every frozen golden is untouched.
    protein = str(Seq(seq.replace("-", "N")).translate())

    res = _residue_of(pos)
    os_ = pos % 3
    # Codon nt bounds (1-based inclusive in R); convert to 0-based slice.
    if os_ == 0:
        first = (pos - 3 + 1) - 1
        last = pos            # slice end (exclusive) == R's inclusive `last`
    else:
        first = (pos - os_ + 1) - 1
        last = (pos - os_ + 3)
    codon_seq = seq[first:last]

    # AA at the residue (1-based -> 0-based index). Guard against out-of-range.
    aa = protein[res - 1] if 0 <= res - 1 < len(protein) else "X"
    return {"codon": codon_seq, "resPos": res, "AA": aa}


def annotate_mutation(variant: str, seq: str) -> dict:
    """Annotate one mutation token like '13A' against reference `seq`.

    Returns dict with pos, base, ref (codon dict), mut (codon dict), subName
    (e.g. 'R5S'), subClass ('Syn'/'Non-Syn'/'X'). The 'reference' sentinel maps
    to the WT row, matching the R.
    """
    if variant == "reference":
        return {"pos": 0, "base": None, "ref": None, "mut": None,
                "subName": "ref", "subClass": "WT"}

    nt = variant[-1].upper()
    pos = int(variant[:-1])

    ref = codon(seq, pos)
    # Build the single-substitution mutant sequence (1-based pos).
    mut_seq = seq[:pos - 1] + nt + seq[pos:]
    mut = codon(mut_seq, pos)

    if mut["AA"] is not None and mut["AA"] != "X":
        sub_class = "Syn" if ref["AA"] == mut["AA"] else "Non-Syn"
    else:
        sub_class = "X"

    sub_name = f"{ref['AA']}{mut['resPos']}{mut['AA']}"
    return {"pos": pos, "base": nt, "ref": ref, "mut": mut,
            "subName": sub_name, "subClass": sub_class}


# --------------------------------------------------------------------------- #
# Haplotype analysis
# --------------------------------------------------------------------------- #
def haplo_analysis(cons: pd.DataFrame) -> pd.DataFrame:
    """Unroll genotype strings into per-mutation records with counts/frequencies.

    Port of R haploanalysis(). Input needs CBC_ID and genotype columns.
    Returns a long table: one row per (cell, mutation), plus explicit reference
    rows for cells with empty genotype. Columns include mutants, pos, base,
    CBC_ID, genotype, count, freq, total, BCMutCount.
    """
    depth = cons["CBC_ID"].nunique()

    rows = []
    reference_cells = []
    for _, r in cons.iterrows():
        geno = r["genotype"]
        if geno == "" or pd.isna(geno):
            reference_cells.append(r["CBC_ID"])
            continue
        for tok in str(geno).split("_"):
            base = tok[-1].upper()
            if base == "-":
                continue  # R drops deletion tokens (base == "-")
            pos = int(tok[:-1])
            rows.append({"mutants": tok, "pos": pos, "base": base,
                         "CBC_ID": r["CBC_ID"], "genotype": geno})

    table = pd.DataFrame(rows)

    # BCMutCount per cell; filter out hypermutated cells (< 200), matching R.
    if not table.empty:
        table["BCMutCount"] = table.groupby("CBC_ID")["mutants"].transform("size")
        table = table[table["BCMutCount"] < 200]

    # Explicit reference rows so ref cells appear downstream.
    ref_rows = pd.DataFrame([{
        "mutants": "", "pos": np.nan, "base": np.nan,
        "CBC_ID": c, "genotype": "reference", "BCMutCount": 0,
    } for c in reference_cells])

    table = pd.concat([table, ref_rows], ignore_index=True)

    # Frequencies across cells: count is number of rows per mutant token.
    table["total"] = depth
    table["count"] = table.groupby("mutants")["mutants"].transform("size")
    table["freq"] = table["count"] / table["total"]
    return table


# --------------------------------------------------------------------------- #
# Network generation
# --------------------------------------------------------------------------- #
def hap_network_gen(haplocounts: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Genotype-overlap network. Port of R hapNetworkGen().

    Returns (single_steps, all_entries) as DataFrames matching the columns of
    the R _epistaticNetwork.csv and _genotypeNetwork.csv respectively.

    QUIRKS reproduced faithfully from the R (not bugs to fix here):
      - the `count` column is rowSums(countMatrix)/rowSums(binaryMatrix) with the
        count-matrix aggregation defaulting to occurrence counts (the R dcast
        'fun.aggregate defaulting to length()' warning), so it is a ratio, not a
        simple cell count.
      - reference genotypes get a duplicated self/edge row with mutNumTarget=0.
    """
    # binary presence matrix: genotype x mutants, 1 if freq>0 for that pair.
    hc = haplocounts.copy()
    binary = (hc.assign(present=(hc["freq"] > 0).astype(int))
                .pivot_table(index="genotype", columns="mutants",
                             values="present", aggfunc=lambda x: int(any(v > 0 for v in x)),
                             fill_value=0))
    # count matrix: occurrences per (genotype, mutants) -- matches R's dcast
    # defaulting to length() (the source of the warning), giving per-pair counts.
    count_mat = hc.pivot_table(index="genotype", columns="mutants",
                               values="count", aggfunc="size", fill_value=0)

    # align columns
    count_mat = count_mat.reindex(columns=binary.columns, fill_value=0)

    genotypes = list(binary.index)

    # per-genotype count metadata: rowSums(count)/rowSums(binary)
    row_bin = binary.sum(axis=1).replace(0, np.nan)
    row_cnt = count_mat.sum(axis=1)
    counts = pd.DataFrame({
        "genotype": binary.index,
        "count": (row_cnt / row_bin).values,
    })
    # Render whole-number ratios as ints (matches R's integer formatting, e.g.
    # "2" not "2.0"), so the network CSV is byte-identical to the R golden.
    counts["count"] = counts["count"].apply(
        lambda x: int(x) if pd.notna(x) and float(x).is_integer() else x)

    # pairwise overlap, R's upper-triangular traversal (removes i from future j's)
    entries = []
    geno_list = list(genotypes)
    for i in genotypes:
        for j in geno_list:
            src = binary.loc[i]
            tgt = binary.loc[j]
            overlap = int(((src == 1) & (tgt == 1)).sum())
            if overlap > 0:
                entries.append({
                    "source": i, "target": j, "overlap": overlap,
                    "mutNumSource": len(str(i).split("_")),
                    "mutNumTarget": len(str(j).split("_")),
                })
        geno_list = [g for g in geno_list if g != i]

    all_entries = pd.DataFrame(entries)
    # attach count metadata (merge on source == genotype)
    all_entries = counts.merge(all_entries, left_on="genotype",
                               right_on="source").drop(columns=["source"])

    # single-step edges: differ by exactly one mutation
    ss = all_entries[
        ((all_entries["mutNumSource"] == all_entries["overlap"]) &
         (all_entries["mutNumTarget"] == all_entries["overlap"] + 1)) |
        ((all_entries["mutNumTarget"] == all_entries["overlap"]) &
         (all_entries["mutNumSource"] == all_entries["overlap"] + 1))
    ]
    self_steps = all_entries[all_entries["genotype"] == all_entries["target"]]

    # reference edges: single-mutation self-steps re-pointed at "reference"
    ref_edges = self_steps[
        (self_steps["genotype"] == self_steps["target"]) &
        (self_steps["mutNumTarget"] == 1)
    ].copy()
    ref_edges["target"] = "reference"
    ref_edges["mutNumTarget"] = 0

    single_steps = pd.concat([ss, self_steps, ref_edges], ignore_index=True)

    # Rename the source-node column from "genotype" to "source" so Cytoscape
    # auto-detects the source/target roles on import (no manual column mapping).
    # This intentionally diverges from the R output's "genotype" header -- a
    # deliberate downstream-usability change, scoped to the network CSVs only
    # (the annotation table keeps "genotype", where that name is correct).
    single_steps = single_steps.rename(columns={"genotype": "source"})
    all_entries = all_entries.rename(columns={"genotype": "source"})
    return single_steps, all_entries


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def run(filt_consensus_csv: str, reference_file: str, out_prefix: str,
        network: bool = True) -> dict:
    """Read genotype table + reference, annotate, optionally build network, write CSVs.

    Outputs (matching the R naming):
      {out_prefix}_annot_v3.csv
      {out_prefix}_epistaticNetwork.csv, {out_prefix}_genotypeNetwork.csv (if network)
    """
    reference = Path(reference_file).read_text().strip().upper()
    cons = pd.read_csv(filt_consensus_csv)
    cons["genotype"] = cons["genotype"].fillna("")

    haplocounts = haplo_analysis(cons)

    # annotate each distinct variant token (excluding the reference sentinel)
    variants = sorted({t for g in haplocounts["genotype"].unique()
                       for t in str(g).split("_")} - {"reference", ""})
    annos = [annotate_mutation(v, reference) for v in variants]
    anno_rows = [{
        "pos": a["pos"], "base": a["base"],
        "ref.codon": a["ref"]["codon"] if a["ref"] else None,
        "ref.resPos": a["ref"]["resPos"] if a["ref"] else None,
        "ref.AA": a["ref"]["AA"] if a["ref"] else None,
        "mut.codon": a["mut"]["codon"] if a["mut"] else None,
        "mut.resPos": a["mut"]["resPos"] if a["mut"] else None,
        "mut.AA": a["mut"]["AA"] if a["mut"] else None,
        "subName": a["subName"], "subClass": a["subClass"],
    } for a in annos]
    anno = pd.DataFrame(anno_rows)

    merged = haplocounts.merge(anno, how="left", on=["pos", "base"])

    # genotypeName = "_".join(unique subName) per genotype; geno/haplo freqs
    merged["genotypeName"] = merged.groupby("genotype")["subName"].transform(
        lambda s: "_".join(pd.unique(s.dropna())))
    merged["genoFreq"] = merged.groupby("genotypeName")["CBC_ID"].transform("nunique") / merged["total"]
    merged["haploFreq"] = merged.groupby("genotype")["CBC_ID"].transform("nunique") / merged["total"]

    written = {}
    merged.to_csv(f"{out_prefix}_annot_v3.csv", index=False)
    written["annot"] = f"{out_prefix}_annot_v3.csv"

    if network:
        single_steps, all_entries = hap_network_gen(haplocounts)
        single_steps.to_csv(f"{out_prefix}_epistaticNetwork.csv", index=False)
        all_entries.to_csv(f"{out_prefix}_genotypeNetwork.csv", index=False)
        written["epistatic"] = f"{out_prefix}_epistaticNetwork.csv"
        written["genotype"] = f"{out_prefix}_genotypeNetwork.csv"

    return {"annot": merged, "written": written}
