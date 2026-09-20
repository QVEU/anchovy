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

from anchovy import regions as regionmod


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
    protein = str(Seq(seq).translate())

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
# Region-aware annotation (GFF3-driven) -- the generalized path
# --------------------------------------------------------------------------- #
def _region_annotations_for_variants(variants: list[str], reference: str,
                                     regions: list) -> tuple[pd.DataFrame, dict]:
    """Build the long-format region-annotation table for a set of variant tokens.

    Args:
        variants: mutation tokens like "13A" (1-based genome pos + mutant base).
        reference: full reference sequence.
        regions: list of Region objects (GFF order preserved).

    Returns:
        (region_table, primary) where region_table is the long-format DataFrame
        (one row per mutation x containing region) and primary maps each variant
        token to its PRIMARY-region annotation (first containing region in GFF
        order), used to fill the existing _annot_v3.csv columns.
    """
    long_rows = []
    primary: dict[str, dict] = {}
    for tok in variants:
        mut_base = tok[-1].upper()
        pos = int(tok[:-1])
        wt_base = reference[pos - 1].upper() if 0 <= pos - 1 < len(reference) else "N"
        rows = regionmod.annotate_mutation(regions, pos, wt_base, mut_base, reference)
        for r in rows:
            r = dict(r)
            r["mutants"] = tok           # tie back to the genotype token
            long_rows.append(r)
        # primary = first containing region (GFF order); annotate_mutation
        # returns rows in classify() order, which preserves GFF order.
        primary[tok] = rows[0]
    return pd.DataFrame(long_rows), primary


def _primary_subname(ann: dict) -> tuple[str, str]:
    """Reconstruct the legacy (subName, subClass) from a primary-region row.

    Coding -> "{wtAA}{residue}{mutAA}" / Syn|Non-Syn (matches the old frame-1
    columns exactly when the primary region is a single CDS). Non-coding ->
    a nucleotide-level name and a "non-coding" class, so the column is populated
    sensibly for UTR/intergenic variants (which the old code couldn't annotate).
    """
    if ann["region_type"] == "coding" and ann["wt_aa"] is not None:
        return (f"{ann['wt_aa']}{ann['residue']}{ann['mut_aa']}",
                ann["sub_class"])
    # non-coding / non-translatable: name by nucleotide change, class label
    return (f"{ann['wt_base']}{ann['genome_pos']}{ann['mut_base']}", "non-coding")
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

    # Handle the no-mutations case gracefully: with only reference cells (or no
    # cells), there are no mutation rows to group on. Build a well-formed empty
    # table with the expected columns rather than crashing on groupby.
    if table.empty:
        table = pd.DataFrame(columns=[
            "mutants", "pos", "base", "CBC_ID", "genotype", "BCMutCount"])
    else:
        # BCMutCount per cell; filter out hypermutated cells (< 200), matching R.
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
    if len(table) and table["mutants"].notna().any():
        table["count"] = table.groupby("mutants")["mutants"].transform("size")
    else:
        table["count"] = 0
    table["freq"] = table["count"] / table["total"] if depth else 0
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
        network: bool = True, gff: str | None = None) -> dict:
    """Read genotype table + reference, annotate, optionally build network, write CSVs.

    Two annotation modes:
      - gff is None (legacy): frame-1 amino-acid annotation against the single
        reference, exactly as before. Preserves existing behavior/goldens.
      - gff provided (region-aware): annotation flows through regions.py. Writes
        an additional long-format {out_prefix}_regionAnnotations.csv (one row per
        mutation x containing region), and fills the legacy _annot_v3.csv
        annotation columns from each mutation's PRIMARY region (first in GFF
        order). Network and frequency logic are unchanged in both modes.

    Outputs:
      {out_prefix}_annot_v3.csv
      {out_prefix}_epistaticNetwork.csv, {out_prefix}_genotypeNetwork.csv (if network)
      {out_prefix}_regionAnnotations.csv (if gff)
    """
    reference = Path(reference_file).read_text().strip().upper()
    cons = pd.read_csv(filt_consensus_csv)
    cons["genotype"] = cons["genotype"].fillna("")

    haplocounts = haplo_analysis(cons)

    variants = sorted({t for g in haplocounts["genotype"].unique()
                       for t in str(g).split("_")} - {"reference", ""})

    written = {}

    if gff is not None:
        # --- region-aware path ------------------------------------------------
        regions = regionmod.parse_gff3(gff)
        for w in regionmod.validate_regions(regions):
            print(f"WARNING: {w}")

        region_table, primary = _region_annotations_for_variants(
            variants, reference, regions)

        # Write the full long-format region annotation table.
        region_csv = f"{out_prefix}_regionAnnotations.csv"
        region_table.to_csv(region_csv, index=False)
        written["regions"] = region_csv

        # Fill the legacy annotation columns from the primary region so
        # _annot_v3.csv keeps its schema (now frame-correct).
        anno_rows = []
        for tok in variants:
            p = primary[tok]
            sub_name, sub_class = _primary_subname(p)
            anno_rows.append({
                "pos": p["genome_pos"], "base": p["mut_base"],
                "ref.codon": p["codon_wt"], "ref.AA": p["wt_aa"],
                "mut.codon": p["codon_mut"], "mut.AA": p["mut_aa"],
                "mut.resPos": p["residue"],
                "subName": sub_name, "subClass": sub_class,
            })
        anno = pd.DataFrame(anno_rows)
    else:
        # --- legacy frame-1 path (unchanged) ---------------------------------
        annos = [annotate_mutation(v, reference) for v in variants]
        anno = pd.DataFrame([{
            "pos": a["pos"], "base": a["base"],
            "ref.codon": a["ref"]["codon"] if a["ref"] else None,
            "ref.resPos": a["ref"]["resPos"] if a["ref"] else None,
            "ref.AA": a["ref"]["AA"] if a["ref"] else None,
            "mut.codon": a["mut"]["codon"] if a["mut"] else None,
            "mut.resPos": a["mut"]["resPos"] if a["mut"] else None,
            "mut.AA": a["mut"]["AA"] if a["mut"] else None,
            "subName": a["subName"], "subClass": a["subClass"],
        } for a in annos])

    merged = haplocounts.merge(anno, how="left", on=["pos", "base"])

    merged["genotypeName"] = merged.groupby("genotype")["subName"].transform(
        lambda s: "_".join(pd.unique(s.dropna())))
    merged["genoFreq"] = merged.groupby("genotypeName")["CBC_ID"].transform("nunique") / merged["total"]
    merged["haploFreq"] = merged.groupby("genotype")["CBC_ID"].transform("nunique") / merged["total"]

    merged.to_csv(f"{out_prefix}_annot_v3.csv", index=False)
    written["annot"] = f"{out_prefix}_annot_v3.csv"

    if network:
        single_steps, all_entries = hap_network_gen(haplocounts)
        single_steps.to_csv(f"{out_prefix}_epistaticNetwork.csv", index=False)
        all_entries.to_csv(f"{out_prefix}_genotypeNetwork.csv", index=False)
        written["epistatic"] = f"{out_prefix}_epistaticNetwork.csv"
        written["genotype"] = f"{out_prefix}_genotypeNetwork.csv"

    return {"annot": merged, "written": written}
