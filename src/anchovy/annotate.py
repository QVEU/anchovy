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
  annotate_variants_by_region() -- GFF3 region-aware annotation (long format)
  run()              -- orchestrate + write the CSVs

REGION AWARENESS (added later, opt-in via run(gff=...)) layers a GFF3-driven
model over the above without disturbing it. See the run() docstring and
regions.py. With no GFF supplied the legacy single-ORF, frame-1 path runs
exactly as before, which is why the R-frozen goldens still pass.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Seq import Seq

from anchovy import regions as regions_mod


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
# Region-aware annotation (GFF3-driven)
# --------------------------------------------------------------------------- #
# The legacy path above assumes ONE reference read in frame 1 from its first
# base. That is wrong in two ways for a real viral genome: it cannot annotate
# non-coding regions at all, and it mis-numbers residues whenever the ORF does
# not start at position 1. The functions below replace that assumption with the
# GFF3 region model in regions.py, while leaving the legacy path byte-identical
# for callers that pass no GFF (gff=None).

REGION_COLUMNS = [
    "mutation_id", "region", "region_type", "strand", "genome_pos",
    "wt_base", "mut_base", "residue", "codon_wt", "codon_mut",
    "wt_aa", "mut_aa", "sub_class",
]


def _token_parts(token: str) -> tuple[int, str]:
    """Split a genotype token like '201A' into (position, base)."""
    return int(token[:-1]), token[-1].upper()


def annotate_variants_by_region(variants: list[str], reference: str,
                                gff: str) -> pd.DataFrame:
    """Annotate genotype tokens against every containing GFF3 region (long format).

    One row per (mutation x containing region): a variant inside both a
    polyprotein and a mature protein yields two rows, each numbered in its own
    frame of reference. Positions in no region yield a single 'intergenic' row.

    `variants` are genotype tokens whose positions are GENOME coordinates, which
    is only true when the consensus stage ran in whole-reference mode
    (consensus.run(trim=False)). Feeding it legacy trimmed tokens would silently
    annotate the wrong positions.
    """
    region_list = regions_mod.parse_gff3(gff)
    for message in regions_mod.validate_regions(region_list):
        warnings.warn(f"{gff}: {message}", stacklevel=2)

    rows: list[dict] = []
    for token in variants:
        pos, mut_base = _token_parts(token)
        # The WT base is the reference base at that genome position. This is
        # exact, not inferred: a token exists only where the cell differs from
        # the reference the genotypes were called against.
        wt_base = reference[pos - 1] if 1 <= pos <= len(reference) else "N"
        rows.extend(regions_mod.annotate_mutation(
            region_list, pos, wt_base, mut_base, reference))

    table = pd.DataFrame(rows, columns=REGION_COLUMNS)
    if not table.empty:
        # STABLE sort, and deliberately NOT by region name. regions.annotate_mutation
        # emits a mutation's rows in GFF3 file order, and "the primary region is the
        # first one listed in the GFF" is the contract _legacy_rows_from_regions
        # relies on to back-fill the legacy columns. A non-stable sort, or one that
        # included `region` as a key, would silently reorder those rows and make the
        # primary region whichever one sorts first alphabetically instead.
        table = table.sort_values(
            ["genome_pos", "mut_base"], kind="stable").reset_index(drop=True)
    return table


def _legacy_rows_from_regions(region_table: pd.DataFrame) -> pd.DataFrame:
    """Back-fill the legacy _annot_v3 columns from each mutation's PRIMARY region.

    The primary region is the first one in GFF3 order that contains the position
    (regions.annotate_mutation preserves file order, so it is the first row for
    that mutation). Keeping the legacy columns populated means the annotated
    table, the genotype names, and every downstream consumer keep working
    unchanged when a GFF is supplied -- region awareness ADDS a file, it does not
    reshape the existing one.

    Verified: for a single + strand CDS starting at genome position 1 with phase
    0, this reproduces the legacy frame-1 calls exactly (see
    test_annotate_region_matches_legacy_on_single_cds).
    """
    if region_table.empty:
        return pd.DataFrame(columns=[
            "pos", "base", "ref.codon", "ref.resPos", "ref.AA",
            "mut.codon", "mut.resPos", "mut.AA", "subName", "subClass"])

    primary = region_table.drop_duplicates(subset=["genome_pos", "mut_base"],
                                           keep="first")
    rows = []
    for r in primary.itertuples(index=False):
        # pd.notna, NOT `is not None`: building the frame coerces the None cells of
        # a mixed column to NaN, and `NaN is not None` is True. Without this, a
        # coding row that has no complete codon (a position inside the feature's
        # phase lead-in) reaches int(NaN) and raises.
        if (r.region_type == "coding"
                and pd.notna(r.wt_aa) and pd.notna(r.residue)):
            sub_name = f"{r.wt_aa}{int(r.residue)}{r.mut_aa}"
        else:
            # Non-coding (and intergenic) mutations have no amino-acid name. Use
            # the genome-anchored mutation id so genotype names stay meaningful
            # and distinct instead of collapsing to a shared blank.
            sub_name = r.mutation_id
        rows.append({
            "pos": r.genome_pos,
            "base": r.mut_base,
            "ref.codon": r.codon_wt,
            "ref.resPos": r.residue,
            "ref.AA": r.wt_aa,
            "mut.codon": r.codon_mut,
            "mut.resPos": r.residue,
            "mut.AA": r.mut_aa,
            "subName": sub_name,
            "subClass": r.sub_class if pd.notna(r.sub_class) else "non-coding",
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def run(filt_consensus_csv: str, reference_file: str, out_prefix: str,
        network: bool = True, gff: str | None = None) -> dict:
    """Read genotype table + reference, annotate, optionally build network, write CSVs.

    Outputs (matching the R naming):
      {out_prefix}_annot_v3.csv
      {out_prefix}_epistaticNetwork.csv, {out_prefix}_genotypeNetwork.csv (if network)
      {out_prefix}_regionAnnotations.csv                                  (if gff)

    TWO MODES, and the default is the old one:

    gff=None -- LEGACY. Annotate every mutation against the single reference read
        in frame 1 from its first base. Unchanged, still golden-tested.

    gff=<path> -- REGION-AWARE. Mutation positions are treated as GENOME
        coordinates and annotated against every GFF3 region containing them, in
        long format, written to the extra _regionAnnotations.csv. The legacy
        _annot_v3.csv is still produced, with its columns back-filled from each
        mutation's primary region, so nothing downstream has to change.

        This requires the consensus stage to have run in whole-reference mode
        (consensus.run(trim=False)); legacy trimmed tokens are numbered from the
        ORF start and would annotate the wrong genome positions.
    """
    reference = Path(reference_file).read_text().strip().upper()
    cons = pd.read_csv(filt_consensus_csv)
    cons["genotype"] = cons["genotype"].fillna("")

    haplocounts = haplo_analysis(cons)

    # annotate each distinct variant token (excluding the reference sentinel)
    variants = sorted({t for g in haplocounts["genotype"].unique()
                       for t in str(g).split("_")} - {"reference", ""})

    written = {}
    region_table = None

    if gff is None:
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
    else:
        region_table = annotate_variants_by_region(variants, reference, gff)
        anno = _legacy_rows_from_regions(region_table)

    merged = haplocounts.merge(anno, how="left", on=["pos", "base"])

    # genotypeName = "_".join(unique subName) per genotype; geno/haplo freqs
    merged["genotypeName"] = merged.groupby("genotype")["subName"].transform(
        lambda s: "_".join(pd.unique(s.dropna())))
    merged["genoFreq"] = merged.groupby("genotypeName")["CBC_ID"].transform("nunique") / merged["total"]
    merged["haploFreq"] = merged.groupby("genotype")["CBC_ID"].transform("nunique") / merged["total"]

    merged.to_csv(f"{out_prefix}_annot_v3.csv", index=False)
    written["annot"] = f"{out_prefix}_annot_v3.csv"

    if region_table is not None:
        region_table.to_csv(f"{out_prefix}_regionAnnotations.csv", index=False)
        written["regions"] = f"{out_prefix}_regionAnnotations.csv"

    if network:
        single_steps, all_entries = hap_network_gen(haplocounts)
        single_steps.to_csv(f"{out_prefix}_epistaticNetwork.csv", index=False)
        all_entries.to_csv(f"{out_prefix}_genotypeNetwork.csv", index=False)
        written["epistatic"] = f"{out_prefix}_epistaticNetwork.csv"
        written["genotype"] = f"{out_prefix}_genotypeNetwork.csv"

    return {"annot": merged, "regions": region_table, "written": written}
