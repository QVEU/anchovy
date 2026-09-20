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
  genotype_nodes()   -- per-genotype NODE attributes for those edges (pure)
  annotate_variants_by_region() -- GFF3 region-aware annotation (long format)
  run()              -- orchestrate + write the CSVs

REGION AWARENESS (added later, opt-in via run(gff=...)) layers a GFF3-driven
model over the above without disturbing it. See the run() docstring and
regions.py. With no GFF supplied the legacy single-ORF, frame-1 path runs
exactly as before, which is why the R-frozen goldens still pass.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from Bio.Seq import Seq

from anchovy import regions as regions_mod
from anchovy.io import read_reference_sequence


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

    # Name the columns explicitly. Built from an empty `rows` list, a bare
    # DataFrame() has no columns at all, and the groupby below then fails with
    # KeyError: 'mutants' -- which is what a run where every cell was too sparse
    # to yield a consensus used to produce, several stages after the real
    # problem. An empty table with the right shape flows through instead.
    table = pd.DataFrame(rows, columns=["mutants", "pos", "base",
                                        "CBC_ID", "genotype"])

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
def hap_network_gen(haplocounts: pd.DataFrame,
                    self_edges: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Genotype-overlap network. Port of R hapNetworkGen().

    Returns (single_steps, all_entries) as DataFrames matching the columns of
    the R _epistaticNetwork.csv and _genotypeNetwork.csv respectively.

    Args:
        haplocounts: the long table from haplo_analysis().
        self_edges: whether to keep edges from a genotype to ITSELF in the
            single-step (epistatic) network. Default False, since Cytoscape
            draws each one as a loop on the node and they carry no information
            there. Pass True for output byte-comparable with the R.

    SELF-EDGES, and why only the epistatic network drops them
    ---------------------------------------------------------
    The pairwise traversal compares every genotype with itself, so each one
    gets an i == j row. In the epistatic network those are pure decoration:
    every genotype is also reachable by a step edge or by the reference edge
    below, so dropping them loses no node. Measured on the annotation fixture:
    6 of 10 rows are self-edges and removing them loses nothing.

    In all_entries they are LOAD-BEARING and are therefore always kept. A
    genotype sharing no mutation with any other appears there ONLY as its own
    self-edge, so dropping them deletes it from the graph outright -- on the
    same fixture that is 3 of 5 genotypes, including the reference. If you want
    a self-edge-free genotype network, take the node list from
    {prefix}_genotypeNodes.csv, which is complete either way, and filter the
    edges yourself knowing what it costs.

    THE `count` COLUMN is the number of cells carrying that genotype, despite
    being computed the R's roundabout way as
    rowSums(countMatrix)/rowSums(binaryMatrix) -- the count matrix aggregation
    defaults to occurrence counts (the R dcast 'fun.aggregate defaulting to
    length()' warning). It reduces exactly: every cell with genotype g
    contributes one row per mutation in g, so for n cells and k distinct
    mutations it is (n*k)/k = n. An earlier version of this docstring called it
    "a ratio, not a simple cell count", which was wrong. It is pinned against
    the node table's nCells in tests/test_network_nodes.py.

    QUIRK still reproduced faithfully from the R (not a bug to fix here):
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
    if not self_edges:
        # Filter the RESULT, not its inputs. Self-loops arrive by two routes:
        # self_steps, and the reference's own ref_edges row (the R quirk noted
        # above re-points a 1-mutation self-step at "reference", which for the
        # reference genotype is itself). Omitting self_steps from the concat
        # would leave that second one behind.
        single_steps = single_steps[
            single_steps["genotype"] != single_steps["target"]
        ].reset_index(drop=True)

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
# Genotype node table (for network visualization)
# --------------------------------------------------------------------------- #
NODE_COLUMNS = ["genotype", "genotypeName", "nMutations", "nCells",
                "genoFreq", "haploFreq"]


def genotype_nodes(annotated: pd.DataFrame) -> pd.DataFrame:
    """One row per genotype: the NODE attributes for the network CSVs (pure).

    WHY THIS EXISTS. Both network files are entirely EDGE-level -- source,
    target, overlap, mutNumSource, mutNumTarget. Cytoscape therefore draws nodes
    with no attributes at all, so a genotype cannot be sized by how common it is
    or labelled by what it does to the protein without hand-joining CSVs inside
    Cytoscape. Every one of those values already exists in the annotated table;
    it just never reached a file keyed by genotype.

    `genotype` matches the source/target values in the network CSVs exactly,
    including the "reference" node, so Cytoscape can key a node-table import on
    it directly (File -> Import -> Table from File).

    genotypeName carries the amino-acid-level name ("R5S", "D3V_R5S"), which is
    what makes a rendered network readable -- and which is frame-correct and
    region-aware when annotate ran with a GFF, so a node can read "5UTR:A121C"
    rather than a bare nucleotide token.
    """
    if annotated.empty:
        return pd.DataFrame(columns=NODE_COLUMNS)

    rows = []
    for genotype, group in annotated.groupby("genotype", dropna=False):
        name = group["genotypeName"].dropna()
        # Reference cells have no substitutions, so no amino-acid name; label the
        # node by what it is rather than leaving the field blank.
        label = name.iloc[0] if len(name) and name.iloc[0] else None
        if genotype == "reference" or not label:
            label = "reference" if genotype == "reference" else str(genotype)

        # Mutation count from the genotype string itself, so it agrees with the
        # network's mutNumSource/mutNumTarget rather than being recomputed from
        # a different source.
        n_mutations = 0 if genotype == "reference" else len(str(genotype).split("_"))

        rows.append({
            "genotype": genotype,
            "genotypeName": label,
            "nMutations": n_mutations,
            "nCells": group["CBC_ID"].nunique(),
            "genoFreq": group["genoFreq"].iloc[0] if "genoFreq" in group else None,
            "haploFreq": group["haploFreq"].iloc[0] if "haploFreq" in group else None,
        })

    return (pd.DataFrame(rows, columns=NODE_COLUMNS)
            .sort_values(["nMutations", "genotype"]).reset_index(drop=True))


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def run(filt_consensus_csv: str, reference_file: str, out_prefix: str,
        network: bool = True, gff: str | None = None,
        self_edges: bool = False) -> dict:
    """Read genotype table + reference, annotate, optionally build network, write CSVs.

    Outputs (matching the R naming):
      {out_prefix}_annot_v3.csv
      {out_prefix}_epistaticNetwork.csv, {out_prefix}_genotypeNetwork.csv (if network)
      {out_prefix}_genotypeNodes.csv                                      (if network)
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
    # Shared reader: accepts FASTA or a raw sequence file. Reading this by hand
    # would splice a '>' header into the sequence and shift every genome
    # coordinate, silently, which is exactly the kind of wrong this stage must
    # not be (see io.read_reference_sequence).
    reference = read_reference_sequence(reference_file).upper()
    cons = pd.read_csv(filt_consensus_csv)

    # No cells at all is a failed run, not an empty one, and writing a set of
    # empty CSVs would hide that. It happens when every barcode was too sparse
    # to reach the consensus stage's depth thresholds -- common on a small slice
    # of a real run, where the reads spread thinly over thousands of barcodes.
    if cons.empty:
        raise ValueError(
            f"{filt_consensus_csv} contains no cells, so there is nothing to "
            f"annotate.\n"
            f"  Every barcode was filtered out before this point. The usual "
            f"causes, in order:\n"
            f"    - too few reads overall (are you running on a subsample?)\n"
            f"    - cons_min_depth too high: a cell needs that much coverage at "
            f"a position\n"
            f"      for sam2consensus to call it at all\n"
            f"    - depth_min too high: it drops whole cells below that "
            f"coverage\n"
            f"  Lower those thresholds, or use more reads.")

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
    # genotypeName joins the per-mutation names in GENOME-POSITION order.
    #
    # The port previously joined them in the order the tokens happened to appear
    # in the genotype string, so "13A_8T" became "R5S_D3V" where the R produced
    # "D3V_R5S". Same genotype, different label. Nothing caught it because the
    # golden test compares only the per-mutation (subName, subClass) calls, not
    # this derived column. Sorting by position matches the R and reads along the
    # genome, which is what you want on a network node.
    #
    # Missing names collapse to "" rather than NaN, exactly as the old transform
    # did: reference cells have no substitutions, and genoFreq groups on this
    # column, so NaN here would silently drop those rows out of the frequency.
    _names = (merged.dropna(subset=["subName"])
                    .sort_values("pos", kind="stable")
                    .groupby("genotype")["subName"]
                    .apply(lambda s: "_".join(pd.unique(s))))
    merged["genotypeName"] = merged["genotype"].map(_names).fillna("")
    merged["genoFreq"] = merged.groupby("genotypeName")["CBC_ID"].transform("nunique") / merged["total"]
    merged["haploFreq"] = merged.groupby("genotype")["CBC_ID"].transform("nunique") / merged["total"]

    merged.to_csv(f"{out_prefix}_annot_v3.csv", index=False)
    written["annot"] = f"{out_prefix}_annot_v3.csv"

    if region_table is not None:
        region_table.to_csv(f"{out_prefix}_regionAnnotations.csv", index=False)
        written["regions"] = f"{out_prefix}_regionAnnotations.csv"

    nodes = None
    if network:
        single_steps, all_entries = hap_network_gen(haplocounts,
                                                    self_edges=self_edges)
        single_steps.to_csv(f"{out_prefix}_epistaticNetwork.csv", index=False)
        all_entries.to_csv(f"{out_prefix}_genotypeNetwork.csv", index=False)
        written["epistatic"] = f"{out_prefix}_epistaticNetwork.csv"
        written["genotype"] = f"{out_prefix}_genotypeNetwork.csv"

        # Node attributes for the two edge tables above. Written with them
        # because it is only useful alongside them.
        nodes = genotype_nodes(merged)
        nodes.to_csv(f"{out_prefix}_genotypeNodes.csv", index=False)
        written["nodes"] = f"{out_prefix}_genotypeNodes.csv"

    return {"annot": merged, "regions": region_table, "nodes": nodes,
            "written": written}
