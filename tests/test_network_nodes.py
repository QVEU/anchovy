"""
test_network_nodes.py -- the genotype NODE table that accompanies the network CSVs.

WHY IT EXISTS. Both network files are entirely edge-level (source, target,
overlap, mutNumSource, mutNumTarget). Cytoscape therefore draws nodes carrying
no attributes, so a genotype cannot be sized by how common it is or labelled by
what it does to the protein without hand-joining CSVs inside Cytoscape -- even
though every one of those values is already computed.

The contract that matters is that `genotype` here matches the source/target
values in the network CSVs exactly, so a node-table import keys straight onto
them. Most of these tests are about that join holding.
"""

from __future__ import annotations

import pandas as pd
import pytest

from anchovy.annotate import NODE_COLUMNS, genotype_nodes, run


@pytest.fixture
def annot_fixture(request):
    data = request.path.parent / "data" / "annot"
    for name in ("filtConsensus.csv", "reference.txt"):
        if not (data / name).exists():
            pytest.skip(f"missing {name}; run tests/make_annot_fixtures.py first.")
    return data


# --------------------------------------------------------------------------- #
# The join onto the network files
# --------------------------------------------------------------------------- #
def test_node_ids_match_the_network_edge_ids_exactly(annot_fixture, tmp_path):
    """Every node id appears in the network, and every network id has a node.

    This is the whole point of the file: if the two sides disagree, Cytoscape
    silently imports attributes for nothing, or leaves real nodes bare.
    """
    out = str(tmp_path / "x")
    result = run(str(annot_fixture / "filtConsensus.csv"),
                 str(annot_fixture / "reference.txt"), out, network=True)

    nodes = set(result["nodes"]["genotype"])
    for network_file in ("genotypeNetwork", "epistaticNetwork"):
        edges = pd.read_csv(f"{out}_{network_file}.csv")
        edge_ids = set(edges["source"]) | set(edges["target"])
        assert edge_ids == nodes, (
            f"{network_file}: node/edge id mismatch\n"
            f"  only in edges: {edge_ids - nodes}\n"
            f"  only in nodes: {nodes - edge_ids}")


def test_reference_is_a_node(annot_fixture, tmp_path):
    """The reference is a real node in the network, so it needs a row too."""
    result = run(str(annot_fixture / "filtConsensus.csv"),
                 str(annot_fixture / "reference.txt"),
                 str(tmp_path / "x"), network=True)
    ref = result["nodes"].set_index("genotype").loc["reference"]
    # It has no substitutions, so it gets a label saying so rather than a blank.
    assert ref["genotypeName"] == "reference"
    assert ref["nMutations"] == 0


def test_one_row_per_genotype(annot_fixture, tmp_path):
    result = run(str(annot_fixture / "filtConsensus.csv"),
                 str(annot_fixture / "reference.txt"),
                 str(tmp_path / "x"), network=True)
    nodes = result["nodes"]
    assert not nodes["genotype"].duplicated().any()


# --------------------------------------------------------------------------- #
# The attributes themselves
# --------------------------------------------------------------------------- #
def test_node_attributes_are_the_ones_worth_styling_by(annot_fixture, tmp_path):
    """Amino-acid label, mutation count, cell count and both frequencies."""
    result = run(str(annot_fixture / "filtConsensus.csv"),
                 str(annot_fixture / "reference.txt"),
                 str(tmp_path / "x"), network=True)
    nodes = result["nodes"].set_index("genotype")
    assert list(result["nodes"].columns) == NODE_COLUMNS

    # Labels are amino-acid level, not bare nucleotide tokens. This is what
    # makes a rendered network readable.
    assert nodes.loc["13A", "genotypeName"] == "R5S"
    assert nodes.loc["13A_8T", "genotypeName"] == "D3V_R5S"

    assert nodes.loc["13A", "nMutations"] == 1
    assert nodes.loc["13A_8T", "nMutations"] == 2
    assert nodes.loc["13A", "nCells"] == 2
    assert nodes.loc["13A_8T", "nCells"] == 1


def test_nmutations_agrees_with_the_networks_mutnum(annot_fixture, tmp_path):
    """nMutations must equal the edge tables' mutNumSource for the same node.

    Two files describing the same genotype disagreeing about how many mutations
    it carries would be worse than not shipping the column.
    """
    out = str(tmp_path / "x")
    result = run(str(annot_fixture / "filtConsensus.csv"),
                 str(annot_fixture / "reference.txt"), out, network=True)
    nodes = result["nodes"].set_index("genotype")["nMutations"]
    edges = pd.read_csv(f"{out}_genotypeNetwork.csv")

    for _, e in edges.iterrows():
        if e["source"] == "reference":
            continue           # the network counts "reference" as one token
        assert nodes[e["source"]] == e["mutNumSource"], e["source"]


def test_ncells_matches_the_networks_count_column(annot_fixture, tmp_path):
    """The network's `count` is the number of cells carrying that genotype.

    Its computation in hap_network_gen is a faithfully-ported R expression,
    rowSums(countMatrix)/rowSums(binaryMatrix), which reads like a ratio but
    reduces to the cell count: every cell with genotype g contributes one row
    per mutation in g, so it is (n*k)/k = n. Pinning the two together here means
    the roundabout version cannot drift away from the plain one unnoticed.
    """
    out = str(tmp_path / "x")
    result = run(str(annot_fixture / "filtConsensus.csv"),
                 str(annot_fixture / "reference.txt"), out, network=True)
    nodes = result["nodes"].set_index("genotype")["nCells"]
    edges = pd.read_csv(f"{out}_genotypeNetwork.csv")

    for _, e in edges[["source", "count"]].drop_duplicates().iterrows():
        assert nodes[e["source"]] == e["count"], (
            f"{e['source']}: node table says {nodes[e['source']]} cells, "
            f"network `count` says {e['count']}")


def test_frequencies_come_through(annot_fixture, tmp_path):
    result = run(str(annot_fixture / "filtConsensus.csv"),
                 str(annot_fixture / "reference.txt"),
                 str(tmp_path / "x"), network=True)
    nodes = result["nodes"].set_index("genotype")
    # 7 cells in the fixture; 13A is carried by 2 of them.
    assert nodes.loc["13A", "genoFreq"] == pytest.approx(2 / 7)
    assert nodes.loc["13A", "haploFreq"] == pytest.approx(2 / 7)


# --------------------------------------------------------------------------- #
# Wiring
# --------------------------------------------------------------------------- #
def test_file_is_written_with_the_networks(annot_fixture, tmp_path):
    out = str(tmp_path / "x")
    result = run(str(annot_fixture / "filtConsensus.csv"),
                 str(annot_fixture / "reference.txt"), out, network=True)
    assert result["written"]["nodes"] == f"{out}_genotypeNodes.csv"

    written = pd.read_csv(f"{out}_genotypeNodes.csv")
    assert list(written.columns) == NODE_COLUMNS
    assert len(written) == len(result["nodes"])


def test_not_written_when_networks_are_off(annot_fixture, tmp_path):
    """It is node attributes FOR the edge tables, so it follows them."""
    out = str(tmp_path / "x")
    result = run(str(annot_fixture / "filtConsensus.csv"),
                 str(annot_fixture / "reference.txt"), out, network=False)
    assert result["nodes"] is None
    assert "nodes" not in result["written"]
    assert not (tmp_path / "x_genotypeNodes.csv").exists()


def test_region_aware_names_reach_the_nodes(tmp_path):
    """With a GFF, node labels carry the region-aware names.

    This is the payoff of region awareness for visualization: a node can read
    "5UTR:A121C" or a frame-correct residue, rather than a nucleotide token.
    """
    reference = "ACGT" * 50
    ref_file = tmp_path / "reference.txt"
    ref_file.write_text(reference)
    gff = tmp_path / "regions.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "r\ta\tfive_prime_UTR\t1\t60\t.\t+\t.\tName=5UTR\n"
        "r\ta\tCDS\t61\t180\t.\t+\t0\tName=polyprotein\n")

    cons = tmp_path / "filtConsensus.csv"
    pd.DataFrame([
        {"CBC_ID": "c1", "genotype": "10G", "sequence": reference,
         "description": "coverage:50"},
        {"CBC_ID": "c2", "genotype": "", "sequence": reference,
         "description": "coverage:50"},
    ]).to_csv(cons, index=False)

    result = run(str(cons), str(ref_file), str(tmp_path / "out"),
                 network=True, gff=str(gff))
    nodes = result["nodes"].set_index("genotype")
    assert nodes.loc["10G", "genotypeName"] == "5UTR:C10G"


def test_empty_input_gives_an_empty_table():
    empty = pd.DataFrame(columns=["genotype", "genotypeName", "CBC_ID",
                                  "genoFreq", "haploFreq"])
    out = genotype_nodes(empty)
    assert out.empty and list(out.columns) == NODE_COLUMNS


# --------------------------------------------------------------------------- #
# Nodes the edges name but no cell carries
# --------------------------------------------------------------------------- #
# hap_network_gen re-points every single-mutation self-step at "reference", so
# "reference" is an edge endpoint whenever any cell carries exactly one
# mutation -- whether or not any cell IS wild-type. Those come apart on a
# passaged population with a mutation fixed relative to the supplied genome:
# nothing is wild-type, yet the hub every one-step edge points at is still
# "reference". Cytoscape then drew that hub unlabelled and unsized.
def _cons(tmp_path, rows):
    import pandas as pd
    p = tmp_path / "filtConsensus.csv"
    pd.DataFrame(rows).to_csv(p, index=False)
    return p


def test_reference_node_exists_when_no_cell_is_wild_type(tmp_path):
    import pandas as pd
    from anchovy import annotate

    ref = tmp_path / "reference.txt"
    ref.write_text("ACGTACGTACGT")
    cons = _cons(tmp_path, [
        {"CBC_ID": "c2", "genotype": "5T", "sequence": "", "description": "coverage:50"},
        {"CBC_ID": "c3", "genotype": "9A", "sequence": "", "description": "coverage:50"},
    ])
    result = annotate.run(str(cons), str(ref), str(tmp_path / "out"), network=True)

    nodes = result["nodes"].set_index("genotype")
    assert "reference" in nodes.index
    # No cell carries it, and the frequencies say so rather than going blank --
    # Cytoscape sizes on them and a missing value is not a small node.
    assert nodes.loc["reference", "nCells"] == 0
    assert nodes.loc["reference", "genoFreq"] == 0.0
    assert nodes.loc["reference", "nMutations"] == 0
    assert nodes.loc["reference", "genotypeName"] == "reference"


def test_every_edge_endpoint_has_a_node_row(tmp_path):
    """The invariant the fix is really about, checked against the edge files."""
    import pandas as pd
    from anchovy import annotate

    ref = tmp_path / "reference.txt"
    ref.write_text("ACGTACGTACGT")
    cons = _cons(tmp_path, [
        {"CBC_ID": "c1", "genotype": "5T", "sequence": "", "description": "coverage:50"},
        {"CBC_ID": "c2", "genotype": "5T_9A", "sequence": "", "description": "coverage:50"},
        {"CBC_ID": "c3", "genotype": "9A", "sequence": "", "description": "coverage:50"},
    ])
    out = tmp_path / "out"
    annotate.run(str(cons), str(ref), str(out), network=True)

    nodes = set(pd.read_csv(f"{out}_genotypeNodes.csv")["genotype"].astype(str))
    for name in ("epistaticNetwork", "genotypeNetwork"):
        edges = pd.read_csv(f"{out}_{name}.csv")
        endpoints = set(edges["source"].astype(str)) | set(edges["target"].astype(str))
        assert endpoints <= nodes, f"{name} names nodes absent from the node table"


def test_wild_type_cells_still_count_toward_the_reference_node(tmp_path):
    # The fill-in must not shadow a real reference row: when cells ARE
    # wild-type, nCells is their count, not 0.
    import pandas as pd
    from anchovy import annotate

    ref = tmp_path / "reference.txt"
    ref.write_text("ACGTACGTACGT")
    cons = _cons(tmp_path, [
        {"CBC_ID": "c1", "genotype": "", "sequence": "", "description": "coverage:50"},
        {"CBC_ID": "c2", "genotype": "5T", "sequence": "", "description": "coverage:50"},
    ])
    result = annotate.run(str(cons), str(ref), str(tmp_path / "out"), network=True)
    nodes = result["nodes"].set_index("genotype")
    assert nodes.loc["reference", "nCells"] == 1
    assert nodes.loc["reference", "genoFreq"] > 0


# --------------------------------------------------------------------------- #
# Self-edges: dropping them must not delete a node
# --------------------------------------------------------------------------- #
# A genotype with no single-step neighbour appears in the epistatic network
# ONLY as its own self-edge, so stripping self-edges removes the genotype from
# the file rather than tidying it. The docstring used to claim otherwise,
# generalising from the fully-connected fixture. On the EV-A71 run it cost 69
# of 300 genotypes carrying 74 of 544 cells -- the isolated lineages.
def _network(tmp_path, genotypes, **kw):
    import pandas as pd
    from anchovy import annotate

    ref = tmp_path / "reference.txt"
    # Length a multiple of 3: a partial trailing codon makes Biopython warn on
    # every translate, which is noise rather than signal for a network test.
    ref.write_text("ACG" * 100)
    cons = tmp_path / "filtConsensus.csv"
    pd.DataFrame([{"CBC_ID": f"c{i}", "genotype": g, "sequence": "",
                   "description": "coverage:50"}
                  for i, g in enumerate(genotypes)]).to_csv(cons, index=False)
    out = tmp_path / "out"
    annotate.run(str(cons), str(ref), str(out), network=True, **kw)
    ep = pd.read_csv(f"{out}_epistaticNetwork.csv")
    nodes = pd.read_csv(f"{out}_genotypeNodes.csv")
    return ep, nodes


# "5T" and "5T_9A" are one step apart; "100G_200C" shares nothing with either.
ISOLATED_CASE = ["5T", "5T_9A", "100G_200C"]


def test_isolated_genotype_keeps_its_self_edge(tmp_path):
    ep, _ = _network(tmp_path, ISOLATED_CASE)
    endpoints = set(ep["source"].astype(str)) | set(ep["target"].astype(str))
    assert "100G_200C" in endpoints, (
        "a genotype with no single-step neighbour was dropped from the "
        "epistatic network along with its self-edge")


def test_connected_genotypes_still_lose_their_self_edges(tmp_path):
    # The clutter the default exists to avoid must still be avoided: only the
    # node that would otherwise vanish keeps its loop.
    ep, _ = _network(tmp_path, ISOLATED_CASE)
    loops = set(ep[ep["source"] == ep["target"]]["source"].astype(str))
    assert loops == {"100G_200C"}


def test_node_table_joins_completely_against_the_epistatic_network(tmp_path):
    """The property the Cytoscape workflow actually depends on."""
    ep, nodes = _network(tmp_path, ISOLATED_CASE)
    endpoints = set(ep["source"].astype(str)) | set(ep["target"].astype(str))
    assert set(nodes["genotype"].astype(str)) == endpoints


def test_self_edges_true_keeps_every_loop(tmp_path):
    # Unchanged: the opt-in still reproduces the R, loops on everything.
    ep, _ = _network(tmp_path, ISOLATED_CASE, self_edges=True)
    loops = set(ep[ep["source"] == ep["target"]]["source"].astype(str))
    assert {"5T", "5T_9A", "100G_200C"} <= loops
