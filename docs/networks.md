# Interpreting genotype networks

[Back to the README](../README.md)

## Viewing the genotype networks

Two of anchovy's output files describe how the viral genotypes relate to each
other as a network:

- `<sample>_epistaticNetwork.csv` — **"single-step" links**: pairs of genotypes
  differing by exactly one mutation, plus links back to the reference
- `<sample>_genotypeNetwork.csv` — **shared-mutation links**: every pair of
  genotypes with at least one mutation in common, however far apart they are
- `<sample>_genotypeNodes.csv` — one row per genotype, describing the genotypes
  themselves rather than the links between them

**Which one you want is probably `_epistaticNetwork.csv`, despite the names.**
The networks in Dábilla & Dolan (2024) — where "edges represent single-nucleotide
substitutions linking individual genotypes" — are the *single-step* network.
`_genotypeNetwork.csv` joins any two genotypes sharing a mutation, so it is much
denser and its edges do not mean one mutational step. The file names are kept
for compatibility with the original R output.

You can explore these visually in **Cytoscape**, a free tool for viewing and
analyzing networks that's widely used in biology. Download it from
[cytoscape.org](https://cytoscape.org/), or use the no-install browser version at
[web.cytoscape.org](https://web.cytoscape.org/).

To load an anchovy network:

1. In Cytoscape, choose **File → Import → Network from File** and pick one of the
   network CSVs.
2. Cytoscape will recognize the `source` and `target` columns automatically and
   draw the network — no manual setup needed. The other columns (`overlap`,
   `mutNumSource`, `count`, and so on) come in as **edge attributes**, properties
   of each link you can use for styling: for example, make links thicker when
   more mutations are shared.
3. Then choose **File → Import → Table from File** and pick
   `<sample>_genotypeNodes.csv`. Cytoscape matches its `genotype` column to the
   genotypes already in your network and attaches the rest as **node
   attributes**.

Each point (node) in the resulting picture is a genotype; each line (edge) is a
relationship between two genotypes. This is a quick way to see, for instance,
which mutations tend to build on one another.

### Making the picture readable

That third file is what turns the network from a set of unlabeled dots into
something you can interpret. Each row describes one genotype:

| Column | What it is |
|--------|------------|
| `genotype` | The genotype's identifier — this is what Cytoscape matches on |
| `genotypeID` | The amino acid change(s), like `R5S` or `D3V_R5S`. If you annotated with a GFF3, non-coding changes appear here too, like `5UTR:A121C` |
| `nMutations` | How many mutations the genotype carries |
| `nCells` | How many cells carry it |
| `genoFreq` | What fraction of cells that is |
| `idFreq` | What fraction of cells carry *any* genotype with this `genotypeID` |

**`genotype` and `genotypeID` are not the same thing, and neither are their
frequencies.** `genotype` is the nucleotide haplotype and is what a node *is* —
one node per distinct nucleotide sequence. `genotypeID` is its translation, and
the mapping is many-to-one: a synonymous change is written `X_n_X`, so several
distinct genotypes can share one ID. Where that happens, `genoFreq` counts the
cells carrying that one nucleotide genotype and `idFreq` counts every cell whose
genotype translates the same way — for three synonymous genotypes in one cell
each out of four, `genoFreq` is 0.25 and `idFreq` is 0.75.

**Size nodes on `genoFreq` or `nCells`, not `idFreq`**, since a node is one
nucleotide genotype; `idFreq` would size each member of a synonymous group by
the whole group.

Once it's imported, the useful moves in Cytoscape's **Style** panel are:

- Set node **Label** to `genotypeID`, so each point is named by the amino acid
  change rather than by its nucleotide haplotype.
- Map node **Size** to `nCells` or `genoFreq` (continuous mapping), so common
  genotypes are visibly bigger. The two are the same quantity, counted and as a
  fraction.
- Map node **Fill Color** to `nMutations` (continuous mapping), so how far a
  genotype has drifted from the reference reads at a glance.

**A note on the two network files.** In the epistatic network, links from a
genotype to itself are left out — they'd draw as a small loop on every point and
tell you nothing. The genotype network keeps them, because there they do carry
information: a genotype that shares no mutation with any other appears *only* as
its own self-link, so removing them would make it vanish from the picture
entirely. If you want the original R behavior in both, add `self_edges: true` to
your settings file.
