# anchovy
![tests](https://github.com/QVEU/anchovy/actions/workflows/tests.yml/badge.svg)
![anchovies](assets/northern-anchovies-rw07-130.webp)
anchovy is an analysis pipeline designed for use with barcoded single-cell sequencing data to reconstruct viral haplotypes from individual cells.

You provide sequencing reads that have already been mapped to a reference genome (in the form of SAM or BAM), anchovy:
- sorts the reads by cell
- builds a consensus genome for each cell, and produces tables describing the mutations and genotypes observed
- how the different genotypes are related.

The entire thing can be run as one automated pipeline (recommended), or run each step by hand.

## What the pipeline does, step by step

1. **extract** — sorts the reads by the cell barcode they carry
2. **fasta** — writes out the reads for each cell into its own file
3. **map** — lines up each cell's reads against the reference genome
4. **cell consensus** — builds one consensus genome per cell from its reads
5. **merge** — collects all the per-cell genomes into one file
6. **consensus** — narrows to the region you care about and lists each cell's mutations
7. **annotate** — works out which mutations change the protein, and builds tables
   showing how genotypes relate. Given a GFF3 file describing your genome's
   regions, it can also annotate non-coding parts and number amino acids
   correctly per region (see below)

The final results are spreadsheet-style files (CSV) you can open in Excel or load
into other analysis tools. Two of them describe the **genotype network** — how
the different viral genotypes relate to one another — and can be opened directly
in Cytoscape to view and explore the network visually (see below).

## Setting it up

anchovy relies on a few external programs (for mapping and handling sequence
files) plus its own code. The easiest way to get everything at once is with
**conda**, an environment manager widely used in bioinformatics. If you don't
have conda, install
[Miniconda](https://docs.conda.io/en/latest/miniconda.html) first.


Install from `anchovy` GitHub repo: 
```bash
git clone https://github.com/QVEU/anchovy.git
```

Setup is two commands. The first creates a self-contained environment with all
the required programs; the second installs anchovy itself into it. 

Navigate to your `anchovy/` directory, where `environment.yml` lives and run:

```bash
conda env create -f environment.yml
conda activate anchovy
```

That first command may take a few minutes. Then install anchovy itself into it:

```bash
pip install -e .
```

Check that it worked:

```bash
anchovy --version
snakemake --version
```

**Already have an anchovy environment from an earlier version?** New releases
sometimes add programs to `environment.yml`, and an environment created before
that won't have them — you'd see a "command not found" for something the
pipeline expects. Bring yours up to date with:

```bash
conda env update -f environment.yml --prune
```

If you have pytest installed you can also run: 
```bash
pytest -v
```
to test the install. 

(Two steps rather than one: the first command installs the ready-made external programs, and the second installs anchovy's own code in a
way that lets you edit it.)

## Running the whole pipeline

You control a run through a **settings file** — a plain text file that lists
where your data is and a few options. You never need to edit anchovy's code; you
just point it at a settings file. An example, `workflow/config.yaml`, is included.

To run everything:

```bash
snakemake -s workflow/Snakefile --configfile workflow/config.yaml --cores 8
```

`--cores 8` lets it use 8 processor cores to work on multiple cells at once — set
this to however many your machine has.

**Tip:** before a real run, add `-n` to do a *preview*. This shows you exactly
what it would do without actually doing it, so you can check everything looks
right first:

```bash
snakemake -s workflow/Snakefile --configfile workflow/config.yaml --cores 8 -n
```

Another handy feature: if a run stops partway (or you change one setting), you
can just run it again and it will pick up where it left off, redoing only the
steps that actually need it — not starting over from scratch.

## Pointing it at your own data

Put your FASTQs in a folder, copy the example settings file, and point it at
them. Every FASTQ in the folder is run through the whole pipeline — mapping,
barcode extraction, per-cell consensus, genotypes, networks and a rendered
report — with its outputs named after it, in `results/<name>/`.

```yaml
input_dir: "path/to/fastqs"         # a FOLDER of reads; one run covers all of them
template: "path/to/reference.fasta" # the reference genome
gff: "path/to/reference.gff3"       # region model, for frame-correct annotation
chemistry: "v3"                     # v2 or v3 -- sets the barcode signature AND
                                    # the whitelist, which have to agree
```

That is the whole required set. The sample names come from the FASTQ filenames,
the reference name is read from the FASTA header, and the barcode whitelist is
downloaded for you — matched to `chemistry` and checked against its expected
barcode count before anything uses it.

Optional settings — analysis window, depth and breadth filters, allele
frequencies — are documented in `workflow/config.yaml`.

Plan the run first, then do it:

```bash
snakemake -s workflow/Snakefile --configfile my_settings.yaml --cores 8 -n
snakemake -s workflow/Snakefile --configfile my_settings.yaml --cores 8
```

`-n` is a dry run: it lists the jobs and the samples it found without doing any
work, which is the quickest way to catch a wrong `input_dir`. On a cluster node
with more cores, raise `--cores` to match:

```bash
snakemake -s workflow/Snakefile --configfile workflow/config_cluster.yaml --cores 64 -n
snakemake -s workflow/Snakefile --configfile workflow/config_cluster.yaml --cores 64
```

### If you already have per-cell sequence files

If you've already split your reads into one file per cell (from an earlier run or
another tool), you can skip the first steps. Add this line to your settings file
pointing at that folder:

```yaml
cells_dir: "path/to/your/per-cell-files"
```

The included `workflow/config_test.yaml` is a working example of this, set up to
run against the small test dataset that ships with anchovy.

### Annotating by genome region

By default, anchovy assumes your genome is one long protein-coding stretch that
starts at the very first base. That works for a single trimmed ORF, but it gets
two things wrong on a real viral genome: mutations in non-coding parts (like the
UTRs at either end) are treated as though they coded for protein, and amino acid
numbering starts from the wrong place whenever the coding region doesn't begin at
position 1.

You can fix both by describing your genome's regions in a **GFF3 file** — a
standard, widely used format for listing the features of a genome. Point at it
from your settings file:

```yaml
gff: "path/to/regions.gff3"
```

#### Writing your own

Here is an example GFF3 for dengue virus 1, which can be found at: `tests/data/mapping/regions.gff3`. Copy it and edit for your own run. 

```
##gff-version 3
NC_001477.1	anchovy	five_prime_UTR	1	94	.	+	.	Name=5UTR
NC_001477.1	anchovy	CDS	95	10270	.	+	0	Name=polyprotein
NC_001477.1	anchovy	mature_protein_region	95	436	.	+	0	Name=capsid
NC_001477.1	anchovy	mature_protein_region	2422	3477	.	+	0	Name=NS1
NC_001477.1	anchovy	mature_protein_region	7568	10267	.	+	0	Name=NS5
NC_001477.1	anchovy	three_prime_UTR	10271	10735	.	+	.	Name=3UTR
```

**Every column must be separated by a single TAB, not spaces.** This is a common
mistake. It will get flagged in your anchovy run, but many text editors silently convert tabs to spaces when
you type, so it's worth turning that off, or editing the file in a plain-text
editor. Copying the block above preserves the tabs.

The nine columns, left to right:

| # | Column | What to put there |
|---|--------|-------------------|
| 1 | sequence name | The name of your reference, exactly as it appears after `>` in your reference FASTA |
| 2 | source | Free text saying where the annotation came from. Put anything; `.` is fine |
| 3 | **type** | What kind of region this is — see the list below |
| 4 | **start** | First base of the region, counting the genome's first base as 1 |
| 5 | **end** | Last base of the region, included |
| 6 | score | Not used by anchovy. Put `.` |
| 7 | **strand** | `+` or `-` |
| 8 | **phase** | For coding regions, how many bases to skip before the first whole codon: `0`, `1` or `2`. Use `0` unless you know otherwise. Put `.` for non-coding regions |
| 9 | attributes | `Name=...` is what anchovy labels the region with in the results |

The types anchovy understands:

- **Coding** (mutations get amino acid annotation): `CDS`, `mature_protein_region`
- **Non-coding** (mutations reported, no amino acids): `five_prime_UTR`,
  `three_prime_UTR`, `UTR`, `stem_loop`, `ncRNA`, `misc_feature`, `region`

Any other type — `gene`, `mRNA`, `exon` and so on — is ignored, so you can leave
an annotation file from elsewhere largely as it is and anchovy will pick out the
parts it can use.

**Regions are allowed to overlap.** In the example above
the capsid, NS1 and NS5 all sit inside the polyprotein, so a mutation in NS5 gets
one row numbered within NS5 and another numbered within the polyprotein. Annotate
at whichever levels are useful to you.

**A check worth knowing about.** A coding region's length should divide exactly
by three, into codons. If it doesn't, anchovy warns you and names
the region, that almost always means the coordinates are off. The
warning doesn't stop the run, but you should double check the coordinates are correct, genomic nucleotide positions. 

If you'd rather start from something known to work, `tests/data/mapping/regions.gff3`
in this repository is a small, complete file used by the test suite.

With that in place you get an extra results file,
`{sample}_regionAnnotations.csv`, with one row per mutation **per region it falls
in**. So a mutation inside NS5, which also sits inside the polyprotein, gets two
rows — numbered correctly in each one (residue 2 of NS5 is also residue 2493 of
the polyprotein). Mutations in the UTRs get rows too, simply with the amino acid
columns left empty, and anything outside every listed region is labelled
`intergenic`.

Positions in this file are always counted from the start of the whole genome, so
they mean the same thing no matter which region a row is about.

A few things worth knowing:

- Adding `gff` also keeps the **whole** reference genome through the pipeline
  instead of trimming it, because region lookup needs real genome positions. When
  you do that, `orf_start` and `orf_end` stop trimming anything and instead just
  mark the stretch that's well covered by your reads, so ragged ends don't get
  called as mutations.
- Your existing results files are still produced exactly as before, so anything
  you already do with them keeps working.
- Both strands are handled. Which strand a region is on, and its reading frame,
  are taken from columns 7 and 8 of the GFF3 itself, so anchovy never has to
  guess. For a region on the `-` strand, amino acids are numbered from the
  region's END and read in the reverse-complement direction, while the
  nucleotide change is still reported in ordinary forward-strand genome terms.
- If you want genome positions but don't have a GFF3, use `whole_reference: true`
  on its own.

### A complete worked example

`examples/eva71_sra/` runs the whole pipeline on real public data, starting from
nothing but an SRA accession: it downloads an enterovirus A71 single-cell run and
its reference genome, builds the region file from the reference's own GenBank
annotation, and hands everything to the workflow. It's written to be copied and
pointed at your own virus. See the README in that folder.

## Running one step at a time

If you'd rather run steps individually instead of the full pipeline, each is its
own command:

```bash
anchovy extract   reads.sam whitelist.txt -o out_anchovy.csv
anchovy fasta     out_anchovy.csv cells/
anchovy consensus allConsensus.fasta 96 10272 --out-prefix results/sample
anchovy annotate  filtConsensus.csv reference.txt results/sample
```

To annotate by genome region (see above), the last two steps become:

```bash
anchovy consensus allConsensus.fasta 96 10272 --whole-reference --out-prefix results/sample
anchovy annotate  filtConsensus.csv reference.txt results/sample --gff regions.gff3
```

Add `--help` to any command (e.g. `anchovy extract --help`) to see its options.

## Viewing the genotype networks

Two of anchovy's output files describe how the viral genotypes relate to each
other as a network:

- `<sample>_genotypeNetwork.csv` — all the relationships between genotypes
- `<sample>_epistaticNetwork.csv` — just the "single-step" links (genotypes that
  differ by exactly one mutation), plus links back to the reference
- `<sample>_genotypeNodes.csv` — one row per genotype, describing the genotypes
  themselves rather than the links between them

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
| `genotypeName` | The amino acid change(s), like `R5S` or `D3V_R5S`. If you annotated with a GFF3, non-coding changes appear here too, like `5UTR:A121C` |
| `nMutations` | How many mutations the genotype carries |
| `nCells` | How many cells carry it |
| `genoFreq`, `haploFreq` | What fraction of cells that is |

Once it's imported, the useful moves in Cytoscape's **Style** panel are:

- Set node **Label** to `genotypeName`, so each point is named by the amino acid
  change rather than an internal identifier.
- Map node **Size** to `nCells` or `genoFreq` (continuous mapping), so common
  genotypes are visibly bigger.
- Map node **Fill Color** to `nMutations` (continuous mapping), so how far a
  genotype has drifted from the reference reads at a glance.

**A note on the two network files.** In the epistatic network, links from a
genotype to itself are left out — they'd draw as a small loop on every point and
tell you nothing. The genotype network keeps them, because there they do carry
information: a genotype that shares no mutation with any other appears *only* as
its own self-link, so removing them would make it vanish from the picture
entirely. If you want the original R behavior in both, add `self_edges: true` to
your settings file.

## For developers

anchovy comes with an automated test suite that checks each step still produces
the expected results. To run it:

```bash
pytest                              # run all tests
pytest tests/test_annotate.py -v   # run the tests for one step
```

Test data lives in `tests/data/`. The `tests/make_*_fixtures.py` scripts
regenerate it if ever needed.

## Good to know

- The mapping/consensus step uses **sam2consensus**, a small existing tool by
  Edgardo Ortiz
  ([original here](https://github.com/edgardomortiz/sam2consensus)). A copy,
  updated to run on modern Python, is included in `workflow/scripts/`.
- The original version of this analysis also made plots. This version produces the
  data tables only; you can make figures from those in your tool of choice.
- anchovy currently assumes you're mapping against a single reference sequence.
  A segmented genome (several reference pieces) would need a small extension.

## License

MIT — free to use and modify. See the LICENSE file.
