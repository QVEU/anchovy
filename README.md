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

Setup is two commands. The first creates a self-contained environment with all
the required programs; the second installs anchovy itself into it. 

Navigate to your `anchovy/` directory, where `environment.yml` lives and run:

```bash
# 1. create the environment (this may take a few minutes)
conda env create -f environment.yml
conda activate anchovy

# 2. install anchovy
pip install -e .
```

Check that it worked:

```bash
anchovy --version
snakemake --version
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

Copy the example settings file and change the paths to match your data:

```yaml
sample: "my_sample"                 # your reads should be at {data_dir}/{sample}.sam
data_dir: "path/to/data"
whitelist: "path/to/10x_whitelist.txt"    # your list of valid cell barcodes
signature: "CTACACGACGCTCTTCCGATCT..."    # the barcode signature for your 10X kit
template: "path/to/reference.fasta"       # the reference genome
reference_name: "your_ref_name"           # the name written after ">" in that file
orf_start: 96                             # start of the region you want analysed
orf_end: 10272                            # end of that region
```

Then run with your file instead of the example:

```bash
snakemake -s workflow/Snakefile --configfile my_settings.yaml --cores 8
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

A small one looks like this (columns are separated by tabs):

```
##gff-version 3
dengue	.	five_prime_UTR	1	96	.	+	.	Name=5UTR
dengue	.	CDS	97	10272	.	+	0	Name=polyprotein
dengue	.	mature_protein_region	7570	10269	.	+	0	Name=NS5
dengue	.	three_prime_UTR	10273	10727	.	+	.	Name=3UTR
```

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
  are taken from the GFF3 file itself, so anchovy never has to guess.
- If you want genome positions but don't have a GFF3, use `whole_reference: true`
  on its own.

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
   more mutations are shared, or color nodes by how many mutations a genotype
   carries.

Each point (node) in the resulting picture is a genotype; each line (edge) is a
relationship between two genotypes. This is a quick way to see, for instance,
which mutations tend to build on one another.

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
