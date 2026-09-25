# anchovy
![tests](https://github.com/QVEU/anchovy/actions/workflows/tests.yml/badge.svg)
![anchovies](assets/northern-anchovies-rw07-130.webp)
anchovy is an analysis pipeline designed for use with barcoded single-cell sequencing data to reconstruct viral haplotypes from individual cells.

It is the genotype-reconstruction half of **SEARCHLIGHT** (scRNAseq-Enabled
Acquisition of mRNA and Consensus Haplotypes Linking Individual Genotypes and
Host Transcriptomes), the method described in [Dábilla & Dolan
(2024)](https://doi.org/10.1126/sciadv.ado1693). The wet-lab side is a 10x
Genomics 5′ run with virus-specific reverse-transcription primers tiled across
the viral genome, sequenced long-read; anchovy takes those reads and returns one
consensus genome per cell, plus the genotype networks built from them. Host
transcriptomes come from the matched short-read libraries through Cell Ranger
and Seurat, which are outside this repository.

You provide a folder of sequencing reads (FASTQ, plain or gzipped) and a
reference genome. anchovy:
- maps the reads and sorts them by cell
- builds a consensus genome for each cell, and produces tables describing the
  mutations and genotypes observed
- describes how the different genotypes are related.

The entire thing can be run as one automated pipeline (recommended), or run each
step by hand.

## What the pipeline does, step by step

1. **map** — lines up the reads against the reference genome
2. **extract** — sorts the reads by the cell barcode they carry
3. **fasta** — writes out the reads for each cell into its own file
4. **map cell** — lines up each cell's reads against the reference
5. **cell consensus** — builds one consensus genome per cell from its reads
6. **merge** — collects all the per-cell genomes into one file
7. **consensus** — narrows to the region you care about and lists each cell's mutations
8. **annotate** — works out which mutations change the protein, and builds tables
   showing how genotypes relate. Given a GFF3 file describing your genome's
   regions, it can also annotate non-coding parts and number amino acids
   correctly per region (see below)
9. **frequencies** — allele frequencies over *every* mapped cell, with a
   per-position denominator, rather than only the cells that passed filtering
10. **report** — renders it all as one HTML page you can read

The tables are spreadsheet-style files (CSV) you can open in Excel or load into
other analysis tools. Two of them describe the **genotype network** — how the
different viral genotypes relate to one another — and can be opened directly in
Cytoscape to view and explore the network visually (see below).

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
just point it at a settings file.

`workflow/config.yaml` is included and is **ready to run as-is**. It points at a
small bundled dataset, so this is the quickest way to check your installation
works — it finishes in seconds and exercises every stage:

```bash
snakemake -s workflow/Snakefile --configfile workflow/config.yaml --cores 8
```

You should get six cells, two mutations (one in the polyprotein, one in the
5'UTR), a genotype network and a rendered report, in
`tests/data/fastqs/results/example/`. Open `example_report.html` and you have
seen everything the pipeline produces. Delete that directory to start over;
`tests/make_fastq_fixtures.py` documents what the bundled data contains.

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

Put your FASTQs in a folder, copy the example settings file **into that same
folder**, and point it at them. Every FASTQ in the folder is run through the
whole pipeline — mapping, barcode extraction, per-cell consensus, genotypes,
networks and a rendered report — with its outputs named after it, in
`results/<name>/` *inside that folder*.

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

### Every other setting

`workflow/config_cluster.yaml` is the fully commented reference — it sets most
of these and says why. The complete list:

| Setting | Default | What it does |
|---|---|---|
| **Input and output** | | |
| `input_dir` | — | The folder of FASTQs. Required, unless resuming with `samples` + `cells_dir` |
| `samples` | from the FASTQ names | Name the samples explicitly. Only for the `cells_dir` resume path, where there are no FASTQs to take a name from |
| `cells_dir` | `<results>/cells` | Resume from per-cell FASTAs you already have |
| `results_dir` | `<input_dir>/results` | Move the outputs |
| `resources_dir` | `<input_dir>/resources` | Move the whitelist cache |
| **Reference** | | |
| `template` | — | Reference FASTA. Required |
| `reference` | unset | Call genotypes against the genome rather than the consensus across cells. Set it to the same file as `template` |
| `reference_name` | from the FASTA header | Override the contig name |
| `gff` | unset | Region model (see below). Implies `whole_reference` |
| `whole_reference` | `false` | Keep genome coordinates without a region model |
| `minimap_preset` | `map-hifi` | `map-ont` for Nanopore |
| `map_threads` | `8` | Threads for the initial mapping |
| **Barcodes** | | |
| `chemistry` | `v3` | `v2` or `v3`. Sets the read signature *and* the whitelist together |
| `whitelist` | downloaded | A local or run-specific barcode list |
| `whitelist_url` | 10X's | Fetch the whitelist from your own mirror |
| `whitelist_barcodes` | per chemistry | Override the barcode count a download is checked against |
| `signature` | from `chemistry` | For an assay whose handles differ from 10X's |
| `max_barcode_errors` | unset | Errors a barcode may carry. Unset, every read is assigned to its nearest entry with no floor. **The stage's biggest runtime lever** — see below |
| `max_distance` | `42` | How far the *signature* match may be before a read is dropped |
| `extract_threads` | `16` | Worker pool inside `extract`. Not `--cores`. **Sets the stage's memory** — see below |
| `extract_chunk_size` | `100000` | Reads held at once. The other half of the memory setting |
| `min_reads` | `5` | Reads a barcode needs to become a cell. This decides the size of the run |
| **Calling and filtering** | | |
| `window` | unset | `cds` takes the analysis window from the GFF |
| `orf_start`, `orf_end` | — | The window explicitly, as 0-based slice bounds |
| `cons_min_depth` | `5` | Reads a position needs before a base is called |
| `cons_threshold` | `0.5` | Consensus threshold. Raising it emits *more* ambiguity codes, not fewer |
| `min_breadth` | unset | Fraction of the window a cell must have called. `1.0` means a complete coding sequence |
| `min_depth_called` | unset | Mean depth at the positions actually called |
| `depth_min` | `10` | Legacy overall filter. Superseded by `min_breadth` — see the note in `config_cluster.yaml` |
| `keep_ambiguous` | `false` | Keep IUPAC codes as genotype tokens |
| `max_mutations_per_cell` | unset | Drop cells carrying more mutations than this |
| **Output** | | |
| `allele_frequencies` | `false` | Also write the two allele-frequency tables |
| `min_alt_reads`, `min_alt_freq` | — | Thresholds for the per-cell frequency table |
| `min_cells_per_allele` | — | Threshold for the population frequency table |
| `self_edges` | `false` | Keep self-loops in the epistatic network |
| `report` | `true` | Render the HTML report |
| `report_rmd` | the bundled one | Use your own report template |

### Where everything goes

A run is a **self-contained folder**. You start with reads and a settings file;
everything the pipeline makes is written underneath, so the same command gives
the same result from any shell, and the whole run can be archived or handed to a
colleague as one directory:

```
my_experiment/                     <- input_dir
    config.yaml                    <- the settings file, kept with the data
    sample_A.fastq.gz              <- your reads; never modified
    sample_B.fastq.gz              <- add a second sample by dropping it in
    resources/                     <- the 10X whitelist, downloaded once
    results/
        sample_A/
            sample_A.sam
            cells/  work/          <- per-cell intermediates, safe to delete
            sample_A_anchovy.csv
            sample_A_allConsensus.fasta
            sample_A_filtConsensus.csv
            sample_A_annot_v3.csv
            sample_A_genotypeNetwork.csv
            sample_A_genotypeNodes.csv
            sample_A_epistaticNetwork.csv
            sample_A_alleleFrequencies.csv
            sample_A_cellAlleleFreq.csv
            sample_A_report.html   <- start here
        sample_B/
```

Nothing is ever written back over your reads. The reference, GFF and any
whitelist you supply can live anywhere — give their full paths in the settings
file — and the checkout itself stays clean, so you can delete and re-clone
anchovy without touching a run.

Two settings move things if you need them to:

| Setting | Default | When to change it |
|---|---|---|
| `results_dir` | `<input_dir>/results` | Reads are on a read-only mount, or you want outputs on a faster disk |
| `resources_dir` | `<input_dir>/resources` | Share one whitelist cache across runs — the v3 list is ~100 MB and identical every time |

Both take a full path. Snakemake also writes its own bookkeeping to
`.snakemake/` in whatever directory you launch from; that one is Snakemake's,
not anchovy's, and is safe to delete between runs.

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

Here is an example GFF3 for dengue virus 1. Copy it and edit it for your own
run — or start from `tests/data/mapping/regions.gff3`, a smaller complete file
the test suite uses.

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
anchovy extract     reads.sam whitelist.txt -o out_anchovy.csv
anchovy fasta       out_anchovy.csv cells/
anchovy consensus   allConsensus.fasta 96 10272 --out-prefix results/sample
anchovy annotate    filtConsensus.csv reference.txt results/sample
anchovy frequencies allConsensus.fasta --reference ref.fasta --out-prefix results/sample
```

`96` and `10272` are the analysis window — the first and last positions to call
variants over, as 0-based slice bounds. The mapping steps are plain `minimap2`
and are not anchovy subcommands.

To annotate by genome region (see above), the last two steps become:

```bash
anchovy consensus allConsensus.fasta 96 10272 --whole-reference --out-prefix results/sample
anchovy annotate  filtConsensus.csv reference.txt results/sample --gff regions.gff3
```

Add `--help` to any command (e.g. `anchovy extract --help`) to see its options.

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

#### If `extract` is killed

A run that stops with `died with <Signals.SIGKILL: 9>` and nothing else was
killed by the kernel for using too much memory. SIGKILL cannot be caught, so the
stage gets no chance to say so itself.

`extract` processes the SAM a chunk at a time, so **its footprint does not
depend on how big the run is** — only on the chunk and the worker count:

```
0.3 GB  +  extract_threads x 4.6 KB x extract_chunk_size
```

At the defaults (16 workers, 100,000 reads per chunk) that is about 7 GB,
whether the file holds one million reads or twenty. Lower either knob on a tight
node: `extract_threads` costs wall time in proportion, `extract_chunk_size`
costs almost nothing until the chunks get small enough that dispatch shows up.

The stage prints what it is about to need before it starts, so a subsequent kill
is at least diagnosable.

#### If `extract` is slow

Barcode assignment is about **90% of the stage's time** on real data — the
signature search is noise beside it — and nearly all of that is spent on reads
whose barcode is *not* an exact whitelist hit. So `max_barcode_errors` is the
lever, not the hardware. Measured on 40,000 reads with 44% of them inexact, 4
workers:

| `max_barcode_errors` | 20,000-barcode list | reads kept | 2,000-barcode list |
|---|---:|---:|---:|
| `0` | 1.4 s | 56% | |
| `1` | 1.9 s | 85% | 1.8 s |
| `2` | 14.4 s | 100% | 10.8 s |
| unset | 105.9 s | 100% | 8.4 s |

`2` costs roughly **7.6× what `1` does** to recover the last 15% of reads: the
2-error neighbourhood of a 16-mer is ~1,128 candidates against 48 for one error.
Whether those reads are worth the time is a judgement about your data, not a
performance question — but if a run is taking hours, start here.

Note also that the bounded search is **not** unconditionally faster than the
scan it replaces. Its cost is fixed whatever the whitelist holds, while the scan
is linear in the list — so against 10x's 737,280 entries the bound wins
enormously, but against a run-specific Cell Ranger list of a couple of thousand
barcodes, leaving it unset was *faster* than setting it to `2`. Check `wc -l` on
your whitelist first.

**More cores help; more memory and more nodes do not.** Memory is bounded (see
above) and was never the time constraint, and `extract` is a single Snakemake
job using one machine's process pool — `--executor slurm` distributes the
per-cell jobs, not this one. Raise `extract_threads` and lower
`extract_chunk_size` to pay for it.

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

## For developers

anchovy comes with an automated test suite that checks each step still produces
the expected results. To run it:

```bash
pytest                              # run all tests
pytest tests/test_annotate.py -v   # run the tests for one step
```

Test data lives in `tests/data/`. The `tests/make_*_fixtures.py` scripts
regenerate it if ever needed.

`anchovy_v3/` is the archived original this version was migrated from. It is
kept because `tests/data/golden/` is frozen from its output, so the port can be
checked against it. Nothing in the pipeline runs it.

## Good to know

- The mapping/consensus step uses **sam2consensus**, a small existing tool by
  Edgardo Ortiz
  ([original here](https://github.com/edgardomortiz/sam2consensus)). A copy,
  updated to run on modern Python, is included in `workflow/scripts/`.
- Every run ends with a rendered HTML report (`{sample}_report.html`) built from
  `visualization/anchovy_report.rmd`. Set `report: false` in your settings file
  to skip it — that is also what removes R from the requirements for a complete
  run. The CSVs are written either way, so you can make your own figures from
  them instead.
- anchovy currently assumes you're mapping against a single reference sequence.
  A segmented genome (several reference pieces) would need a small extension.

## Citation

N. Dábilla, P. T. Dolan, Structure and dynamics of enterovirus genotype networks. **Sci Adv** 10, eado1693 (2024).

The data behind that paper:

| | |
|---|---|
| Raw sequencing | GEO [GSE260709](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE260709), SRA [PRJNA1082267](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1082267) |
| Processed data and analysis code | Dryad [10.5061/dryad.6hdr7sr76](https://doi.org/10.5061/dryad.6hdr7sr76) |

anchovy also builds on two other tools, both of which should be cited if you use
it: **minimap2** (H. Li, *Bioinformatics* 34, 3094–3100, 2018) for mapping, and
**sam2consensus** (E. M. Ortiz,
[github.com/edgardomortiz/sam2consensus](https://github.com/edgardomortiz/sam2consensus))
for the per-cell consensus.
  

## License

MIT — free to use and modify. See the LICENSE file.
