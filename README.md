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

## Quickstart

Install [Miniforge](https://conda-forge.org/download/) if you do not already have
Conda. Run these commands in a shell where `conda` is available:

```bash
git clone https://github.com/QVEU/anchovy.git
cd anchovy
conda env create -f environment.yml
conda activate anchovy
python -m pip install -e . --no-deps
anchovy --version
snakemake --version
snakemake -s workflow/Snakefile --configfile workflow/config.yaml --cores 2
```

This installs Python 3.12, minimap2, samtools, Snakemake, R, the report packages,
and the developer test tools. Conda supplies dependencies; the editable pip
install connects the `anchovy` command to this checkout. Allow a few minutes
for environment creation.

The bundled example needs no downloads or credentials. It should produce:

- **Six cells and four genotypes:** `reference`, `121C`, `201A`, `121C_201A`.
- **Two variant sites:** one in the 5′UTR and one in the polyprotein.
- **An HTML report:** `tests/data/fastqs/results/example/example_report.html`.
- Consensus sequences, annotations, genotype networks and allele-frequency tables
  beside the report.

Open the HTML report to inspect the results. `--cores 2` fits a small machine;
raise it for larger runs after checking [memory requirements](docs/performance.md).
Running the command again resumes the workflow and leaves up-to-date outputs
alone. To rerun from scratch without deleting results, choose a new output path:

```bash
snakemake -s workflow/Snakefile --configfile workflow/config.yaml --cores 2 \
  --config results_dir=/absolute/path/to/new-results
```

For exact Linux x86_64 dependencies, use
`conda create -n anchovy --file locks/linux-64/explicit.txt` **instead of** the
`conda env create` command above. Use `--no-build-isolation` on the pip install
to use the locked build tools. See [reproducibility](docs/development.md#dependency-environments)
for platform limits and lock-refresh instructions.

To check the installation, run `pytest -v`. The suite includes the bundled
workflow and report rendering; CI requires their external tools.

## Use your own data

Inputs are long-read FASTQs (plain or gzipped), a closely matched single-sequence
reference FASTA, and preferably a GFF3 region model. Each FASTQ filename identifies
one sample. anchovy currently assumes a single reference sequence; segmented
genomes need additional support.

Create a settings file such as `my_settings.yaml`. Use absolute paths so its
meaning is independent of where you store the file:

```yaml
input_dir: "/data/experiment/fastqs"
template: "/data/experiment/reference.fasta"
reference: "/data/experiment/reference.fasta"
gff: "/data/experiment/regions.gff3"
chemistry: "5p-v3"  # v2, v3 or 5p-v3; must match the library
```

`reference` calls genotypes against that genome. If omitted, genotypes are called
against the consensus across cells, which can absorb population-wide variants.
The reference must closely match the sequenced virus: divergence is otherwise
reported as mutations in every cell and can obscure variation within the sample.
A consensus from the infected stock is a useful starting point.

The GFF3 supplies coding frames and non-coding regions; see the
[GFF3 guide and tab-separated example](docs/usage.md#annotating-by-genome-region).
Without a GFF3, the default assumes one coding region beginning at base one.

A chemistry-matched 10X barcode whitelist downloads automatically unless you
supply `whitelist`. Run-specific lists of called cells are also supported, either
one shared by all samples or one per sample; see
[whitelists and multiplexed samples](#whitelists-per-sample-lists-and-multiplexed-samples)
below. When
adapting `workflow/config.yaml`, replace its bundled `whitelist`, `template`,
`reference`, GFF and analysis-window settings; those values belong to the example.

### Whitelists, per-sample lists and multiplexed samples

The whitelist is the set of cell barcodes that reads may be assigned to. There
are three ways to supply one, and which is right depends on how the PacBio
library relates to the Cell Ranger run.

**The full 10x list (default).** Leave `whitelist` unset and the list matching
`chemistry` downloads once. Every barcode in the library can be matched,
including background barcodes outside called cells, so expect more cells than
Cell Ranger called. Restrict to Cell Ranger's cells afterward by barcode.

**One shared list for all samples.** `whitelist: "/path/cells.txt"`, used for
every FASTQ in `input_dir`, with one 16 nt
barcode per line. Cell Ranger writes barcodes with a `-1` GEM-well suffix;
strip it (`zcat barcodes.tsv.gz | cut -d- -f1 > cells.txt`).

**A list per sample.** When each FASTQ has its own called cells, map sample
names to lists. A sample's name is its FASTQ filename without `.fastq`, `.fq`
and `.gz`; every sample must be listed, and the run stops before any job if
one is missing.

```yaml
whitelist:
  P10_1: "/data/cellranger/p10_1/barcodes.txt"
  P10_2: "/data/cellranger/p10_2/barcodes.txt"
```

#### Multiplexed samples: one list per library, not per sample

A list must cover the cells **in the library being sequenced**, not the cells
of one sample. The distinction matters whenever several samples were pooled
into one 10x channel (one GEM well), using hashtag antibodies, CellPlex, on-chip
multiplexing or similar, and demultiplexed by `cellranger multi`. Its output
then has a `per_sample_outs/<sample>/` folder for each sample.

In that case:

- Every cell in the channel gets a 10x barcode, whichever sample it came from.
  The long-read library made from that channel therefore holds reads from all
  of its samples' cells.
- `per_sample_outs/<sample>/.../barcodes.tsv.gz` holds only that sample's
  cells. Used as the whitelist for the pooled library, it covers a fraction of
  the reads: with four samples, roughly a quarter at best.
- `extract` checks that at least `min_exact_rate` (default 5%) of reads carry
  a barcode that is exactly in the list, because a list from the wrong kit
  looks the same: a run that matches almost nothing. A single sample's list
  on a pooled library can fall below that floor with nothing wrong in the
  data, and the stage stops with "only X% of the first N reads carry a barcode
  that is an EXACT entry in this whitelist".

So treat each **channel** as one anchovy sample, give it the combined list for
every sample in the channel, and split the results by sample afterward:

```bash
# One whitelist for the channel: all samples' called cells, suffix stripped.
for d in cellranger/p10_5-8/outs/per_sample_outs/*/; do
  zcat "$d"count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz
done | cut -d- -f1 | sort -u > p10_5-8_cells.txt

# A barcode -> sample table, to split anchovy's per-cell outputs (CBC_ID) later.
for d in cellranger/p10_5-8/outs/per_sample_outs/*/; do
  s=$(basename "$d")
  zcat "$d"count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz \
    | cut -d- -f1 | sed "s/$/,$s/"
done > p10_5-8_barcode_to_sample.csv

# Barcodes in more than one sample cannot be split by barcode; check for any.
cut -d, -f1 p10_5-8_barcode_to_sample.csv | sort | uniq -d | wc -l
```

Then name the channel's FASTQ for the channel (for example `p10_5-8.fastq.gz`)
and point `whitelist` at `p10_5-8_cells.txt`, or at the full 10x list. Joining
`_filtConsensus.csv` or `_annot_v3.csv` on `CBC_ID` against the table assigns
each cell to its sample.

#### On-Chip Multiplexing (OCM)

With 10x On-Chip Multiplexing (GEM-X Universal 3′ v4 or 5′ v3 4-plex), each
sample's inlet is written into the cell barcode itself, so anchovy can split
the channel FASTQ into its samples without any join. Name the channel FASTQ
and each sample's inlet, copied from the `cellranger multi` config:

```yaml
ocm:
  p10_5-8: {P10_5: OB1, P10_6: OB2, P10_7: OB3, P10_8: OB4}
```

The channel is mapped and extracted once, and every sample then gets its own
results under `results/<sample>/`. The
[On-Chip Multiplexing guide](docs/ocm.md) covers how to tell whether a run is
OCM, the built-in OB table and checking it with `anchovy ocm-table`, naming, what runs, and
troubleshooting.

Only if a library genuinely contains just some of the channel's cells, or you
deliberately want one sample's cells out of a pooled library, lower the floor
with `min_exact_rate: 0.01` (or `0` to turn the check off). Reads from other
samples' cells are then dropped when `max_barcode_errors` is set, and assigned
to the nearest listed cell when it is not, so set it.

Run from the Anchovy checkout, using your settings file:

```bash
snakemake -s workflow/Snakefile --configfile my_settings.yaml --cores 2 -n
snakemake -s workflow/Snakefile --configfile my_settings.yaml --cores 2 --keep-going
```

The first command previews the jobs. The second runs them; `--keep-going` allows
healthy samples to finish if another sample fails. Set filter thresholds and the
analysis window for your assay using the [configuration reference](docs/usage.md).

By default, results go to `<input_dir>/results/<sample>/` and downloaded
whitelists to `<input_dir>/resources/`. `results_dir` and `resources_dir` override
those locations. Input reads are never overwritten. Snakemake keeps bookkeeping
in `.snakemake/` in the directory where you launch it.

## Read the results

Start with `<sample>_report.html`, which summarizes coverage, filtering, mutations,
and genotype networks. Reports use R Markdown and Pandoc; set `report: false` to
skip reporting and run the remaining pipeline without R.

| Output | Meaning |
|---|---|
| `_allConsensus.fasta` | Per-cell consensus sequences before genotype filtering |
| `_filtConsensus.csv` | Cells retained for genotype analysis |
| `_annot_v3.csv` | Annotated mutations; a cell with several mutations has several rows |
| `_regionAnnotations.csv` | One row per variant per overlapping GFF3 region, when a GFF is supplied |
| `_genotypeNodes.csv` | One row per nucleotide genotype, with cell counts and frequencies |
| `_epistaticNetwork.csv` | Single-step links between genotypes, including links to the reference |
| `_epistaticNetwork.cyjs` | The same single-step network with node attributes attached, for a one-step Cytoscape import |
| `_genotypeNetwork.csv` | Links between genotypes sharing mutations; these need not be single steps |
| `_alleleFrequencies.csv`, `_cellAlleleFreq.csv` | Population and per-cell frequencies when `allele_frequencies: true` |

`genotype` identifies the nucleotide haplotype. `genotypeID` is its annotation
label, which several synonymous genotypes can share. Use `nCells` or `genoFreq`
to size network nodes; `idFreq` pools cells sharing an annotation label.
See [network interpretation and Cytoscape import](docs/networks.md) for details.

## Guides

- [Configuration, output layout, GFF3 and individual commands](docs/usage.md)
- [Memory, barcode matching and performance troubleshooting](docs/performance.md)
- [Cluster allocations and SLURM execution](docs/cluster.md)
- [On-Chip Multiplexing (OCM) samples](docs/ocm.md)
- [Genotype networks and Cytoscape](docs/networks.md)
- [Development, legacy fixtures and dependency locks](docs/development.md)
- [Worked example using public enterovirus A71 data](examples/eva71_sra/README.md)

The per-cell consensus step uses a modernized copy of
[sam2consensus](https://github.com/edgardomortiz/sam2consensus), included under
`workflow/scripts/`. Host transcriptome processing with Cell Ranger/Seurat is
outside this repository.

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
