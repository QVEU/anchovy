# Configuration and detailed usage

[Back to the README](../README.md)

Run commands from the repository root. Relative paths in configuration are resolved from the working directory.

## Configuration reference

`workflow/config_cluster.yaml` is the fully commented reference — it sets most
of these and says why. `workflow/config_5p_V3.yaml` is the opposite: a
minimally commented 5' v3 cluster config with five paths to fill in and
working defaults for the rest. The complete list:

| Setting | Default | What it does |
|---|---|---|
| **Input and output** | | |
| `input_dir` | — | The folder of FASTQs. Required, unless resuming with `samples` + `cells_dir` |
| `samples` | from the FASTQ names | Name the samples explicitly. Only for the `cells_dir` resume path, where there are no FASTQs to take a name from |
| `cells_dir` | `<results>/cells` | Resume from per-cell FASTAs you already have |
| `results_dir` | `<input_dir>/results` | Move the outputs |
| `resources_dir` | `<input_dir>/resources` | Move the whitelist cache |
| **Reference** | | |
| `template` | — | Reference FASTA. Required. Must closely match the sequenced virus; see [reference choice](../README.md#use-your-own-data) |
| `reference` | unset | Call genotypes against the genome rather than the consensus across cells. Set it to the same file as `template` |
| `reference_name` | from the FASTA header | Override the contig name |
| `gff` | unset | Region model (see below). Implies `whole_reference` |
| `whole_reference` | `false` | Keep genome coordinates without a region model |
| `minimap_preset` | `map-hifi` | `map-ont` for Nanopore |
| `map_threads` | `--cores` | Threads for the initial mapping. Follows the run's core count unless set |
| **Barcodes** | | |
| `chemistry` | `v3` | `v2`, `v3` or `5p-v3`. Sets the read signature *and* the whitelist together |
| `whitelist` | downloaded | A local or run-specific barcode list |
| `whitelist_url` | 10X's | Fetch the whitelist from your own mirror |
| `whitelist_barcodes` | per chemistry | Override the barcode count a download is checked against |
| `signature` | from `chemistry` | For an assay whose handles differ from 10X's |
| `max_barcode_errors` | unset | Errors a barcode may carry. Unset, every read is assigned to its nearest entry with no floor. **The stage's biggest runtime lever** — see [performance tuning](performance.md) |
| `max_distance` | `42` | How far the *signature* match may be before a read is dropped |
| `extract_threads` | `--cores` | Worker pool inside `extract`. Follows the run's core count unless set. **Sets the stage's memory** — see [performance tuning](performance.md) |
| `downstream_window` | `12` | How far past the alignment offset to search for the signature — see [performance tuning](performance.md) |
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
snakemake -s workflow/Snakefile --configfile my_settings.yaml --cores 8 --keep-going
```

`-n` is a dry run: it lists the jobs and the samples it found without doing any
work, which is the quickest way to catch a wrong `input_dir`. On a cluster node
with more cores, raise `--cores` to match:

```bash
snakemake -s workflow/Snakefile --configfile workflow/config_cluster.yaml --cores 64 -n
snakemake -s workflow/Snakefile --configfile workflow/config_cluster.yaml --cores 64 --keep-going
```

**Use `--keep-going` on multi-sample runs.** Without it, snakemake stops the
whole run at the first failed job, so one sample that yields no usable cells
(a near-empty library, or every cell dropped by `min_breadth`) also halts every
other sample partway. With it, the healthy samples run to their reports and
the failed one is listed at the end; fix or drop it from the config and rerun,
and only its jobs are redone.

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
pointed at your own virus. See the [worked-example README](../examples/eva71_sra/README.md).

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
