# Worked example: EV-A71 single-cell data from the SRA

This runs anchovy end to end on real public data, starting from nothing but an
accession number. It downloads a single-cell sequencing run, fetches the
reference genome, works out the genome's regions from the reference's own
annotation, and hands the result to the normal workflow.

It's here to be copied. The two scripts are short and the settings file is
commented, so the usual way to run anchovy on your own virus is to copy this
directory and change the accessions.

## Before you start

Two things you have to supply or check yourself.

**A 10X barcode whitelist.** anchovy matches each read's cell barcode against a
list of the barcodes that exist in the chemistry you used. That list ships with
Cell Ranger and isn't ours to redistribute, so you point at your copy:

- v2 chemistry: `cellranger-x.y.z/lib/python/cellranger/barcodes/737K-august-2016.txt`
- v3 chemistry: the same folder's `3M-february-2018.txt.gz` — unzip it first

The signature in `config.yaml` expects 16 barcode bases plus a 10-base UMI,
which is the v2 layout. If the run used v3, use the v3 whitelist.

> **Two settings in these files could not be verified when they were written,
> because the machine had no access to NCBI. Please check both before you trust
> any results.**
>
> - **`REFERENCE_ACC=AF304458`** in `fetch.sh`, meant to be enterovirus A71
>   strain Tainan/4643/98. If it's the wrong accession you won't get an error —
>   you'll get a full set of confident variant calls against the wrong genome.
>   Check it at <https://www.ncbi.nlm.nih.gov/nuccore/AF304458>.
> - **`MINIMAP_PRESET=map-hifi`**, which assumes PacBio HiFi reads. For Oxford
>   Nanopore use `map-ont`. The wrong one quietly costs you alignments rather
>   than failing.

## Running it

```bash
conda activate anchovy

# 1. Download the data and build the input file.
WHITELIST=/path/to/737K-august-2016.txt bash examples/eva71_sra/fetch.sh

# 2. Put the whitelist path and the reference name it prints into config.yaml.

# 3. See what the workflow plans to do, then do it.
snakemake -s workflow/Snakefile --configfile examples/eva71_sra/config.yaml --cores 8 -n
snakemake -s workflow/Snakefile --configfile examples/eva71_sra/config.yaml --cores 8
```

This is a real sequencing run, so step 1 takes a while. To try the example on a
slice of it first:

```bash
MAX_READS=200000 WHITELIST=/path/to/whitelist.txt bash examples/eva71_sra/fetch.sh
```

Every step is skipped if its output is already there, so if something fails you
can fix it and re-run without starting over. Delete a file to redo that step.

## What `fetch.sh` does

1. **Downloads the reference genome** as FASTA.
2. **Downloads its GenBank record and turns it into a region file**
   (`genbank_to_gff3.py`). This is the interesting part — see below.
3. **Downloads the sequencing reads** with `fasterq-dump`.
4. **Optionally takes a subsample**, if you set `MAX_READS`.
5. **Maps the reads to the reference** with minimap2.

What comes out is `SRR28178313.sam`, which is exactly what the pipeline's first
stage expects. From there it's the ordinary workflow.

## Why the region file is generated, not written by hand

anchovy can annotate mutations by genome region, which needs a GFF3 file saying
where each region is. The obvious thing would be to ship one for EV-A71 with the
coordinates typed in.

We don't, because a wrong coordinate in that file doesn't cause an error. It
silently renumbers every amino acid after it, and the results look completely
normal. Deriving the coordinates from the GenBank record the pipeline already
downloads removes that risk, and means the example still works when you point it
at a different virus.

You also get the mature peptides for free. Picornavirus records annotate the
polyprotein *and* each protein cut out of it — VP1 to VP4, 2A to 2C, 3A to 3D —
so a mutation in VP1 comes out described twice:

```
mutation_id            region       region_type  residue
polyprotein:C2891T     polyprotein  coding           731
VP1:C2891T             VP1           coding           156
```

Residue 156 of VP1 is residue 731 of the polyprotein. Both are correct; which
one you want depends on what you're comparing against. Papers usually number
within the mature protein, so having both saves converting by hand.

## Adapting it to your own virus

Copy this directory and change:

- `SRR` in `fetch.sh` — your run accession
- `REFERENCE_ACC` in `fetch.sh` — your reference
- `MINIMAP_PRESET` — `map-hifi` for PacBio, `map-ont` for Nanopore
- `sample`, `template`, `reference_name`, `gff` in `config.yaml` to match

The region file is generated from whatever reference you name, so nothing there
needs editing by hand. If your reference has no annotation, write the GFF3
yourself — the main README's "Annotating by genome region" section explains the
format.

## Setting an analysis window

Real runs have ragged coverage at the ends of the genome, where only a few cells
have any reads. Those positions produce noisy variant calls.

`orf_start` and `orf_end` in `config.yaml` restrict which positions get called
without renumbering anything. They're left unset here on purpose: the right
window depends on where *this* run actually has coverage. Run it once, look at
the per-cell consensus sequences, then set the window to the well-covered core
and run again.
