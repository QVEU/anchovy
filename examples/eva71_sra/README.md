# Worked example: EV-A71 single-cell data from the SRA

This runs anchovy end to end on real public data, starting from nothing but an
accession number. It downloads a single-cell sequencing run, fetches the
reference genome, works out the genome's regions from the reference's own
annotation, and hands the result to the normal workflow.

It's here to be copied. The two scripts are short and the settings file is
commented, so the usual way to run anchovy on your own virus is to copy this
directory and change the accessions.

## What it needs

Nothing but the anchovy conda environment. The reads, the reference genome and
the 10X barcode whitelist are all downloaded for you.

| Setting | Value | |
|---|---|---|
| Reads | `SRR28178313` | from the SRA |
| Reference | `AF304458` | EV-A71 Tainan/4643/98 |
| Read technology | `map-hifi` | PacBio |
| Chemistry | 26-base signature | 10X v2 |
| Whitelist | `737K-august-2016.txt` | 10X v2, 737,280 barcodes |

The whitelist comes from 10X's open-source
[supernova](https://github.com/10XGenomics/supernova) repository. It's the same
file Cell Ranger ships, so if you already have a copy you can use it instead:

```bash
WHITELIST=/path/to/737K-august-2016.txt bash examples/eva71_sra/fetch.sh
```

Either way `fetch.sh` checks the file really is a whitelist before going on —
barcodes must be a uniform 16 bases, and it tells you if the count isn't the
737,280 a v2 whitelist should have. A partial download would otherwise just
fail to match anything, which looks like bad data rather than a bad file.

Note that v2 and v3 whitelists *both* use 16-base barcodes, so barcode length
can't tell them apart — the entry count is what distinguishes them (~3,000,000
for v3).

### Checking the chemistry, if you're unsure

Using the wrong chemistry doesn't produce an error. The signature search still
runs, it just matches badly, and you end up with fewer cells than you should
have. You can see it in the extract stage's own output:

```bash
anchovy extract examples/eva71_sra/data/SRR28178313.sam \
    /path/to/737K-august-2016.txt -o /tmp/check.csv

python -c "
import pandas as pd
d = pd.read_csv('/tmp/check.csv')
print(f'reads assigned : {len(d)}')
print(f'distinct cells : {d.CBC.nunique()}')
print(d.minD.describe())
"
```

`minD` is how far each read's barcode was from its best match in the whitelist.
On the right chemistry most reads sit at a low distance and you get a sensible
number of distinct cells. On the wrong one the distances shift high and the cell
count collapses — try it with the other whitelist and compare if the numbers
look off.

## Running it

First, download everything and build the input file:

```bash
conda activate anchovy
bash examples/eva71_sra/fetch.sh
```

Then see what the workflow plans to do:

```bash
snakemake -s workflow/Snakefile --configfile examples/eva71_sra/config.yaml --cores 8 -n
```

and if that looks right, do it:

```bash
snakemake -s workflow/Snakefile --configfile examples/eva71_sra/config.yaml --cores 8
```

### Start small

The reads dominate step 1 — the whitelist is 12 MB and the reference is
trivial, but the sequencing run is a real one and downloading all of it takes a
long time. For a first pass, pull a slice instead:

```bash
MAX_SPOTS=50000 bash examples/eva71_sra/fetch.sh
```

`MAX_SPOTS` limits the **download**: only that many spots come down the wire.
That's the knob that makes this quick to iterate on. Fifty thousand is enough to
see cells, genotypes and annotations appear; drop to a few thousand if you just
want to watch the pipeline run end to end.

There's a second, different knob:

```bash
MAX_READS=200000 bash examples/eva71_sra/fetch.sh
```

`MAX_READS` trims a FASTQ you've **already downloaded** before mapping. It saves
mapping time only — by the time it applies, the whole run has come down. Use it
when you have the full data and want a faster mapping pass, not to shorten the
download.

When you want the real thing, delete the FASTQ and re-run with neither set.

`fetch.sh` prints the reference name it found when it finishes. It should match
`reference_name` in `config.yaml`; the version suffix on an accession can change
(`.1` vs `.2`), so it's worth a glance rather than an assumption.

Every step is skipped if its output is already there, so if something fails you
can fix it and re-run without starting over. Delete a file to redo that step.

## What `fetch.sh` does

1. **Downloads the 10X v2 barcode whitelist** and checks it really is one.
2. **Downloads the reference genome** as FASTA.
3. **Downloads its GenBank record and turns it into a region file**
   (`genbank_to_gff3.py`). This is the interesting part — see below.
4. **Downloads the sequencing reads** with `fasterq-dump`.
5. **Optionally takes a subsample**, if you set `MAX_READS`.
6. **Maps the reads to the reference** with minimap2.

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
