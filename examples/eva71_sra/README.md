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

To run on reads you already have, point `FASTQ` at them and no download
happens. `.gz` is fine, and the sample name is derived from the filename:

```bash
FASTQ=/path/to/5_EVA71_6h_P5.ccs.fastq bash examples/eva71_sra/fetch.sh
```

`fetch.sh` prints the exact config lines to use when it finishes.

### A run with the paths already filled in

`workflow/config_cluster.yaml` holds the QVEU EV-A71 paths, so the whole thing
is one command:

```bash
snakemake -s workflow/Snakefile --configfile workflow/config_cluster.yaml --cores 64 -n
snakemake -s workflow/Snakefile --configfile workflow/config_cluster.yaml --cores 64
```

There is no run script. There used to be `run_cluster.sh`, which mapped the
reads and derived the sample name to hand to the workflow — both of which the
workflow now does itself from `input_dir`, so what was left was a single
`snakemake` line. Keeping it meant two places had to agree about how a run
works, which is exactly how the sample name came to be wrong when you pointed
it at different reads.

Run it from the repo root so `results/` lands there.

### Keep real-run data outside the repository

`DATA_DIR` defaults to `examples/eva71_sra/data`, **inside the checkout**. That
suits the example, whose inputs are all re-downloadable, and is a trap for a
real run: the mapped SAM lands there too, so an `rm -rf anchovy/` or a fresh
clone takes 11 GB of mapping with it. For anything you would rather not redo,
point it beside the repo and set the config to match:

```bash
DATA=/path/to/project/anchovy_run
FASTQ=/path/to/reads.fastq DATA_DIR="$DATA" THREADS=64 \
  bash examples/eva71_sra/fetch.sh
```

`workflow/config_cluster.yaml` is already written that way. Note `results/` is
still relative to wherever you run snakemake from, so run it from the repo root
or the outputs will follow you around.

### Use the whole machine

`snakemake --cores N` sets how many **jobs** run at once. It does not size the
pool inside `extract`, which is one job doing its own multiprocessing — that is
`extract_threads` in the config, and it sits at 16 unless you set it. On a
64-core node, leaving it unset costs roughly four hours on this dataset.

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

### Then go bigger

A subsample is for checking the pipeline runs, not for reading the biology off.
Spots are spread across ~750,000 possible barcodes, so cutting the download cuts
*per-cell depth*, and per-cell depth is what every downstream number rests on.
At 500,000 spots this run yields on the order of five cells past `depth_min: 3`,
each covered thinly enough that over half their genotype tokens come back as
ambiguity codes rather than called bases — see
[Ambiguity codes](#ambiguity-codes-and-why-they-arent-mutations) below.

Lowering `depth_min` does not fix that; it admits more cells at the same thin
coverage. The fix is more reads per barcode: raise `MAX_SPOTS` by an order of
magnitude, or drop it entirely.

When you want the real thing, delete the FASTQ and re-run with neither set.

`fetch.sh` prints the reference name it found when it finishes. It should match
the reference FASTA's header; the version suffix on an accession can change
(`.1` vs `.2`), so it's worth a glance rather than an assumption.

Every step is skipped if its output is already there, so if something fails you
can fix it and re-run without starting over. Delete a file to redo that step.

## What `fetch.sh` does

1. **Downloads the 10X v2 barcode whitelist** and checks it really is one.
2. **Downloads the reference genome** as FASTA.
3. **Downloads its GenBank record and turns it into a region file**
   (`genbank_to_gff3.py`). This is the interesting part — see below.
4. **Gets the sequencing reads** — your own file if you set `FASTQ`, otherwise
   downloaded with `fasterq-dump`.
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

**Not every record has them.** `AF304458` annotates only the polyprotein CDS,
which gives you one region and residue numbers against all 2,194 residues of it.
`fetch.sh` notices this and transfers the cleavage sites from an annotated
relative (`MATPEP_DONOR_ACC`, default `NC_001612`) by aligning the two
polyproteins — so the coordinates are derived, not typed, and a donor that is
not the same virus is refused rather than silently renumbering everything.

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

- `FASTQ` — your own reads, or `SRR` in `fetch.sh` for a run accession
- `REFERENCE_ACC` in `fetch.sh` — your reference
- `minimap_preset` in the config — `map-hifi` for PacBio, `map-ont` for Nanopore
  (it moved there when mapping became a pipeline stage)
- `input_dir`, `template`, `gff`, `chemistry` in `config.yaml` to match

`reference_name` is no longer a key: it is read from the reference FASTA's own
header, so there is nothing to keep in sync with it.

The region file is generated from whatever reference you name, so nothing there
needs editing by hand. If your reference has no annotation, write the GFF3
yourself — the main README's "Annotating by genome region" section explains the
format.

## What the variants are measured against

By default anchovy works out its own reference: the consensus of whatever cells
survived filtering. That answers "which cells differ from the crowd", which is
often what you want — but it needs a crowd. With a single surviving cell, the
reference *is* that cell, so its genotype comes out empty no matter what it
carries. An empty genotype then means "there was nothing to compare against",
which reads identically to "matches the virus".

This example sets `reference` to the genome instead:

```yaml
reference: "examples/eva71_sra/data/AF304458.fasta"
```

Each cell is then compared to EV-A71 itself. Two things change:

- **A single cell can carry variants.** Useful whenever filtering is strict or
  coverage is thin.
- **Mutations fixed across every cell are reported.** Against a computed
  reference they vanish, because they *become* the reference — the more
  completely a mutation has swept your population, the more certainly it
  disappears. On a passaged or lab-adapted stock that can be most of what you
  care about.

It has to be the same genome the reads were mapped to, since positions are
compared by index. anchovy checks the lengths match and stops if they don't.

Delete the line to go back to the cross-cell consensus.

## If you get fewer cells than you expected

Two filters drop cells, and both are set low in this example so it produces
something to look at:

- `cons_min_depth` — coverage a cell needs *at a position* for sam2consensus to
  call it. A cell that never reaches it produces no consensus at all.
- `depth_min` — coverage a cell needs overall to enter the consensus stage.

On a real run most barcodes carry few reads, so the package defaults (5 and 10)
can leave very few cells standing on a subsample. Raise them for a deeply
sequenced run, where you can afford to demand more evidence per cell.

anchovy warns if fewer than two cells survive with no reference supplied, since
that combination cannot produce a genotype.

## Ambiguity codes, and why they aren't mutations

sam2consensus does not only emit A/C/G/T. Where a position's reads disagree it
emits an IUPAC ambiguity code — `R` for A-or-G, `Y` for C-or-T, `S` for G-or-C —
and the lowercase form when a gap or an `N` was among the observed bases. Those
mark positions where the reads **disagreed**, not positions where the cell
carries a mutation.

anchovy does not turn them into genotype tokens. A genotype token is a claim
about a base, and the token machinery compares cells by token identity — so two
cells both reading `R` at 3185 would be grouped as *sharing a mutation*, when
what they actually share is not having enough reads to call one. The consensus
stage prints how many calls it dropped, so the filtering is never silent.

On this example's data at `depth_min: 3` it mattered a lot. Five cells gave
eleven distinct tokens, six of which were ambiguity codes:

```
1294C_2873C_3176S_3185R_3587R_3751C_6058R   →   1294C_2873C_3751C
1294C_3751C                                 →   1294C_3751C
2873C                                       →   2873C
3288T_5059C_5431Y                           →   3288T_5059C
3288T_6202R                                 →   3288T
```

The filtered column is not just tidier, it is more informative. Single-step
edges in the genotype network went from **0 to 2**, because they are defined by
one genotype's mutation set being a subset of another's, and the uncalled
positions were padding those sets apart. `1294C_3751C` sits inside
`1294C_2873C_3751C`; with four phantom mutations in the way, nothing sat inside
anything.

The real fix is coverage. These are self-inflicted at 3–5× depth, where a single
discordant read stops any base reaching the consensus threshold. **Raising
`cons_threshold` makes this worse, not better** — the algorithm accumulates
bases most-frequent-first until their *combined* coverage reaches the threshold,
so a higher threshold pulls in more bases and emits more ambiguity. 3A/2G calls
`A` at 0.5 and `R` at 0.75. More reads per cell is the only thing that resolves
them.

Set `keep_ambiguous: true` if you want them kept — reasonable if you are
deliberately after genuine within-cell mixed populations, which for a viral
quasispecies can be real signal rather than noise.

## Setting an analysis window

Real runs have ragged coverage at the ends of the genome, where only a few cells
have any reads. Those positions produce noisy variant calls.

`orf_start` and `orf_end` in `config.yaml` restrict which positions get called
without renumbering anything. They're left unset here on purpose: the right
window depends on where *this* run actually has coverage. Run it once, look at
the per-cell consensus sequences, then set the window to the well-covered core
and run again.
