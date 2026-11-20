# On-Chip Multiplexing (OCM)

[Back to the README](../README.md)

This page covers running a 10x **On-Chip Multiplexing** library through
anchovy. With OCM, up to four samples share one 10x channel and so one
long-read FASTQ. anchovy splits that FASTQ into its samples itself, so each
sample gets its own cells, consensus sequences, genotypes, networks and
frequencies.

If your samples were pooled some other way (hashtag antibodies, CellPlex), or
each sample has its own channel, this page does not apply. See
[whitelists and cell lists](../README.md#whitelists-and-cell-lists)
in the README instead.

## Is my run OCM?

Open the `cellranger multi` config CSV for the channel. An OCM run has an
`ocm_barcode_ids` column in its `[samples]` section:

```
[samples]
sample_id,ocm_barcode_ids
P10_5,OB1
P10_6,OB2
P10_7,OB3
P10_8,OB4
```

OCM kits are the GEM-X Universal 3′ v4 and 5′ v3 "4-plex" assays. One channel
has four inlets, and each inlet gets its own set of gel beads, labelled OB1 to
OB4.

## How the samples are told apart

Every read's 16 nt cell barcode carries its inlet. **Bases 8–9 of the barcode
are the OCM "overhang"**. Cell Ranger's chemistry definitions put it there
for every OCM chemistry (3′ v3/v4, 5′ v3, V(D)J): an `overhang` barcode at
offset 7, length 2, inside the gel-bead barcode, checked against one shared
`overhang` list
([chemistry_defs.json](https://github.com/10XGenomics/cellranger/blob/main/lib/rust/cr_types/src/chemistry/chemistry_defs.json)).
Each OB set uses its own pair of bases:

| Inlet | Barcode bases 8–9 |
|---|---|
| OB1 | GT |
| OB2 | CA |
| OB3 | TC |
| OB4 | AG |

This table was read from a real GEM-X 5′ v3 4-plex run, where every called
cell of each sample carried its inlet's pair and no other. It is built into
anchovy (`anchovy.ocm.OB_OVERHANGS`). 10x does not publish the table, so if
a run ever disagrees, `anchovy ocm-table` (Step 1) will say so.

Two consequences:

- **The barcode alone identifies the sample.** No tag reads and no Cell Ranger
  output are needed at run time.
- **Which bases mean OB1, OB2, OB3 and OB4 is the same in every OCM run.**
  The only thing that changes from run to run is which sample went into
  which inlet, and that is the multi config's `ocm_barcode_ids`.

The full 10x 5′ v3 whitelist (`3M-5pgex-jan-2023.txt`) holds every OB set, so
leave `whitelist` unset and let `chemistry: "5p-v3"` pick it.

## Step 1: get the sample → OB lines

All a run adds is which sample went into which inlet. Copy that from the
multi config's `[samples]` section into `ocm:` (Step 2). That is all most runs
need.

To check a run against the built-in table, or to write the `ocm:` block for
you, `anchovy ocm-table` reads the channel's `cellranger multi` output. The
folders under `per_sample_outs` are named by `sample_id`, and each holds only
that sample's called cells. Give it the multi config too:

```bash
anchovy ocm-table cellranger/p10_5-8/outs/per_sample_outs \
    --multi-config p10_5-8_multi_config.csv
```

It prints:

- the share of each sample's cells that carry each pair of bases, as
  evidence;
- whether this run matches the built-in OB table;
- an `ocm_overhangs:` block (OB → bases, as seen in this run);
- an `ocm:` block, which is sample → OB, ready to paste.

```yaml
# Overhang (barcode bases 8-9) share of each sample's called cells:
#   P10_5 (1,776 cells): GT 100.0%
#   ...
#
# Built-in OB table: MATCHES this run.
ocm_overhangs:
  OB1: [GT]
  OB2: [CA]
  OB3: [TC]
  OB4: [AG]
ocm:
  <channel FASTQ name>:
    P10_5: OB1
    P10_6: OB2
    P10_7: OB3
    P10_8: OB4

cell_lists:                          # optional: Cell Ranger / Seurat cells, per sample
  P10_5: /data/cr/p10_5-8/outs/per_sample_outs/P10_5/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz
```

You don't need to paste `ocm_overhangs:`. The built-in table covers it. If
`ocm-table` ever says the built-in table **differs**, please report it, and
paste its `ocm_overhangs:` block into the settings file in the meantime.

`ocm-table` refuses to print a table if a pair of bases carries at least 1%
of two samples' cells. That would mean the barcode position is wrong for this
chemistry, or the run was not OCM. Without `--multi-config` it prints each
sample's bases directly instead of OB ids.

## Step 2: the settings file

```yaml
input_dir: "/data/p10/pacbio"        # holds p10_5-8.fastq.gz, and any other FASTQs
template: "/data/p10/reference.fasta"
reference: "/data/p10/reference.fasta"
gff: "/data/p10/regions.gff3"
chemistry: "5p-v3"                   # leave `whitelist` unset: the full list holds every OB set
# max_barcode_errors: 1              # optional; the default keeps exact barcode matches only

ocm:
  p10_5-8:                           # the channel: a FASTQ in input_dir
    P10_5: OB1                       # sample name: its inlet
    P10_6: OB2
    P10_7: OB3
    P10_8: OB4
```

### The two kinds of name

| Where | What it is | Rules |
|---|---|---|
| Key under `ocm:` (`p10_5-8`) | The **channel FASTQ**: its filename without `.fastq`, `.fq` or `.gz` | Must match a FASTQ in `input_dir` |
| Name on each OB (`P10_5`) | A **sample** | Any name. It becomes the results folder and file prefix. It must be unique in the run (not another FASTQ's or OCM sample's name) and usable as a folder name (no `/`; avoid spaces) |

The sample names do not have to match Cell Ranger's, but using the
`sample_id`s keeps the later join with Seurat simple. `ocm-table` fills them
in for that reason.

### Other forms

- **A sample loaded into two inlets:** `P10_5: OB1|OB2`, as in the multi
  config.
- **Bases instead of an OB id:** `P10_5: [GT]`.
- **A different OB table:** `ocm_overhangs: {OB1: [..], ...}` replaces the
  built-in one.
- **An unused inlet:** leave its OB out. Its reads go to no sample.
- **Several channels:** one entry each under `ocm:`. FASTQs in `input_dir`
  that are not named under `ocm:` run as ordinary samples.

## What runs

```
p10_5-8.fastq.gz ─ map_reads ─ extract ─► results/p10_5-8/p10_5-8_anchovy.csv        once per channel
                                             │
                                ocm_split: keep reads whose barcode bases 8-9
                                belong to the sample's OB
            ┌───────────────┬───────────────┼───────────────┐
   P10_5_anchovy.csv  P10_6_anchovy.csv  P10_7_anchovy.csv  P10_8_anchovy.csv
            │  cells → consensus → genotypes → networks → frequencies → report
   results/P10_5/      results/P10_6/      results/P10_7/      results/P10_8/
```

- Mapping and extract, the expensive stages, run **once per channel**.
- Each read goes to **at most one** sample. Reads whose bases belong to no
  listed sample are dropped. They come from an unused inlet, from background
  barcodes, or from a sequencing error in those two bases.
- Everything after the split is an ordinary per-sample run. The channel
  itself has no results beyond its SAM, its anchovy CSV, and
  `<channel>_barcodes.csv`: every matched barcode in the channel, with an
  `ocm_sample` column naming the sample its bases belong to.
- Each sample gets `results/<sample>/<sample>_barcodes.csv`, and one
  `results/<sample>/cell_lists/<sample>_<label>.csv` per entry under
  `cell_lists`. A sample with a cell list runs consensus and the networks on
  those cells only, in `results/<sample>_<label>/`, instead of on every
  barcode. A later Cell Ranger or Seurat list is matched with
  `anchovy match-cells results/<sample>/<sample>_barcodes.csv <list>`.

Each sample's split prints the share of the channel's reads it kept:

```
P10_5: kept 412,330 of 1,604,221 reads (25.7%) whose barcode carries overhang GT.
```

With four samples loaded evenly, expect roughly a quarter each, a little
less after background. A sample far below that, or at zero, usually has the
wrong OB in the config.

## Checks before anything runs

The workflow stops with a message, before any job, when:

- an `ocm:` key does not match a FASTQ in `input_dir`;
- an OB id is unknown;
- a pair of bases is not exactly two of A/C/G/T;
- the same bases, or the same OB, are given to two samples;
- a sample name is already used by another FASTQ or OCM sample;
- `ocm:` is combined with `samples` / `cells_dir`.

## Running the split by hand

The split is one command per sample, on the channel's anchovy CSV:

```bash
anchovy ocm-split results/p10_5-8/p10_5-8_anchovy.csv P10_5_anchovy.csv \
    --overhangs GT --sample P10_5
anchovy fasta P10_5_anchovy.csv P10_5_cells/
```

It reads the channel CSV in chunks, so memory stays small on large runs. It
always writes a header, so a sample with no reads still gives the next stage
a readable, empty table.

## Troubleshooting

**`ocm-table` says a pair of bases is in two samples.** The run is probably
not OCM. Check the multi config for `ocm_barcode_ids`. For hashtag or
CellPlex pooling, split with the barcode-to-sample join in the README.

**A sample keeps far fewer reads than expected.** Compare its OB in `ocm:`
with the multi config's `[samples]`. Then run `ocm-table` and check that
sample's bases.

**`extract` stops on a low exact-match rate.** That check runs on the whole
channel before the split. With the full 5′ v3 list it should pass. If it
fails, check `chemistry` (5′ and 3′ kits use different lists) and that
`whitelist` is unset or points at the full list. Give each sample's Cell
Ranger cells as `cell_lists`, never as `whitelist`.
