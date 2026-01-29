# anchovy

![anchovies](https://github.com/QVEU/anchovy/blob/main/assets/northern-anchovies-rw07-130.webp)

anchovy is a small pipeline to generate consensus viral sequences for individual cells in single-cell experiments. It was developed for long-read cDNA data (Oxford Nanopore or PacBio) and is intended to be run on mapped read files (SAM/BAM) produced by minimap2 against a viral reference.

## Features
- Extracts 10X cell barcodes and UMI signatures from mapped reads
- Groups reads by cell barcode and collapses UMIs to per-cell FASTAs
- Re-maps per-cell FASTAs and generates a consensus sequence per cell
- Filters consensus sequences by genomic region and prepares outputs for downstream annotation and network analysis

## Repository contents
- `anchovy.py` — main script that processes mapped SAM files and extracts reads per cell
- `UMItoFasta.py` — converts `anchovy.csv` into per-cell FASTA files (one FASTA per CBC)
- `sam2consensus.py` — small script to generate a consensus for CBC-specific SAM files (requires external dependency)
- `ConsensusTool.py` — filters and compiles consensus sequences across cells for a target ORF/region
- `Consensus_Annotation.R` — accessory R script to translate and annotate consensus sequences and generate genotype networks
- `AnchovyJob_script.sh` — example pipeline wrapper scripts (PB and NP variants included)

## Prerequisites
To run the core Anchovy pipeline (steps 1–5):
- Python 3.8+
- minimap2 (for mapping FASTA -> SAM)
- sam2consensus.py (from https://github.com/edgardomortiz/sam2consensus) — used for consensus generation

For downstream annotation:
- Python for `ConsensusTool.py` (notebook / scripts; see file header for specifics)
- R >= 3.6 (recommended R 4.x) with ggplot2 for `Consensus_Annotation.R`
- RStudio (optional)

## Installation (local)
Clone the repository and change into the project directory:

```bash
# Option A: GitHub CLI
gh repo clone QVEU/anchovy

# Option B: HTTPS
git clone https://github.com/QVEU/anchovy.git

cd anchovy
```

The scripts are also available on the QVEU shared space at:
`/data/lvd_qve/QVEU_code/Anchovy/` (internal)

## Quickstart / Example workflow

Typical flow (high level):
1. Map raw reads to a viral reference with minimap2 to produce a SAM file.
2. Run `anchovy.py` on the SAM file to extract the 10X barcode signature and produce `anchovy.csv`.
3. Convert `anchovy.csv` to per-cell FASTAs using `UMItoFasta.py`.
4. Re-map per-cell FASTAs to the reference (minimap2) to produce per-cell SAM files.
5. Generate consensus sequences per cell using `Sam2Consensus.py` / `sam2consensus.py`.
6. Concatenate per-cell consensus sequences into an `allConsensus.fasta`.
7. Filter and process consensus sequences with `ConsensusTool.py`.
8. Annotate and analyze with `Consensus_Annotation.R`.

### Example: Interactive Anchovy run
Edit `AnchovyJob_script.sh` for your paths (examples for PacBio and Nanopore are included), then run interactively:

```bash
sh AnchovyJob_script.sh /path/to/whitelist_file /path/to/input_dir /path/to/template.fasta inputfilename.sam
```

or submit as a Slurm job:

```bash
sbatch AnchovyJob_script.sh /path/to/whitelist_file /path/to/input_dir /path/to/template.fasta inputfilename.sam
```

Anchovy pipeline outputs (typical):
- `anchovy.csv` — read-to-barcode/UMI mapping
- per-CBC FASTA files (one file per cell barcode)
- per-CBC SAM files (reads re-mapped to reference)
- per-CBC consensus FASTA sequences
- `allConsensus.fasta` — concatenated per-cell consensus sequences

## Consensus filtering — ConsensusTool.py

ConsensusTool.py extracts and filters consensus sequences that fully cover a genomic region (e.g., an ORF) and prepares outputs for annotation.

Usage:
```bash
python ConsensusTool.py <NAME>_allConsensus.fasta <ORF_start_nt> <ORF_end_nt>
```

Example:
```bash
python ConsensusTool.py /Volumes/LVD_qve/Projects/DENV_SEARCHLIGHT/Anchovy_2/DENV_6dpi_allConsensus.fasta 96 10272
```

Outputs:
- `<NAME>_consensus_reference.txt` — reference consensus used for annotation
- `<NAME>_filtConsensus.csv` — table of filtered consensus sequences and metadata
- `<NAME>_filtConsensus.fasta` — filtered consensus sequences in FASTA format

## Post-analysis — Consensus_Annotation.R

This R script translates consensus sequences, annotates mutations, and generates genotype/epistatic networks.

Usage:
```bash
Rscript /path/to/anchovy/Consensus_Annotation_v3.R /path/to/<filtConsensus_reference.txt> /path/to/<filtConsensus.csv> /path/to/output_prefix
```

Example:
```bash
Rscript ~/Documents/GitHub/anchovy/Consensus_Annotation.R ~/reference_consensus.txt filtConsensus.csv /path/to/output/prefix
```

Typical outputs:
- `<prefix>.pdf` — plots and visual summaries
- `<prefix>_annot_v3.csv` — annotated consensus table
- `<prefix>_epistaticNetwork.csv`
- `<prefix>_genotypeNetwork.csv`

## Notes and tips
- Barcode whitelist lists: 10X barcode inclusion lists (whitelists) ship with the CellRanger installation; e.g.:
  `/path/to/cellranger-version/lib/python/cellranger/barcodes/`
  Alternatively, public lists (3M, 3M-3pgex) can be downloaded from sources like the Teichlab site.
- For v4 chemistry (16-nt barcodes), use the appropriate whitelist (e.g., `3M-3pgex-may-2023.txt.gz`) when extracting CBCs.
- Adjust mapping and consensus parameters depending on read error profiles (Nanopore vs PacBio).

## Citation
If you use anchovy in published work, please cite:
N. Dábilla, P. T. Dolan, Structure and dynamics of enterovirus genotype networks. Sci Adv 10, eado1693 (2024).

## Contributing
- Fork the repository
- Create a feature branch: `git checkout -b feature-branch`
- Commit your changes: `git commit -am "Add new feature"`
- Push to your branch: `git push origin feature-branch`
- Open a pull request

## License & Contact
- License: MIT License, Copyright (c) 2026 Patrick Dolan
- Contact / Maintainers: @ptdolan Patrick Dolan, NIAID/NIH
