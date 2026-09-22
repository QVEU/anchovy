#!/usr/bin/env bash
#
# run_cluster.sh -- the QVEU EV-A71 6h P5 run, end to end, with the paths
# filled in. fetch.sh stays generic (FASTQ unset there means "download the
# public SRA example", so a lab path cannot be its default); this holds the
# ones specific to this dataset and this filesystem.
#
#   bash examples/eva71_sra/run_cluster.sh          # fetch, then run
#   DRY_RUN=1 bash examples/eva71_sra/run_cluster.sh   # fetch, then plan only
#
# Override any of the variables below from the environment.

set -euo pipefail

# The reads, where they live. READ IN PLACE -- fetch.sh maps straight from
# here and never copies, so there is no reason to stage an 11 GB FASTQ
# anywhere else.
READS="${READS:-/data/lvd_qve/Sequencing_Data/QVEU_Seq_0065_PacBio_NDAS_SearchSeq-EVA71Passage/Analysis_r54242Ue_20230322_183435_8_H01/Sample_5_EVA71_6h_P5/5_EVA71_6h_P5.ccs.fastq}"

# Where the reference, whitelist and mapped SAM go. OUTSIDE the checkout, so
# deleting or re-cloning the repo cannot take the mapping with it. Must match
# data_dir in the config.
DATA_DIR="${DATA_DIR:-/data/lvd_qve/Projects/PacBio_virus-inclusive_EVA71Passage/anchovy_run}"

CONFIG="${CONFIG:-workflow/config_cluster.yaml}"
THREADS="${THREADS:-64}"

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO"                    # results/ is relative; anchor it to the repo

[ -f "$READS" ] || {
    printf '\nerror: reads not found: %s\n' "$READS" >&2
    printf '  Set READS=/path/to/reads.fastq\n' >&2
    exit 1
}

[ -f "$CONFIG" ] || {
    printf '\nerror: config not found: %s\n' "$CONFIG" >&2
    printf '  Run this from the repo, or set CONFIG=/path/to/config.yaml\n' >&2
    exit 1
}

# CHECKED BEFORE FETCHING, not after. data_dir in the config has to agree with
# DATA_DIR or the workflow hunts for a SAM that fetch.sh just wrote somewhere
# else -- and finding that out after five minutes of mapping is a waste.
CONFIGURED=$(python -c "
import yaml
print(yaml.safe_load(open('$CONFIG'))['data_dir'])
" 2>/dev/null) || {
    printf '\nerror: could not read data_dir from %s\n' "$CONFIG" >&2
    printf '  Is it valid YAML, and is the anchovy environment active?\n' >&2
    exit 1
}
if [ "$CONFIGURED" != "$DATA_DIR" ]; then
    printf '\nerror: %s has\n      data_dir: %s\n' "$CONFIG" "$CONFIGURED" >&2
    printf '  but this run would fetch into\n      %s\n' "$DATA_DIR" >&2
    printf '  Point them at the same directory.\n' >&2
    exit 1
fi

printf '\n\033[1m==> Inputs\033[0m\n'
printf '    reads   %s\n' "$READS"
printf '    data    %s\n' "$DATA_DIR"
printf '    config  %s\n' "$CONFIG"
printf '    threads %s\n\n' "$THREADS"

FASTQ="$READS" DATA_DIR="$DATA_DIR" THREADS="$THREADS" \
    bash examples/eva71_sra/fetch.sh

printf '\n\033[1m==> Workflow\033[0m\n'
if [ -n "${DRY_RUN:-}" ]; then
    snakemake -s workflow/Snakefile --configfile "$CONFIG" --cores "$THREADS" -n
else
    snakemake -s workflow/Snakefile --configfile "$CONFIG" --cores "$THREADS" -n
    printf '\n    plan above; running it now\n\n'
    snakemake -s workflow/Snakefile --configfile "$CONFIG" --cores "$THREADS"
fi
