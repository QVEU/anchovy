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

# THE SAMPLE NAME FOLLOWS THE READS, and has to be handed to BOTH halves of the
# run. fetch.sh names the SAM it writes after the FASTQ; the workflow reads
# `sample` from the config. Left to themselves those disagree the moment READS
# is not the default -- fetch.sh maps the new FASTQ, the workflow reads the
# PREVIOUS sample's SAM still sitting in data_dir, and reports it up to date.
# Nothing fails: the mapping is simply discarded and the results are the old
# sample's under the new sample's name.
#
# So derive it once here, pass it down to fetch.sh, and override the config key
# on the snakemake command line. Set SAMPLE to name it yourself.
. "$REPO/examples/eva71_sra/sample_name.sh"
SAMPLE="${SAMPLE:-$(sample_name_from_fastq "$READS")}"

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
printf '    sample  %s\n' "$SAMPLE"
printf '    config  %s\n' "$CONFIG"
printf '    threads %s\n\n' "$THREADS"

FASTQ="$READS" DATA_DIR="$DATA_DIR" THREADS="$THREADS" SAMPLE="$SAMPLE" \
    bash examples/eva71_sra/fetch.sh

printf '\n\033[1m==> Workflow\033[0m\n'
# --config OVERRIDES the configfile, so `sample` follows READS rather than
# whatever the config was last edited to say. `reads` is provenance only (the
# workflow never reads it), but a stale value there is a trap for anyone asking
# later what produced a results directory, so it is overridden too.
SM=(snakemake -s workflow/Snakefile --configfile "$CONFIG" --cores "$THREADS"
    --config sample="$SAMPLE" reads="$READS")

if [ -n "${DRY_RUN:-}" ]; then
    "${SM[@]}" -n
else
    "${SM[@]}" -n
    printf '\n    plan above; running it now\n\n'
    "${SM[@]}"
fi
