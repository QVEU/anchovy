#!/usr/bin/env bash
#
# fetch.sh -- build the SAM that anchovy starts from.
#
# Gets reads (local file or SRA download), a reference genome, and a region
# GFF3 derived from that reference's GenBank record, then maps the reads.
# Output is {DATA_DIR}/{SAMPLE}.sam, the input to workflow/Snakefile.
#
#   # your own reads:
#   FASTQ=/path/to/reads.fastq bash examples/eva71_sra/fetch.sh
#
#   # the public example run:
#   bash examples/eva71_sra/fetch.sh
#
# Steps whose output already exists are skipped, so a re-run resumes.

set -euo pipefail

# --------------------------------------------------------------------------- #
# Settings -- override from the environment.
# --------------------------------------------------------------------------- #

# Reads. Set FASTQ to use a local file and skip the SRA download entirely;
# .gz is fine. Leave it empty to download SRR instead.
FASTQ="${FASTQ:-}"
SRR="${SRR:-SRR28178313}"

# NCBI nucleotide accession. AF304458 = enterovirus A71 Tainan/4643/98.
REFERENCE_ACC="${REFERENCE_ACC:-AF304458}"

# minimap2 preset. map-hifi for PacBio (incl. .ccs), map-ont for Nanopore.
# The wrong preset does not fail; it silently costs alignments.
MINIMAP_PRESET="${MINIMAP_PRESET:-map-hifi}"

# 10X v2 barcode whitelist, downloaded if unset. Matches the 26-base signature
# in config.yaml (16 nt barcode + 10 nt UMI). v3 needs 28 Ns and the
# 3M-february-2018.txt whitelist -- both chemistries use 16 nt barcodes, so
# entry count, not barcode length, is what tells them apart.
WHITELIST="${WHITELIST:-}"
WHITELIST_URL="${WHITELIST_URL:-https://raw.githubusercontent.com/10XGenomics/supernova/refs/heads/master/tenkit/lib/python/tenkit/barcodes/737K-august-2016.txt}"
WHITELIST_EXPECTED_BARCODES="${WHITELIST_EXPECTED_BARCODES:-737280}"

# Output name. Defaults to the FASTQ's basename, or the accession.
if [ -n "$FASTQ" ]; then
    _base=$(basename "$FASTQ"); _base="${_base%.gz}"
    SAMPLE="${SAMPLE:-${_base%.fastq}}"; SAMPLE="${SAMPLE%.fq}"
else
    SAMPLE="${SAMPLE:-$SRR}"
fi

DATA_DIR="${DATA_DIR:-examples/eva71_sra/data}"
THREADS="${THREADS:-4}"

# Downloading only. MAX_SPOTS limits the transfer (needs fastq-dump);
# MAX_READS trims an already-downloaded FASTQ, saving mapping time only.
MAX_SPOTS="${MAX_SPOTS:-}"
MAX_READS="${MAX_READS:-}"

# --------------------------------------------------------------------------- #
say() { printf '\n\033[1m==> %s\033[0m\n' "$*"; }
die() { printf '\nerror: %s\n' "$*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "'$1' not found. $2"; }
uncat() { case "$1" in *.gz) zcat "$1";; *) cat "$1";; esac; }

say "Checking prerequisites"
need minimap2 "It is in environment.yml; did you 'conda activate anchovy'?"
need curl     "Required to download the reference."
need python   "Run this inside the anchovy conda environment."

if [ -n "$FASTQ" ]; then
    [ -f "$FASTQ" ] || die "FASTQ does not exist: $FASTQ"
    echo "    using local reads: $FASTQ"
else
    need fasterq-dump "Install sra-tools:
      conda env update -f environment.yml --prune
  or:
      conda install -c bioconda sra-tools"
fi

[ -z "$WHITELIST" ] || [ -f "$WHITELIST" ] \
    || die "WHITELIST was set but does not exist: $WHITELIST"

mkdir -p "$DATA_DIR"

REF_FA="$DATA_DIR/${REFERENCE_ACC}.fasta"
REF_GB="$DATA_DIR/${REFERENCE_ACC}.gb"
REF_GFF="$DATA_DIR/${REFERENCE_ACC}.gff3"
SAM="$DATA_DIR/${SAMPLE}.sam"
EFETCH="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${REFERENCE_ACC}"

# --------------------------------------------------------------------------- #
say "1/6  Barcode whitelist"
if [ -n "$WHITELIST" ]; then
    echo "    using your copy: $WHITELIST"
else
    WHITELIST="$DATA_DIR/737K-august-2016.txt"
    if [ -s "$WHITELIST" ]; then
        echo "    $WHITELIST exists, skipping."
    else
        echo "    downloading the 10X v2 whitelist (~12 MB)"
        curl -fsSL "$WHITELIST_URL" -o "$WHITELIST.part"
        mv "$WHITELIST.part" "$WHITELIST"
    fi
fi

# A truncated download or an error page would not error here -- it would just
# match fewer reads, which reads as bad data rather than a bad file.
WL_COUNT=$(grep -c . "$WHITELIST" || true)
WL_LENGTHS=$(awk 'NF {print length($1)}' "$WHITELIST" | sort -u | tr '\n' ' ')
echo "    $WL_COUNT barcodes, length(s): $WL_LENGTHS"

[ "$WL_LENGTHS" = "16 " ] \
    || die "whitelist barcodes are not a uniform 16 nt (got: $WL_LENGTHS).
  File: $WHITELIST"

if [ "$WL_COUNT" != "$WHITELIST_EXPECTED_BARCODES" ]; then
    echo "    note: expected $WHITELIST_EXPECTED_BARCODES (v2). ~3,000,000 means"
    echo "          the v3 whitelist, which needs a 28-N signature. Much smaller"
    echo "          usually means a truncated download."
fi

# --------------------------------------------------------------------------- #
say "2/6  Reference genome  ($REFERENCE_ACC)"
if [ -s "$REF_FA" ]; then
    echo "    $REF_FA exists, skipping."
else
    curl -fsSL "${EFETCH}&rettype=fasta&retmode=text" -o "$REF_FA.part"
    # A failed efetch can return 200 with an HTML or empty body.
    head -c1 "$REF_FA.part" | grep -q '>' \
        || die "downloaded reference is not FASTA -- is $REFERENCE_ACC valid?
  First line: $(head -1 "$REF_FA.part")"
    mv "$REF_FA.part" "$REF_FA"
    echo "    $(grep -c '^>' "$REF_FA") sequence(s), $(grep -v '^>' "$REF_FA" | tr -d '\n' | wc -c) nt"
fi

# Must match reference_name in the config: sam2consensus prefixes every
# per-cell output file with it.
REF_NAME=$(head -1 "$REF_FA" | sed 's/^>//' | awk '{print $1}')
echo "    reference_name: $REF_NAME"

# --------------------------------------------------------------------------- #
say "3/6  Region annotation  (GenBank -> GFF3)"
if [ -s "$REF_GFF" ]; then
    echo "    $REF_GFF exists, skipping."
else
    curl -fsSL "${EFETCH}&rettype=gb&retmode=text" -o "$REF_GB"
    python examples/eva71_sra/genbank_to_gff3.py "$REF_GB" "$REF_GFF"
fi

# --------------------------------------------------------------------------- #
if [ -n "$FASTQ" ]; then
    say "4/6  Reads  (local file)"
    echo "    $FASTQ"
else
    say "4/6  Reads  ($SRR)"
    FASTQ="$DATA_DIR/${SRR}.fastq"
    if [ -s "$FASTQ" ]; then
        echo "    $FASTQ exists, skipping."
    elif [ -n "$MAX_SPOTS" ]; then
        # fastq-dump -X stops the transfer; fasterq-dump has no supported
        # equivalent (--row-limit is per-thread and undocumented).
        command -v fastq-dump >/dev/null 2>&1 \
            || die "MAX_SPOTS needs fastq-dump (ships with sra-tools).
  Install it, or unset MAX_SPOTS for a full download."
        echo "    downloading the first $MAX_SPOTS spots only"
        fastq-dump -X "$MAX_SPOTS" --stdout "$SRR" > "$FASTQ.part" \
            || die "fastq-dump failed; unset MAX_SPOTS to fall back."
        mv "$FASTQ.part" "$FASTQ"
    else
        echo "    downloading the FULL run; this takes a while."
        echo "    (MAX_SPOTS=50000 pulls a small slice instead)"
        fasterq-dump "$SRR" --concatenate-reads --threads "$THREADS" \
                     --outdir "$DATA_DIR" --outfile "$(basename "$FASTQ")" --progress
    fi
fi

[ -s "$FASTQ" ] || die "no reads at $FASTQ."
FIRST=$(uncat "$FASTQ" | head -1) || true
case "$FIRST" in
    @*) ;;
    *)  die "not FASTQ -- first line: $FIRST" ;;
esac

FASTQ_READS=$(( $(uncat "$FASTQ" | wc -l) / 4 ))
echo "    $FASTQ_READS reads"

if [ -n "$MAX_SPOTS" ] && [ "$FASTQ_READS" -gt "$(( MAX_SPOTS * 2 ))" ]; then
    echo "    note: asked for $MAX_SPOTS spots, got $FASTQ_READS reads -- the"
    echo "          spot limit may not have been applied."
fi

# --------------------------------------------------------------------------- #
say "5/6  Optional subsample"
MAP_INPUT="$FASTQ"
if [ -n "$MAX_READS" ]; then
    SUB="$DATA_DIR/${SAMPLE}.subsample.fastq"
    if [ -s "$SUB" ]; then
        echo "    $SUB exists, skipping."
    else
        uncat "$FASTQ" | head -n "$(( MAX_READS * 4 ))" > "$SUB" || true
        echo "    took the first $MAX_READS reads into $SUB"
    fi
    MAP_INPUT="$SUB"
else
    echo "    MAX_READS unset -- using all reads."
fi

# --------------------------------------------------------------------------- #
say "6/6  Map to the reference  (minimap2 -ax $MINIMAP_PRESET)"
if [ -s "$SAM" ]; then
    echo "    $SAM exists, skipping."
else
    minimap2 -ax "$MINIMAP_PRESET" -t "$THREADS" "$REF_FA" "$MAP_INPUT" > "$SAM.part"
    mv "$SAM.part" "$SAM"
fi
MAPPED=$(awk '$1 !~ /^@/ && $2 != 4' "$SAM" | wc -l)
echo "    $MAPPED mapped alignment records"
[ "$MAPPED" -gt 0 ] || die "nothing mapped. Wrong reference, or wrong
  MINIMAP_PRESET for this read technology (currently $MINIMAP_PRESET)."

# --------------------------------------------------------------------------- #
cat <<EOF

Done. Set these in your config:

    sample: "$SAMPLE"
    data_dir: "$DATA_DIR"
    template: "$REF_FA"
    reference: "$REF_FA"
    reference_name: "$REF_NAME"
    whitelist: "$WHITELIST"
    gff: "$REF_GFF"

Then:

    snakemake -s workflow/Snakefile --configfile <your config> --cores $THREADS -n
EOF
