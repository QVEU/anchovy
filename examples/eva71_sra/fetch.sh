#!/usr/bin/env bash
#
# fetch.sh -- get the example's input data and build the SAM anchovy starts from.
#
# Downloads one public SRA run and one reference genome, derives a region GFF3
# from the reference's own GenBank annotation, and maps the reads. What comes
# out is {DATA_DIR}/{SAMPLE}.sam, which is exactly what the `extract` rule in
# workflow/Snakefile expects as its input -- so this script ends where the normal
# pipeline begins.
#
#   bash examples/eva71_sra/fetch.sh
#   snakemake -s workflow/Snakefile --configfile examples/eva71_sra/config.yaml --cores 8
#
# Every step is skipped if its output already exists, so re-running after a
# failure resumes rather than starting over. Delete a file to force that step.

set -euo pipefail

# --------------------------------------------------------------------------- #
# Settings -- override any of these from the environment, e.g.
#   THREADS=16 MAX_READS=100000 bash examples/eva71_sra/fetch.sh
# --------------------------------------------------------------------------- #

# The SRA run. Long-read single-cell cDNA; anchovy needs the 10X barcode still
# present in the read, so do NOT use a preprocessed/trimmed submission.
SRR="${SRR:-SRR28178313}"

# Reference genome, as an NCBI nucleotide accession.
# AF304458 = enterovirus A71 strain Tainan/4643/98. Confirmed by the lab.
REFERENCE_ACC="${REFERENCE_ACC:-AF304458}"

# minimap2 preset for the read technology. This run is PacBio, so map-hifi.
# Change to map-ont if you point this script at Oxford Nanopore data -- the
# wrong preset does not fail, it just quietly costs you alignments.
MINIMAP_PRESET="${MINIMAP_PRESET:-map-hifi}"

# 10X barcode whitelist. REQUIRED, and not downloadable here -- it ships with
# Cell Ranger. This run is v2 chemistry, so you want:
#
#   cellranger-x.y.z/lib/python/cellranger/barcodes/737K-august-2016.txt
#
# That matches the signature in config.yaml, whose 26-base N-run is 16 nt of
# barcode plus a 10 nt UMI -- the v2 layout. (v3 would be 28: a 12 nt UMI, and
# the 3M-february-2018.txt whitelist.) See "Checking the chemistry" in the
# README for how to confirm it from the data if you are unsure.
WHITELIST="${WHITELIST:-}"

SAMPLE="${SAMPLE:-$SRR}"
DATA_DIR="${DATA_DIR:-examples/eva71_sra/data}"
THREADS="${THREADS:-4}"

# Optional: cap the number of reads, to try the example without pulling the full
# run. Empty means all reads. e.g. MAX_READS=200000
MAX_READS="${MAX_READS:-}"

# --------------------------------------------------------------------------- #
say() { printf '\n\033[1m==> %s\033[0m\n' "$*"; }
die() { printf '\nerror: %s\n' "$*" >&2; exit 1; }

need() {
    command -v "$1" >/dev/null 2>&1 \
        || die "'$1' not found. $2"
}

say "Checking prerequisites"
need fasterq-dump "Install sra-tools: conda install -c bioconda sra-tools (it is in environment.yml)."
need minimap2     "Install minimap2: it is in environment.yml; did you 'conda activate anchovy'?"
need curl         "curl is required to download the reference."
need python       "Run this inside the anchovy conda environment."

[ -n "$WHITELIST" ] \
    || die "WHITELIST is not set. Point it at a 10X barcode whitelist:
    WHITELIST=/path/to/737K-august-2016.txt bash examples/eva71_sra/fetch.sh
  It ships with Cell Ranger -- see the comment at the top of this script."
[ -f "$WHITELIST" ] || die "whitelist not found: $WHITELIST"

mkdir -p "$DATA_DIR"

REF_FA="$DATA_DIR/${REFERENCE_ACC}.fasta"
REF_GB="$DATA_DIR/${REFERENCE_ACC}.gb"
REF_GFF="$DATA_DIR/${REFERENCE_ACC}.gff3"
FASTQ="$DATA_DIR/${SRR}.fastq"
SAM="$DATA_DIR/${SAMPLE}.sam"

EFETCH="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${REFERENCE_ACC}"

# --------------------------------------------------------------------------- #
say "1/5  Reference genome  ($REFERENCE_ACC)"
if [ -s "$REF_FA" ]; then
    echo "    $REF_FA exists, skipping."
else
    curl -fsSL "${EFETCH}&rettype=fasta&retmode=text" -o "$REF_FA.part"
    # A failed efetch can still return 200 with an HTML or empty body, which
    # would otherwise sail on and produce an empty index.
    head -c1 "$REF_FA.part" | grep -q '>' \
        || die "downloaded reference is not FASTA -- check that $REFERENCE_ACC is a valid accession.
  First line: $(head -1 "$REF_FA.part")"
    mv "$REF_FA.part" "$REF_FA"
    echo "    $(grep -c '^>' "$REF_FA") sequence(s), $(grep -v '^>' "$REF_FA" | tr -d '\n' | wc -c) nt"
fi

# The name after '>' has to match reference_name in config.yaml, because
# sam2consensus prefixes its per-cell output files with it.
REF_NAME=$(head -1 "$REF_FA" | sed 's/^>//' | awk '{print $1}')
echo "    reference_name: $REF_NAME"

# --------------------------------------------------------------------------- #
say "2/5  Region annotation  (GenBank -> GFF3)"
if [ -s "$REF_GFF" ]; then
    echo "    $REF_GFF exists, skipping."
else
    curl -fsSL "${EFETCH}&rettype=gb&retmode=text" -o "$REF_GB"
    # Derived from the record rather than hand-written, so the coordinates
    # cannot be mistyped -- and picornavirus records carry a mat_peptide per
    # cleavage product, which is what makes the long-format output interesting.
    python examples/eva71_sra/genbank_to_gff3.py "$REF_GB" "$REF_GFF"
fi

# --------------------------------------------------------------------------- #
say "3/5  Reads  ($SRR)"
if [ -s "$FASTQ" ]; then
    echo "    $FASTQ exists, skipping."
else
    echo "    This is a real sequencing run; expect a long download."
    fasterq-dump "$SRR" --concatenate-reads --threads "$THREADS" \
                 --outdir "$DATA_DIR" --outfile "$(basename "$FASTQ")" --progress
    [ -s "$FASTQ" ] || die "fasterq-dump produced no reads for $SRR."
fi
echo "    $(( $(wc -l < "$FASTQ") / 4 )) reads"

# --------------------------------------------------------------------------- #
say "4/5  Optional subsample"
MAP_INPUT="$FASTQ"
if [ -n "$MAX_READS" ]; then
    SUB="$DATA_DIR/${SRR}.subsample.fastq"
    if [ -s "$SUB" ]; then
        echo "    $SUB exists, skipping."
    else
        head -n "$(( MAX_READS * 4 ))" "$FASTQ" > "$SUB"
        echo "    took the first $MAX_READS reads into $SUB"
    fi
    MAP_INPUT="$SUB"
else
    echo "    MAX_READS unset -- using all reads."
fi

# --------------------------------------------------------------------------- #
say "5/5  Map to the reference  (minimap2 -ax $MINIMAP_PRESET)"
if [ -s "$SAM" ]; then
    echo "    $SAM exists, skipping."
else
    minimap2 -ax "$MINIMAP_PRESET" -t "$THREADS" "$REF_FA" "$MAP_INPUT" > "$SAM.part"
    mv "$SAM.part" "$SAM"
fi
MAPPED=$(awk '$1 !~ /^@/ && $2 != 4' "$SAM" | wc -l)
echo "    $MAPPED mapped alignment records"
[ "$MAPPED" -gt 0 ] || die "nothing mapped. Wrong reference, or the wrong
  MINIMAP_PRESET for this read technology (currently $MINIMAP_PRESET)."

# --------------------------------------------------------------------------- #
cat <<EOF

Done. Inputs are in $DATA_DIR:

    $(basename "$SAM")   the mapped reads anchovy starts from
    $(basename "$REF_FA")   reference genome
    $(basename "$REF_GFF")  region model derived from its GenBank annotation

Before running the workflow, set these in examples/eva71_sra/config.yaml:

    reference_name: "$REF_NAME"
    whitelist: "$WHITELIST"

then:

    snakemake -s workflow/Snakefile --configfile examples/eva71_sra/config.yaml --cores $THREADS -n
    snakemake -s workflow/Snakefile --configfile examples/eva71_sra/config.yaml --cores $THREADS
EOF
