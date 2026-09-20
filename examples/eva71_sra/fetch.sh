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

# 10X barcode whitelist. Downloaded automatically -- this run is v2 chemistry,
# and 10X publish that whitelist in the open-source supernova repository, so the
# example needs no manual setup. Set WHITELIST to use a local copy instead (the
# identical file ships with Cell Ranger as
# cellranger-x.y.z/lib/python/cellranger/barcodes/737K-august-2016.txt).
#
# It matches the signature in config.yaml, whose 26-base N-run is 16 nt of
# barcode plus a 10 nt UMI -- the v2 layout. v3 would be 28 Ns (12 nt UMI) and
# the 3M-february-2018.txt whitelist; note both chemistries use 16 nt barcodes,
# so the ENTRY COUNT is what distinguishes them, not the barcode length.
WHITELIST="${WHITELIST:-}"
WHITELIST_URL="${WHITELIST_URL:-https://raw.githubusercontent.com/10XGenomics/supernova/refs/heads/master/tenkit/lib/python/tenkit/barcodes/737K-august-2016.txt}"
# Expected size of the v2 whitelist, used to verify the download completed.
WHITELIST_EXPECTED_BARCODES="${WHITELIST_EXPECTED_BARCODES:-737280}"

SAMPLE="${SAMPLE:-$SRR}"
DATA_DIR="${DATA_DIR:-examples/eva71_sra/data}"
THREADS="${THREADS:-4}"

# TWO WAYS TO WORK ON LESS DATA. They are not interchangeable:
#
#   MAX_SPOTS  limits the DOWNLOAD. Only this many spots are transferred, so
#              this is the one that makes the example quick to iterate on.
#              Needs fastq-dump (ships with sra-tools alongside fasterq-dump).
#
#   MAX_READS  trims an ALREADY-DOWNLOADED FASTQ before mapping. It saves
#              mapping time only -- the whole run has already come down the
#              wire by then, so it does nothing for the slow part.
#
# Start with MAX_SPOTS. e.g. MAX_SPOTS=50000 bash examples/eva71_sra/fetch.sh
MAX_SPOTS="${MAX_SPOTS:-}"
MAX_READS="${MAX_READS:-}"

# --------------------------------------------------------------------------- #
say() { printf '\n\033[1m==> %s\033[0m\n' "$*"; }
die() { printf '\nerror: %s\n' "$*" >&2; exit 1; }

need() {
    command -v "$1" >/dev/null 2>&1 \
        || die "'$1' not found. $2"
}

say "Checking prerequisites"
need fasterq-dump "Install sra-tools. It is in environment.yml, but an
  environment created before it was added will not have it -- update yours:
      conda env update -f environment.yml --prune
  or install just this one tool:
      conda install -c bioconda sra-tools"
need minimap2     "Install minimap2: it is in environment.yml; did you 'conda activate anchovy'?"
need curl         "curl is required to download the reference."
need python       "Run this inside the anchovy conda environment."

[ -z "$WHITELIST" ] || [ -f "$WHITELIST" ] \
    || die "WHITELIST was set but does not exist: $WHITELIST
  Unset it to download the v2 whitelist automatically."

mkdir -p "$DATA_DIR"

REF_FA="$DATA_DIR/${REFERENCE_ACC}.fasta"
REF_GB="$DATA_DIR/${REFERENCE_ACC}.gb"
REF_GFF="$DATA_DIR/${REFERENCE_ACC}.gff3"
FASTQ="$DATA_DIR/${SRR}.fastq"
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

# Verify it really is a whitelist. A partial download or an error page would
# otherwise sail through and simply fail to match anything, which looks like bad
# data rather than a bad file.
WL_COUNT=$(grep -c . "$WHITELIST" || true)
WL_LENGTHS=$(awk 'NF {print length($1)}' "$WHITELIST" | sort -u | tr '\n' ' ')
echo "    $WL_COUNT barcodes, length(s): $WL_LENGTHS"

[ "$WL_LENGTHS" = "16 " ] \
    || die "whitelist barcodes are not a uniform 16 nt (got: $WL_LENGTHS).
  The signature assumes a 16 nt barcode, so this file will not work.
  File: $WHITELIST"

if [ "$WL_COUNT" != "$WHITELIST_EXPECTED_BARCODES" ]; then
    echo "    note: expected $WHITELIST_EXPECTED_BARCODES barcodes (10X v2)."
    echo "          ~3,000,000 would mean this is the v3 whitelist, which does"
    echo "          NOT match the 26-base signature in config.yaml -- v3 needs 28."
    echo "          A much smaller number usually means a truncated download."
fi

# --------------------------------------------------------------------------- #
say "2/6  Reference genome  ($REFERENCE_ACC)"
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
say "3/6  Region annotation  (GenBank -> GFF3)"
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
say "4/6  Reads  ($SRR)"
if [ -s "$FASTQ" ]; then
    echo "    $FASTQ exists, skipping."
elif [ -n "$MAX_SPOTS" ]; then
    # fastq-dump -X stops the transfer after N spots, so this is genuinely a
    # partial download rather than a full one that is then trimmed.
    #
    # fasterq-dump has no supported equivalent. It does carry --row-limit, but
    # its own source says "do not advertize row-limit" and the limit applies
    # PER THREAD, so the count you get back depends on --threads. Not something
    # to build a documented option on.
    command -v fastq-dump >/dev/null 2>&1 \
        || die "MAX_SPOTS needs fastq-dump, which ships with sra-tools next to
  fasterq-dump. Either install it, or unset MAX_SPOTS to download the full run."
    echo "    downloading the first $MAX_SPOTS spots only"
    fastq-dump -X "$MAX_SPOTS" --stdout "$SRR" > "$FASTQ.part" \
        || die "fastq-dump failed. If it rejected -X, your sra-tools is unusual;
  unset MAX_SPOTS to fall back to a full fasterq-dump download."
    mv "$FASTQ.part" "$FASTQ"
else
    echo "    downloading the FULL run; expect this to take a while."
    echo "    (set MAX_SPOTS=50000 to pull a small slice instead)"
    fasterq-dump "$SRR" --concatenate-reads --threads "$THREADS" \
                 --outdir "$DATA_DIR" --outfile "$(basename "$FASTQ")" --progress
fi

[ -s "$FASTQ" ] || die "no reads were downloaded for $SRR."
head -c1 "$FASTQ" | grep -q '@' \
    || die "downloaded reads are not FASTQ -- first line:
  $(head -1 "$FASTQ")"

FASTQ_READS=$(( $(wc -l < "$FASTQ") / 4 ))
echo "    $FASTQ_READS reads"

# If a spot limit was asked for, say whether it was honoured. A tool that
# ignored -X would otherwise just look like a slow download.
if [ -n "$MAX_SPOTS" ] && [ "$FASTQ_READS" -gt "$(( MAX_SPOTS * 2 ))" ]; then
    echo "    note: asked for $MAX_SPOTS spots but got $FASTQ_READS reads."
    echo "          More than one read per spot is normal; far more than that"
    echo "          suggests the spot limit was not applied."
fi

# --------------------------------------------------------------------------- #
say "5/6  Optional subsample (mapping only -- see MAX_SPOTS for the download)"
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
say "6/6  Map to the reference  (minimap2 -ax $MINIMAP_PRESET)"
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
    $(basename "$WHITELIST")   10X v2 barcode whitelist

Check these two lines in examples/eva71_sra/config.yaml match:

    reference_name: "$REF_NAME"
    whitelist: "$WHITELIST"

then:

    snakemake -s workflow/Snakefile --configfile examples/eva71_sra/config.yaml --cores $THREADS -n
    snakemake -s workflow/Snakefile --configfile examples/eva71_sra/config.yaml --cores $THREADS
EOF
