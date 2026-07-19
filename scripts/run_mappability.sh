#!/usr/bin/env bash
# Compute single-read and multi-read mappability for a (masked) genome with GenMap.
#
# GenMap is the default, modern engine. For k-mer uniqueness with exact matching
# (-E 0) it is numerically equivalent to Umap; it is far faster and trivial to
# containerise. Semantics reproduced to match Umap/Hoffman tracks:
#
#   single-read mappability(pos) = 1 if the k-mer STARTING at pos is unique, else
#                                  1/frequency  (GenMap's native per-k-mer value)
#   multi-read  mappability(pos) = mean over the k k-mers OVERLAPPING pos of the
#                                  single-read uniqueness indicator (Umap's
#                                  definition), computed by genmap_to_umap_tracks.py
#
# Outputs, per k, under $OUTDIR/<label>/k<k>/ :
#   <label>.k<k>.single_read.bw
#   <label>.k<k>.multi_read.bw
#
# Requires: genmap, bedGraphToBigWig (UCSC), python3 (+numpy), the sizes file.
set -euo pipefail

log() { echo "[run_mappability $(date +%H:%M:%S)] $*" >&2; }

usage() {
  cat <<EOF
Usage: $0 -f MASKED_FASTA -s CHROM_SIZES -o OUTDIR -l LABEL [-k "24 36 50 100 150"] [-t THREADS] [-S SAMPLING]
  -f  masked FASTA (output of mask_genome.py)
  -s  chrom.sizes for that FASTA
  -o  output directory
  -l  label used in filenames, e.g. GRCh38.XY
  -k  space-separated k-mer sizes (default: "24 36 50 100 150")
  -t  threads for genmap (default: 4)
  -S  GenMap index sampling (memory saver). Empty = default (fastest, most RAM).
      For whole human/mouse genomes on <~32 GB RAM use -S 20 (up to 64) to avoid
      the "Create fwd Index" segfault. Does NOT change the mappability values,
      only index size/build memory and speed.
EOF
  exit 1
}

KMERS="24 36 50 100 150"
THREADS=4
SAMPLING=""
FASTA=""; SIZES=""; OUTDIR=""; LABEL=""
while getopts "f:s:o:l:k:t:S:h" opt; do
  case $opt in
    f) FASTA=$OPTARG ;; s) SIZES=$OPTARG ;; o) OUTDIR=$OPTARG ;;
    l) LABEL=$OPTARG ;; k) KMERS=$OPTARG ;; t) THREADS=$OPTARG ;;
    S) SAMPLING=$OPTARG ;;
    *) usage ;;
  esac
done
[ -z "$FASTA" ] || [ -z "$SIZES" ] || [ -z "$OUTDIR" ] || [ -z "$LABEL" ] && usage

HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$OUTDIR"
INDEX="$OUTDIR/${LABEL}.genmap_index"

log "label=$LABEL  fasta=$FASTA"
log "sizes=$SIZES  kmers=[$KMERS]  threads=$THREADS  outdir=$OUTDIR"

# GenMap will not read a .gz and needs a valid FASTA extension. If the reference
# is gzipped, decompress to a temporary .fa just for indexing; the map step below
# uses the FM-index (-I), not the FASTA, so we can delete the temp right after.
GENMAP_FASTA="$FASTA"
TMP_FASTA=""
case "$FASTA" in
  *.gz)
    TMP_FASTA="$OUTDIR/${LABEL}.decompressed.fa"
    log "GenMap needs an uncompressed FASTA; decompressing -> $TMP_FASTA"
    zcat "$FASTA" > "$TMP_FASTA"
    GENMAP_FASTA="$TMP_FASTA"
    ;;
  *.fa|*.fasta|*.fna|*.fas|*.faa|*.fsa) : ;;
  *) log "WARNING: '$FASTA' lacks a GenMap-recognised extension; indexing may fail." ;;
esac

# genmap requires the index directory NOT to pre-exist.
rm -rf "$INDEX"
SAMPLING_OPT=""
if [ -n "$SAMPLING" ]; then SAMPLING_OPT="-S $SAMPLING"; fi
log "step 1/2: building GenMap FM-index -> $INDEX ${SAMPLING_OPT:+(sampling $SAMPLING)}"
log "note: GenMap needs ~10x the genome size in RAM (~30 GB for a whole human"
log "      genome). If it segfaults at 'Create fwd Index', add -S 20 and/or give"
log "      the machine more RAM (on WSL, set memory=48GB in %UserProfile%\\.wslconfig)."
# shellcheck disable=SC2086
genmap index -F "$GENMAP_FASTA" -I "$INDEX" $SAMPLING_OPT
log "index done"
if [ -n "$TMP_FASTA" ]; then rm -f "$TMP_FASTA"; log "removed temporary decompressed FASTA"; fi

for K in $KMERS; do
  KDIR="$OUTDIR/$LABEL/k${K}"
  mkdir -p "$KDIR"
  RAW="$KDIR/${LABEL}.k${K}.genmap"
  log "k=$K: computing per-k-mer mappability (genmap map -E 0)"
  # -bg = bedGraph of the native per-k-mer mappability (1/frequency), exact (-E 0)
  genmap map -K "$K" -E 0 -I "$INDEX" -O "$RAW" -bg -T "$THREADS"

  log "k=$K: deriving single/multi-read tracks (Umap semantics)"
  python3 "$HERE/genmap_to_umap_tracks.py" \
    --genmap-bedgraph "${RAW}.bedgraph" \
    --chrom-sizes "$SIZES" \
    --kmer "$K" \
    --single-out "$KDIR/${LABEL}.k${K}.single_read.bedgraph" \
    --multi-out  "$KDIR/${LABEL}.k${K}.multi_read.bedgraph"

  log "k=$K: bedGraph -> bigWig"
  bedGraphToBigWig "$KDIR/${LABEL}.k${K}.single_read.bedgraph" "$SIZES" "$KDIR/${LABEL}.k${K}.single_read.bw"
  bedGraphToBigWig "$KDIR/${LABEL}.k${K}.multi_read.bedgraph"  "$SIZES" "$KDIR/${LABEL}.k${K}.multi_read.bw"
  gzip -f "$KDIR/${LABEL}.k${K}.single_read.bedgraph" "$KDIR/${LABEL}.k${K}.multi_read.bedgraph" "${RAW}.bedgraph"
  log "k=$K: wrote ${LABEL}.k${K}.{single,multi}_read.bw"
done

log "step 2/2 complete -> $OUTDIR/$LABEL"
