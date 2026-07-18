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
Usage: $0 -f MASKED_FASTA -s CHROM_SIZES -o OUTDIR -l LABEL [-k "24 36 50 100 150"] [-t THREADS]
  -f  masked FASTA (output of mask_genome.py)
  -s  chrom.sizes for that FASTA
  -o  output directory
  -l  label used in filenames, e.g. GRCh38.XY
  -k  space-separated k-mer sizes (default: "24 36 50 100 150")
  -t  threads for genmap (default: 4)
EOF
  exit 1
}

KMERS="24 36 50 100 150"
THREADS=4
FASTA=""; SIZES=""; OUTDIR=""; LABEL=""
while getopts "f:s:o:l:k:t:h" opt; do
  case $opt in
    f) FASTA=$OPTARG ;; s) SIZES=$OPTARG ;; o) OUTDIR=$OPTARG ;;
    l) LABEL=$OPTARG ;; k) KMERS=$OPTARG ;; t) THREADS=$OPTARG ;;
    *) usage ;;
  esac
done
[ -z "$FASTA" ] || [ -z "$SIZES" ] || [ -z "$OUTDIR" ] || [ -z "$LABEL" ] && usage

HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$OUTDIR"
INDEX="$OUTDIR/${LABEL}.genmap_index"

log "label=$LABEL  fasta=$FASTA"
log "sizes=$SIZES  kmers=[$KMERS]  threads=$THREADS  outdir=$OUTDIR"

# genmap requires the index directory NOT to pre-exist.
rm -rf "$INDEX"
log "step 1/2: building GenMap FM-index -> $INDEX"
genmap index -F "$FASTA" -I "$INDEX"
log "index done"

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
