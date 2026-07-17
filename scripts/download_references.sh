#!/usr/bin/env bash
# Download a source genome FASTA for a registered assembly, plus (for human) the
# UCSC territory tables and ENCODE blacklist used by build_territory_tables.py.
# Generalises the original 0B_download_references.sh.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REGISTRY="$HERE/../config/assemblies.tsv"

usage() { echo "Usage: $0 -a ASSEMBLY -o OUTDIR   (ASSEMBLY in: GRCh38 GRCh37 T2T-CHM13v2 GRCm39 GRCm38)"; exit 1; }
ASSEMBLY=""; OUTDIR=""
while getopts "a:o:h" opt; do case $opt in a) ASSEMBLY=$OPTARG;; o) OUTDIR=$OPTARG;; *) usage;; esac; done
[ -z "$ASSEMBLY" ] || [ -z "$OUTDIR" ] && usage
mkdir -p "$OUTDIR"

URL=$(awk -F'\t' -v a="$ASSEMBLY" '$1==a{print $8}' <(grep -v '^#' "$REGISTRY"))
[ -z "$URL" ] && { echo "Unknown assembly $ASSEMBLY"; exit 1; }

echo ">> FASTA: $URL"
wget -c -P "$OUTDIR" "$URL"

# Human-only auxiliary tracks for territory tables (edit UCSC db as needed).
case "$ASSEMBLY" in
  GRCh38) DB=hg38 ;; GRCh37) DB=hg19 ;; *) DB="" ;;
esac
if [ -n "$DB" ]; then
  echo ">> UCSC tables + ENCODE blacklist for $DB"
  for t in chromInfo gap cytoBand; do
    wget -c -P "$OUTDIR" "https://hgdownload.soe.ucsc.edu/goldenPath/$DB/database/$t.txt.gz"
  done
  if [ "$DB" = "hg38" ]; then
    wget -c -P "$OUTDIR" "https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k100.Umap.MultiTrackMappability.bw"
    curl -L -o "$OUTDIR/ENCFF356LFX.hg38.blacklist.bed.gz" \
      "https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz"
  fi
fi
echo ">> done -> $OUTDIR"
