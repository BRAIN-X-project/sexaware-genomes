#!/usr/bin/env bash
# Download a source genome FASTA for a registered assembly, choosing among
# providers by preference with graceful fallback.
#
#   --preferred-source ensembl,ucsc,ncbi   (default)
#       Try each source in order; use the first that has a URL for the assembly
#       in config/sources.tsv.
#
# Contig naming is normalised downstream by mask_genome.py, so any source works
# with the same PAR BEDs:
#   * Ensembl -> 1/X/Y            (Ensembl already hard-masks the chrY PAR with Ns;
#                                   the XY step is then idempotent. XX still masks Y.)
#   * UCSC    -> chr1/chrX/chrY
#   * NCBI    -> RefSeq accessions, renamed here to chrN via the assembly_report.
#
# For human it also fetches the UCSC territory tables + ENCODE blacklist used by
# build_territory_tables.py (no clean Ensembl equivalent for gap/cytoBand).
set -euo pipefail

log() { echo "[download_references $(date +%H:%M:%S)] $*" >&2; }

HERE="$(cd "$(dirname "$0")" && pwd)"
SOURCES="$HERE/../config/sources.tsv"

usage() {
  cat <<EOF
Usage: $0 -a ASSEMBLY -o OUTDIR [--preferred-source LIST] [--no-rename]
  -a  assembly: GRCh38 GRCh37 T2T-CHM13v2 GRCm39 GRCm38
  -o  output directory
  --preferred-source  comma list, default "ensembl,ucsc,ncbi" (order = preference)
  --no-rename         keep NCBI accession contig names (skip chrN rename)
EOF
  exit 1
}

ASSEMBLY=""; OUTDIR=""; PREF="ensembl,ucsc,ncbi"; RENAME=1
while [ $# -gt 0 ]; do case "$1" in
  -a) ASSEMBLY=$2; shift 2;;
  -o) OUTDIR=$2; shift 2;;
  --preferred-source) PREF=$2; shift 2;;
  --no-rename) RENAME=0; shift;;
  -h|--help) usage;;
  *) usage;;
esac; done
[ -z "$ASSEMBLY" ] || [ -z "$OUTDIR" ] && usage
mkdir -p "$OUTDIR"

# Look up the URL for a given assembly+source in sources.tsv (empty if absent).
lookup_url() { awk -F'\t' -v a="$1" -v s="$2" '$1==a && $2==s {print $4}' <(grep -v '^#' "$SOURCES"); }
lookup_rename() { awk -F'\t' -v a="$1" -v s="$2" '$1==a && $2==s {print $3}' <(grep -v '^#' "$SOURCES"); }

# Pick the first preferred source that exists for this assembly.
CHOSEN=""; URL=""; NEEDS_RENAME=0
IFS=',' read -ra ORDER <<< "$PREF"
for src in "${ORDER[@]}"; do
  src=$(echo "$src" | tr -d '[:space:]')
  candidate=$(lookup_url "$ASSEMBLY" "$src")
  if [ -n "$candidate" ]; then
    CHOSEN=$src; URL=$candidate; NEEDS_RENAME=$(lookup_rename "$ASSEMBLY" "$src")
    break
  fi
  log "source '$src' has no URL for $ASSEMBLY; trying next"
done
[ -z "$URL" ] && { echo "No source in '$PREF' provides $ASSEMBLY (see config/sources.tsv)"; exit 1; }

log "assembly=$ASSEMBLY  chosen source=$CHOSEN"
log "FASTA URL: $URL"
wget -c -P "$OUTDIR" "$URL"
FASTA="$OUTDIR/$(basename "$URL")"
log "downloaded $(basename "$URL")"
USE_FASTA="$FASTA"

# ---- NCBI: rename RefSeq accessions -> chrN so chrX/chrY PAR BED matches.
if [ "$NEEDS_RENAME" = "1" ] && [ "$RENAME" = "1" ]; then
  REPORT_URL="${URL%_genomic.fna.gz}_assembly_report.txt"
  log "assembly_report: $REPORT_URL"
  wget -c -O "$OUTDIR/${ASSEMBLY}.assembly_report.txt" "$REPORT_URL"
  log "renaming RefSeq accessions to chrN"
  python3 - "$FASTA" "$OUTDIR/${ASSEMBLY}.assembly_report.txt" "$OUTDIR/${ASSEMBLY}.chrNamed.fa.gz" <<'PY'
import sys, gzip
fasta, report, out = sys.argv[1:4]
name_map = {}
for line in open(report, encoding="utf-8", errors="replace"):
    if line.startswith("#") or not line.strip():
        continue
    f = line.rstrip("\n").split("\t")
    role, molecule, refseq = f[1], f[2], f[6]
    if refseq in ("na", ""):
        continue
    if role == "assembled-molecule":
        chrom = "chrM" if molecule in ("MT", "Mitochondrion") else f"chr{molecule}"
    else:
        chrom = refseq  # keep unplaced/unlocalized as their accession
    name_map[refseq] = chrom
op = gzip.open if fasta.endswith(".gz") else open
with op(fasta, "rt") as fh, gzip.open(out, "wt") as oh:
    for line in fh:
        if line.startswith(">"):
            oh.write(">%s\n" % name_map.get(line[1:].split()[0], line[1:].split()[0]))
        else:
            oh.write(line)
print("wrote", out)
PY
  USE_FASTA="$OUTDIR/${ASSEMBLY}.chrNamed.fa.gz"
  log "chr-named FASTA: $USE_FASTA"
fi

log "==> feed this FASTA to mask_genome.py:  $USE_FASTA"

# ---- Human auxiliary tables for territory accounting (UCSC only source).
# These describe assembly STRUCTURE (gaps, cytobands) and problem regions
# (ENCODE blacklist) and are independent of the sex-aware masking, so they are
# safe to reuse. We deliberately DO NOT download UCSC's precomputed Umap
# mappability bigWig: it is built on the unmasked genome and would report
# mappability inside the PAR (as if both X and Y copies existed), which defeats
# the purpose here. build_territory_tables.py must instead be fed the sex-aware
# k100 multi_read bigWig produced by run_mappability.sh on the masked reference.
case "$ASSEMBLY" in GRCh38) DB=hg38;; GRCh37) DB=hg19;; *) DB="";; esac
if [ -n "$DB" ]; then
  log "UCSC structural tables (gap/cytoBand/chromInfo) + ENCODE blacklist ($DB)"
  for t in chromInfo gap cytoBand; do
    wget -c -P "$OUTDIR" "https://hgdownload.soe.ucsc.edu/goldenPath/$DB/database/$t.txt.gz"
  done
  if [ "$DB" = "hg38" ]; then
    curl -L -o "$OUTDIR/ENCFF356LFX.hg38.blacklist.bed.gz" \
      "https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz"
  fi
fi
log "done -> $OUTDIR"
