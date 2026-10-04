#!/usr/bin/env bash
# Download an annotation GTF for a registered assembly + release.
#
# Companion to download_references.sh: that script fetches the FASTA, this one
# fetches the matching GTF. Both read a registry under config/ and print the
# exact local path to hand to the next step.
#
#   scripts/download_annotations.sh -a GRCh38 -o data/GRCh38                 # latest
#   scripts/download_annotations.sh -a GRCh38 -o data/GRCh38 --release 110
#   scripts/download_annotations.sh -a GRCh38 -o data/GRCh38 --release all
#   scripts/download_annotations.sh -a GRCh38 --list                         # no download
#
# Checksums: Ensembl publishes a CHECKSUMS file per directory in BSD `sum`
# format (`<checksum> <1K-blocks> <filename>`), which GNU coreutils `sum`
# reproduces by default. Ensembl publishes no md5 for GTFs, so we compute one
# ourselves and leave it next to the file for downstream verification.
set -euo pipefail

log() { echo "[download_annotations $(date +%H:%M:%S)] $*" >&2; }
die() { echo "[download_annotations] ERROR: $*" >&2; exit 1; }

HERE="$(cd "$(dirname "$0")" && pwd)"
ANNOTATIONS="$HERE/../config/annotations.tsv"

usage() {
  cat <<EOF
Usage: $0 -a ASSEMBLY [-o OUTDIR] [options]
  -a  assembly: GRCh38 GRCh37 GRCm39 GRCm38   (T2T-CHM13v2 has no Ensembl annotation)
  -o  output directory (required unless --list)
  --release N|latest|all  which release(s); default latest
  --flavour primary|chr   GTF flavour; default primary (matches dna.primary_assembly)
  --source NAME           annotation provider; default ensembl
  --list                  print the matching registry rows and exit
  --print-path            print the local GTF path for the match and exit (no download)
  --no-verify             skip CHECKSUMS verification
EOF
  exit 1
}

ASSEMBLY=""; OUTDIR=""; RELEASE="latest"; FLAVOUR="primary"; SOURCE="ensembl"
LIST=0; PRINT_PATH=0; VERIFY=1
while [ $# -gt 0 ]; do case "$1" in
  -a) ASSEMBLY=$2; shift 2;;
  -o) OUTDIR=$2; shift 2;;
  --release) RELEASE=$2; shift 2;;
  --flavour) FLAVOUR=$2; shift 2;;
  --source) SOURCE=$2; shift 2;;
  --list) LIST=1; shift;;
  --print-path) PRINT_PATH=1; shift;;
  --no-verify) VERIFY=0; shift;;
  -h|--help) usage;;
  *) usage;;
esac; done
[ -z "$ASSEMBLY" ] && usage
[ "$LIST" = "0" ] && [ "$PRINT_PATH" = "0" ] && [ -z "$OUTDIR" ] && usage

# Select the matching rows. 'latest' keeps is_latest==1, 'all' keeps every
# release, a number keeps that release.
select_rows() {
  awk -F'\t' -v a="$ASSEMBLY" -v s="$SOURCE" -v f="$FLAVOUR" -v r="$RELEASE" '
    $1==a && $2==s && $4==f {
      if (r=="all") print
      else if (r=="latest") { if ($5=="1") print }
      else if ($3==r) print
    }' <(grep -v '^#' "$ANNOTATIONS" | tail -n +2)
}

ROWS=$(select_rows)
if [ -z "$ROWS" ]; then
  die "no annotation row for assembly=$ASSEMBLY source=$SOURCE flavour=$FLAVOUR release=$RELEASE
     (see config/annotations.tsv; T2T-CHM13v2 is intentionally unannotated)"
fi

if [ "$LIST" = "1" ]; then
  printf 'assembly\tsource\trelease\tflavour\tis_latest\trestrict_to_genome\turl\n'
  echo "$ROWS" | cut -f1-7
  exit 0
fi

if [ "$PRINT_PATH" = "1" ]; then
  echo "$ROWS" | while IFS=$'\t' read -r _a _s _r _f _l _rg url _rest; do
    echo "${OUTDIR:-.}/$(basename "$url")"
  done
  exit 0
fi

mkdir -p "$OUTDIR"

echo "$ROWS" | while IFS=$'\t' read -r asm src rel flav latest restrict url checks notes; do
  # Guard rail: these two flavours are never publishable here, whatever the TSV
  # says. chr_patch_hapl_scaff adds contigs that do not exist in the
  # primary_assembly FASTA; abinitio is not the genebuild. Refuse BEFORE
  # spending a 50-150 MB download that the validator would reject anyway.
  case "$url" in
    *chr_patch_hapl_scaff*) die "refusing $url: chr_patch_hapl_scaff annotates contigs absent from primary_assembly" ;;
    *abinitio*)             die "refusing $url: abinitio is ab-initio prediction, not the genebuild" ;;
  esac

  FILE="$OUTDIR/$(basename "$url")"
  log "assembly=$asm source=$src release=$rel flavour=$flav$([ "$latest" = "1" ] && echo ' (latest)')"
  log "GTF URL: $url"
  wget -c -P "$OUTDIR" "$url"
  log "downloaded $(basename "$FILE") ($(du -h "$FILE" | cut -f1))"

  if [ "$VERIFY" = "1" ]; then
    log "verifying against $checks"
    # BSD sum fields are space-padded, so compare NUMERICALLY, not as strings.
    EXPECTED=$(curl -fsSL "$checks" | awk -v n="$(basename "$FILE")" '$3==n {print $1" "$2}')
    if [ -z "$EXPECTED" ]; then
      log "WARNING: $(basename "$FILE") not listed in CHECKSUMS; skipping verification"
    else
      ACTUAL=$(sum "$FILE")
      exp_a=$(echo "$EXPECTED" | awk '{print $1+0}'); exp_b=$(echo "$EXPECTED" | awk '{print $2+0}')
      act_a=$(echo "$ACTUAL"   | awk '{print $1+0}'); act_b=$(echo "$ACTUAL"   | awk '{print $2+0}')
      if [ "$exp_a" -eq "$act_a" ] && [ "$exp_b" -eq "$act_b" ]; then
        log "CHECKSUMS OK (sum $act_a $act_b)"
      else
        die "CHECKSUMS MISMATCH for $(basename "$FILE"): expected '$exp_a $exp_b', got '$act_a $act_b'.
     The download is incomplete or corrupt; delete it and retry."
      fi
    fi
  fi

  # Ensembl publishes no md5 for GTFs; leave one so downstream users (and
  # build_manifest.sh) can verify the bytes they actually hold.
  ( cd "$OUTDIR" && md5sum "$(basename "$FILE")" > "$(basename "$FILE").md5" )
  log "wrote $(basename "$FILE").md5"
  [ "$restrict" = "1" ] && log "NOTE: this release annotates contigs absent from primary_assembly; \
make_sexaware_gtf.py must be run with --on-unknown-seqname drop"

  log "==> feed this GTF to make_sexaware_gtf.py:  $FILE"
done
