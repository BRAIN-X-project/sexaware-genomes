#!/usr/bin/env bash
# Optional: reproduce classic Umap (Hoffman lab) mappability for continuity with
# UCSC/ENCODE tracks. Legacy Umap needs Python 2.7 + bowtie1 (see
# envs/environment.umap_legacy.yml). GenMap (run_mappability.sh) is the default;
# use this only when you must match pre-existing Umap bigWigs exactly.
#
# This is a parameterised, path-free rewrite of the original 1B_run_umap_command.sh.
set -euo pipefail

usage() {
  cat <<EOF
Usage: $0 -g GENOME_FASTA -s CHROM_SIZES -o OUTDIR -l LABEL -u UMAP_DIR -p PY27_ENV \\
          [-k "24 36 50 100 150"] [-j JOBS] [-q QUEUE]
  -g  masked genome FASTA (from mask_genome.py)
  -s  chrom.sizes
  -o  Umap working/output directory
  -l  label, e.g. GRCh38.XY
  -u  path to the umap/ checkout containing ubismap.py
  -p  path to the Python 2.7 conda env (bin/python, bin/bowtie-build inside)
  -k  k-mer sizes (default "24 36 50 100 150")
  -j  parallel jobs for uint8->bed (default 8)
  -q  cluster queue name passed to ubismap.py (default all.q; local runs ignore it)
EOF
  exit 1
}

KMERS="24 36 50 100 150"; JOBS=8; QUEUE="all.q"
GENOME=""; SIZES=""; OUTDIR=""; LABEL=""; UMAP_DIR=""; PY27=""
while getopts "g:s:o:l:u:p:k:j:q:h" opt; do
  case $opt in
    g) GENOME=$OPTARG ;; s) SIZES=$OPTARG ;; o) OUTDIR=$OPTARG ;; l) LABEL=$OPTARG ;;
    u) UMAP_DIR=$OPTARG ;; p) PY27=$OPTARG ;; k) KMERS=$OPTARG ;; j) JOBS=$OPTARG ;;
    q) QUEUE=$OPTARG ;; *) usage ;;
  esac
done
for v in GENOME SIZES OUTDIR LABEL UMAP_DIR PY27; do [ -z "${!v}" ] && usage; done

PYBIN="$PY27/bin/python"
BOWTIE_BUILD="$PY27/bin/bowtie-build"
mkdir -p "$OUTDIR"
NCHR=$(wc -l < "$SIZES")

for K in $KMERS; do
  KTAG="k${K}"
  KDIR="$OUTDIR/$LABEL/$KTAG"
  GMAP="$KDIR/globalmap_${KTAG}to${KTAG}"
  BG="$KDIR/bedGraph"
  mkdir -p "$KDIR" "$BG"

  echo ">> [Umap] ubismap k=$K ($LABEL)"
  # ubismap.py builds the bowtie index, runs get_kmers.py + run_bowtie.py + unify.
  # On HPC it writes a Slurm script; for local runs feed each generated block.
  ( cd "$UMAP_DIR" && "$PYBIN" ubismap.py "$GENOME" "$SIZES" "$KDIR" "$QUEUE" "$BOWTIE_BUILD" \
      --kmers "$K" -write_script "$KDIR/umap_${KTAG}.generated.sh" )
  echo "   Review/adapt $KDIR/umap_${KTAG}.generated.sh (block-by-block for local), then run it."

  echo ">> [Umap] uint8 -> bedGraph k=$K ($LABEL)"
  seq 1 "$NCHR" | xargs -P "$JOBS" -I{} "$PYBIN" "$UMAP_DIR/uint8_to_bed_parallel.py" \
      "$GMAP" "$BG" "$LABEL" -chrsize_path "$SIZES" -wiggle -kmers "$KTAG" -job_id {}

  echo ">> [Umap] wig -> bigWig k=$K ($LABEL)"
  WIG="$BG/${LABEL}.${KTAG}.MultiReadMappability.wig"
  if ls "$BG"/*.wg.gz >/dev/null 2>&1; then zcat "$BG"/*.wg.gz > "$WIG"; else cat "$BG"/*.wg > "$WIG"; fi
  wigToBigWig "$WIG" "$SIZES" "$BG/${LABEL}.${KTAG}.MultiReadMappability.bw"
done
echo ">> [Umap] done: $OUTDIR/$LABEL"
