#!/usr/bin/env bash
# End-to-end build of every published product, for every registered assembly.
#
# This is the script that produced the Zenodo records. It is linear and verbose
# on purpose: each step prints what it is doing and stops the whole run on the
# first failure, so a half-built reference never silently reaches a release.
#
# Expect many hours and ~1 TB of scratch for the full matrix. To rebuild a
# single assembly, set GENOMES before calling:
#   GENOMES=GRCh38 scripts/run_all_pipeline.sh
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE/.."

GENOMES=${GENOMES:-"GRCh38 GRCh37 T2T-CHM13v2 GRCm38 GRCm39"}
THREADS=${THREADS:-16}
KMERS=${KMERS:-"24 36 50 100 150"}

step() { echo; echo "############ $* ############"; echo; }

# Per-assembly download provider, source-FASTA basename and extra mask_genome
# flags. Keeping these in one place is what went wrong before: the FASTA name
# was hardcoded as Homo_sapiens.* for all five assemblies.
provider_for() {
  case "$1" in
    T2T-CHM13v2) echo "ucsc" ;;
    *)           echo "ensembl" ;;
  esac
}

fasta_for() {
  case "$1" in
    GRCh38|GRCh37)   echo "data/$1/Homo_sapiens.$1.dna.primary_assembly.fa.gz" ;;
    GRCm38|GRCm39)   echo "data/$1/Mus_musculus.$1.dna.primary_assembly.fa.gz" ;;
    T2T-CHM13v2)     echo "data/$1/hs1.fa.gz" ;;
    *) echo "unknown assembly '$1'" >&2; return 1 ;;
  esac
}

# Mouse PAR intervals are approximate (see docs/PAR_coordinates.md), so XY
# masking there is gated behind an explicit opt-in -- both for the FASTA and for
# the annotation, which derives its Y-PAR drop set from the same intervals.
mask_flags_for() {
  case "$1" in
    GRCm38|GRCm39) echo "--allow-approximate-mouse-par" ;;
    *)             echo "" ;;
  esac
}

# 1. source genomes (+ human territory tables / ENCODE blacklist)
step "1. download source genomes"
for genome in $GENOMES; do
  scripts/download_references.sh -a "$genome" -o "data/$genome" \
    --preferred-source "$(provider_for "$genome")"
done

# 2. both SCC references (XX: whole-Y masked, XY: Y-PAR masked)
step "2. build sex-aware reference genomes"
for genome in $GENOMES; do
  # shellcheck disable=SC2046  # mask_flags_for is empty or a single flag
  python scripts/mask_genome.py --assembly "$genome" \
    --fasta "$(fasta_for "$genome")" \
    --complement both --outdir "refs/$genome" --gzip $(mask_flags_for "$genome")
done

# 2b. sex-aware annotation, filtered against the masks step 2 just wrote.
#     T2T-CHM13v2 has no Ensembl annotation and is skipped (see
#     config/annotations.tsv); the loop below is driven by that registry, so it
#     stays correct as releases are added.
step "2b. build sex-aware annotation (GTF)"
for genome in $GENOMES; do
  releases=$(scripts/download_annotations.sh -a "$genome" --release all --list 2>/dev/null \
             | awk 'NR>1 {print $3}' || true)
  [ -z "$releases" ] && { echo "no annotation registered for $genome; skipping"; continue; }
  for rel in $releases; do
    scripts/download_annotations.sh -a "$genome" -o "data/$genome" --release "$rel"
    gtf=$(scripts/download_annotations.sh -a "$genome" -o "data/$genome" --release "$rel" --print-path)
    # shellcheck disable=SC2046
    python scripts/make_sexaware_gtf.py --assembly "$genome" --gtf "$gtf" \
      --variant all --outdir "refs/$genome" --emit-t2g $(mask_flags_for "$genome")
  done
done

# 3. mappability (single + multi-read) for every k, per flavour
step "3. mappability tracks"
for genome in $GENOMES; do
  for c in XX XY; do
    bash scripts/run_mappability.sh -f "refs/$genome/$genome.$c.fa.gz" \
      -s "refs/$genome/$genome.$c.chrom.sizes" -o "map/$genome" -l "$genome.$c" -t "$THREADS"
  done
done

# 4. callable BEDs + ploidy files (extra products)
step "4. callable BEDs + ploidy configs"
for genome in $GENOMES; do
  mkdir -p callable
  for k in $KMERS; do
    python scripts/build_callable_beds.py \
      --bigwig "map/$genome/$genome.XY/k$k/$genome.XY.k$k.single_read.bw" \
      --threshold 1.0 --out "callable/$genome.XY.k$k.callable.bed"
  done
  # Ploidy files depend on the PAR + chrom.sizes only, not on k: build once.
  par_bed=$(awk -F'\t' -v a="$genome" '$1==a {print $3}' <(grep -v '^#' config/assemblies.tsv))
  # ploidy/ at the repo root, NOT config/ploidy: config/ holds hand-curated
  # registries, and a generated ploidy BED sitting next to the curated
  # config/par/*.bed is how the README came to list a config/ploidy/ that had
  # never existed. Generated products stay out of config/.
  [ "$par_bed" = "NA" ] || python scripts/make_ploidy_configs.py --par-bed "$par_bed" \
    --chrom-sizes "refs/$genome/$genome.XY.chrom.sizes" --label "$genome" --outdir ploidy
done

# 5. (optional, human) territory tables -- feed OUR sex-aware k100 multi_read
#    bigWig, NOT UCSC's precomputed track (that one is built on the UNMASKED
#    genome and would report mappability inside the PAR).
step "5. territory tables (GRCh38)"
if [ -f data/GRCh38/chromInfo.txt.gz ]; then
  mkdir -p territory
  python scripts/build_territory_tables.py \
    --chrom-info data/GRCh38/chromInfo.txt.gz --gap data/GRCh38/gap.txt.gz \
    --cytoband data/GRCh38/cytoBand.txt.gz \
    --blacklist data/GRCh38/ENCFF356LFX.hg38.blacklist.bed.gz \
    --umap-k100-bw map/GRCh38/GRCh38.XY/k100/GRCh38.XY.k100.multi_read.bw \
    --xy-regions config/xy_regions.txt \
    --out-prefix territory/GRCh38.XY
fi

# 6. release manifest (checksums + provenance for every published file)
step "6. release manifest"
scripts/build_manifest.sh -o MANIFEST.tsv

echo
echo "done. publishable products: refs/*/  map/*/  MANIFEST.tsv"
