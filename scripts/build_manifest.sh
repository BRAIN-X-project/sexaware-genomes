#!/usr/bin/env bash
# Build the release manifest: one row per published file, with the checksums a
# downloader can actually verify.
#
# The descriptive columns are read from the *.metadata.json each producer wrote
# next to its output, so the manifest cannot drift from the thing it describes.
# The checksums are recomputed here from the bytes on disk rather than copied
# out of the JSON -- the whole point of the file is to be an independent check,
# and Zenodo publishes an md5 per file in its own API, so these cross directly
# against their metadata.
#
#   scripts/build_manifest.sh [-o MANIFEST.tsv] [-d refs -d map]
#                             [--no-sha256] [--include-intermediates]
#
# --no-sha256 skips the second digest; md5 alone still cross-checks against
# Zenodo. Both digests come from a SINGLE read of each file, so the default
# costs no extra I/O over --no-sha256 -- running md5sum then sha256sum would
# read all 47 GB twice.
#
# --include-intermediates keeps GenMap's index directories in the walk. They are
# pruned by default: unpublished build scratch, 2.4 GB x 16, which only ever
# shows up as kind=other rows in a file that is supposed to describe a release.
set -euo pipefail

log() { echo "[build_manifest $(date +%H:%M:%S)] $*" >&2; }

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE/.."

usage() {
  cat <<EOF
Usage: $0 [-o OUT] [-d DIR]... [--no-sha256] [--include-intermediates]
  -o  output TSV (default MANIFEST.tsv; '-' for stdout)
  -d  directory to walk, repeatable (default: refs map)
  --no-sha256  omit the sha256 column (faster on multi-GB FASTAs)
  --include-intermediates  also hash aligner/GenMap index internals (see below)
EOF
  exit 1
}

OUT="MANIFEST.tsv"; DIRS=(); SHA=1; INTERMEDIATES=0
while [ $# -gt 0 ]; do case "$1" in
  -o) OUT=$2; shift 2;;
  -d) DIRS+=("$2"); shift 2;;
  --no-sha256) SHA=0; shift;;
  --include-intermediates) INTERMEDIATES=1; shift;;
  -h|--help) usage;;
  *) usage;;
esac; done
[ ${#DIRS[@]} -eq 0 ] && DIRS=(refs map)

# GenMap index directories are build scratch, never published: 2.4 GB per
# `index.lf.drv` x 16 of them under map/ means a default run otherwise spends
# ~40 GB of double hashing on files that must not appear in a release manifest
# at all (they land as kind=other, which is the giveaway). They are gitignored
# as *_genmap_index/ for the same reason. --include-intermediates is the escape
# hatch for auditing a working tree rather than describing a release.
PRUNE=()
if [ "$INTERMEDIATES" = 0 ]; then
  PRUNE=( '(' -type d -name '*genmap_index' -prune ')' -o )
fi

# Pull the descriptive columns out of one metadata.json. Keys differ between
# producers (mask_genome.py writes 'complement', make_sexaware_gtf.py writes
# 'variant'), so ask for both and take whichever is there.
read_meta() {
  python3 - "$1" "$2" "$3" <<'PY'
import json, re, sys
kind, fname = sys.argv[2], sys.argv[3]
try:
    d = json.load(open(sys.argv[1]))
except Exception:
    print("\t".join([""] * 9)); raise SystemExit
def g(*keys, default=""):
    for k in keys:
        v = d.get(k)
        if v not in (None, ""):
            return str(v)
    return default
# Notes a downloader needs before picking a file. GRCh37's XY-noYPAR is the
# live case: Ensembl annotates nothing on that assembly's Y PAR, so the file is
# content-identical to the standard XY one and nobody should wonder which to use.
notes = []
if d.get("identical_to_standard_xy"):
    notes.append("content-identical to the XY variant (this assembly has no Y-PAR annotation)")
if d.get("restrict_to_genome"):
    notes.append(f"restricted to the primary assembly: {d.get('genes_dropped_off_genome', 0)} "
                 "gene(s) on absent scaffolds removed (not sex-aware filtering)")
if d.get("par_status") and d["par_status"] != "verified":
    notes.append(f"PAR intervals are {d['par_status']}")
if d.get("genes_partially_overlapping"):
    notes.append(f"{d['genes_partially_overlapping']} gene(s) straddle the mask boundary "
                 "and were dropped whole")

tool = g("tool", default="mask_genome.py")
source = g("annotation_source", default="-")
release = g("annotation_release", default="-")
flavour = g("annotation_flavour", default="-")

# A mappability bigWig legitimately inherits assembly and variant from the
# genome product it was computed on, but NOT its provenance: that JSON was
# written by mask_genome.py and says nothing about GenMap or about k. Carrying
# its tool field over would credit the wrong producer for 80 of 265 rows.
if kind == "mappability":
    tool, source, release, flavour = "run_mappability.sh (GenMap)", "-", "-", "-"
    m = re.search(r"\.k(\d+)\.(single_read|multi_read)\.", fname)
    if m:
        notes.insert(0, f"k={m.group(1)}, {m.group(2)}")
    notes.append("computed on the sex-aware FASTA of the same variant")

print("\t".join([
    g("assembly"),
    g("variant", "complement"),
    source,
    release,
    flavour,
    tool,
    g("tool_version", default="-"),
    g("repo_revision", default="-"),
    "; ".join(notes) if notes else "-",
]))
PY
}

# Map a published file to the metadata.json describing it. Both producers name
# their sidecar after the output with the extension replaced, so strip the
# known suffixes rather than guessing.
# Registered assembly keys, used to validate the first filename token. A file
# whose leading token is not a registered assembly gets NO metadata rather than
# a guess -- a wrong assembly in a release manifest is worse than an empty cell.
ASSEMBLIES=$(awk -F'\t' 'NF && $1 !~ /^#/ {print $1}' config/assemblies.tsv)

# Derive the <ASM>.<XX|XY>[.noYPAR][.<source>-<release>] stem that says which
# PRODUCT a file belongs to.
#
# Parsed forwards, not by stripping suffixes: the suffixes are not a closed set
# (.fa.gz, .fa.md5, .chrom.sizes, .masked.bed, .k100.multi_read.bw,
# .dropped_genes.tsv, .t2g.tsv, ...) and any strip-list rots the moment a
# producer emits a new companion file. Consuming stem tokens and stopping at the
# first one that cannot belong to a stem needs no such list.
stem_of() {
  local base; base=$(basename "$1")
  local -a tok; IFS='.' read -ra tok <<< "$base"
  grep -qxF -- "${tok[0]}" <<< "$ASSEMBLIES" || return 1
  case "${tok[1]:-}" in XX|XY) ;; *) return 1;; esac
  local stem="${tok[0]}.${tok[1]}" i=2
  [ "${tok[$i]:-}" = "noYPAR" ] && { stem="$stem.noYPAR"; i=$((i + 1)); }
  # an annotation stem carries a <source>-<release> token, e.g. ensembl-116;
  # .chrom / .k50 / .masked / .fa do not match, so they correctly end the stem.
  case "${tok[$i]:-}" in *-[0-9]*) stem="$stem.${tok[$i]}";; esac
  printf '%s\n' "$stem"
}

# Find the metadata.json for a file's product. Companions and sidecars inherit
# from the product they describe, which is what takes the manifest from 37 of
# 265 rows described to all of them: a .fa.md5, a .chrom.sizes, a .t2g.tsv and a
# mappability bigWig all belong to a product that DID write a metadata.json.
# The bigWigs live under map/<ASM>/..., so refs/<ASM> is searched too.
meta_for() {
  local f=$1 stem asm d
  stem=$(stem_of "$f") || { echo ""; return; }
  asm=${stem%%.*}
  for d in "$(dirname "$f")" "refs/$asm"; do
    [ -f "$d/$stem.metadata.json" ] && { printf '%s\n' "$d/$stem.metadata.json"; return; }
  done
  echo ""
}

kind_of() {
  case "$1" in
    *.gtf.gz|*.gtf)       echo annotation ;;
    *.fa.gz|*.fa)         echo genome ;;
    *.bw|*.bigWig)        echo mappability ;;
    *.chrom.sizes)        echo chrom_sizes ;;
    *.masked.bed)         echo mask_bed ;;
    *.t2g.tsv)            echo t2g ;;
    *.dropped_genes.tsv)  echo dropped_genes ;;
    *.metadata.json)      echo provenance ;;
    *.md5)                echo checksum_sidecar ;;
    *)                    echo other ;;
  esac
}

tmp=$(mktemp); trap 'rm -f "$tmp"' EXIT
hdr="filename\tbytes\tmd5"
[ "$SHA" = 1 ] && hdr="$hdr\tsha256"
hdr="$hdr\tkind\tassembly\tvariant\tsource\trelease\tflavour\ttool\ttool_version\trepo_revision\tnotes\tpath"
printf "$hdr\n" > "$tmp"

n=0
# -print0 and read -d '' so a space in a path cannot split a row.
while IFS= read -r -d '' f; do
  case "$f" in *.fai|*.log) continue;; esac
  bytes=$(stat -c %s "$f")
  # One read, both digests. md5sum followed by sha256sum reads every byte twice,
  # which on the published bigWigs alone (~47 GB) is minutes of pure I/O.
  hashes=$(python3 - "$f" "$SHA" <<'PYHASH'
import hashlib, sys
path, want_sha = sys.argv[1], sys.argv[2] == "1"
m, sh = hashlib.md5(), hashlib.sha256()
with open(path, "rb") as fh:
    for chunk in iter(lambda: fh.read(8 << 20), b""):
        m.update(chunk)
        if want_sha:
            sh.update(chunk)
print(m.hexdigest() + ("\t" + sh.hexdigest() if want_sha else ""))
PYHASH
)
  row="$(basename "$f")\t$bytes\t$hashes"
  m=$(meta_for "$f")
  k=$(kind_of "$f")
  if [ -n "$m" ]; then desc=$(read_meta "$m" "$k" "$(basename "$f")")
  else desc=$(printf -- '-\t-\t-\t-\t-\t-\t-\t-\t-'); fi
  printf "$row\t%s\t%s\t%s\n" "$k" "$desc" "$f" >> "$tmp"
  n=$((n + 1))
done < <(find "${DIRS[@]}" "${PRUNE[@]}" -type f -print0 2>/dev/null | sort -z)

if [ "$OUT" = "-" ]; then cat "$tmp"; else cp "$tmp" "$OUT"; log "wrote $OUT ($n file(s))"; fi
