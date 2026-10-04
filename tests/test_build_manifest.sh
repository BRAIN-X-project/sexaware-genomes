#!/usr/bin/env bash
# Behavioural suite for scripts/build_manifest.sh.
#
# Same conventions as tests/test_make_sexaware_gtf.sh: plain bash asserts, no
# pytest, no downloads, about a second.
#
#   bash tests/test_build_manifest.sh
#
# Worth testing because the manifest is itself a published artefact, and both of
# its interesting behaviours fail SILENTLY: a wrong stem attributes a file to
# the wrong assembly or release (a plausible-looking row that is simply false),
# and a missed prune quietly adds tens of GB of build scratch to a release
# description. Neither shows up as an error.
set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
BM="$ROOT/scripts/build_manifest.sh"

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

PASS=0; FAIL=0
ok()  { PASS=$((PASS+1)); printf '  \033[32mok\033[0m   %s\n' "$1"; }
bad() { FAIL=$((FAIL+1)); printf '  \033[31mFAIL\033[0m %s\n' "$1"; [ $# -gt 1 ] && printf '       %s\n' "$2"; }
is()  { if [ "$2" = "$3" ]; then ok "$1"; else bad "$1" "expected '$3', got '$2'"; fi; }

# Source the script's own helpers rather than reimplementing them: a test that
# reimplements stem_of would pass while the real one is broken.
# Extracted one definition at a time: a single sed with three ranges does NOT
# work here, because /^ASSEMBLIES=/,/^}/ runs from that one-line assignment all
# the way to the next closing brace and swallows a function body with it.
cd "$ROOT"   # ASSEMBLIES reads config/assemblies.tsv relative to the repo root
# shellcheck disable=SC1090
eval "$(sed -n '/^ASSEMBLIES=/p' "$BM")"
eval "$(sed -n '/^stem_of()/,/^}/p' "$BM")"
eval "$(sed -n '/^kind_of()/,/^}/p' "$BM")"
[ -n "${ASSEMBLIES:-}" ] || { echo "could not extract ASSEMBLIES from $BM" >&2; exit 1; }
type stem_of >/dev/null 2>&1 || { echo "could not extract stem_of from $BM" >&2; exit 1; }
type kind_of >/dev/null 2>&1 || { echo "could not extract kind_of from $BM" >&2; exit 1; }

echo
echo "stem_of: every published filename shape maps to its product"
stem() { stem_of "$1" 2>/dev/null || echo "(none)"; }
is "genome FASTA"              "$(stem GRCh38.XX.fa.gz)"                            "GRCh38.XX"
is "its md5 sidecar"           "$(stem GRCh38.XX.fa.md5)"                           "GRCh38.XX"
is "chrom.sizes"               "$(stem GRCh38.XY.chrom.sizes)"                      "GRCh38.XY"
is "masked.bed"                "$(stem GRCh38.XY.masked.bed)"                       "GRCh38.XY"
is "genome metadata"           "$(stem GRCh38.XY.metadata.json)"                    "GRCh38.XY"
is "annotation"                "$(stem GRCh38.XX.ensembl-116.gtf.gz)"               "GRCh38.XX.ensembl-116"
is "annotation md5 sidecar"    "$(stem GRCh38.XX.ensembl-116.gtf.gz.md5)"           "GRCh38.XX.ensembl-116"
is "t2g"                       "$(stem GRCh38.XX.ensembl-98.t2g.tsv)"               "GRCh38.XX.ensembl-98"
is "dropped_genes"             "$(stem GRCh38.XX.ensembl-98.dropped_genes.tsv)"     "GRCh38.XX.ensembl-98"
is "noYPAR variant"            "$(stem GRCh38.XY.noYPAR.ensembl-110.gtf.gz)"        "GRCh38.XY.noYPAR.ensembl-110"
is "noYPAR sidecar"            "$(stem GRCh38.XY.noYPAR.ensembl-110.gtf.gz.md5)"    "GRCh38.XY.noYPAR.ensembl-110"
# The k token must NOT be mistaken for a <source>-<release> token, or 80
# mappability rows would claim a release they have nothing to do with.
is "mappability keeps the genome stem" "$(stem GRCh38.XX.k50.multi_read.bw)"        "GRCh38.XX"
is "k150, single_read"         "$(stem GRCm39.XY.k150.single_read.bw)"              "GRCm39.XY"
# An assembly key that itself contains a hyphen and digits.
is "hyphenated assembly key"   "$(stem T2T-CHM13v2.XY.fa.gz)"                       "T2T-CHM13v2.XY"
is "full path, not basename"   "$(stem map/GRCh38/GRCh38.XX/k24/GRCh38.XX.k24.multi_read.bw)" "GRCh38.XX"

echo
echo "stem_of: refuses to guess"
is "unregistered assembly"     "$(stem GRCh99.XX.fa.gz)"        "(none)"
is "no complement token"       "$(stem GRCh38.fa.gz)"           "(none)"
is "unknown complement"        "$(stem GRCh38.XZ.fa.gz)"        "(none)"
is "a source table"            "$(stem chromInfo.txt.gz)"       "(none)"
is "a doc"                     "$(stem README.md)"              "(none)"

echo
echo "kind_of: every published extension is classified, nothing is 'other'"
is "gtf.gz"          "$(kind_of a.gtf.gz)"               annotation
is "fa.gz"           "$(kind_of a.fa.gz)"                genome
is "bw"              "$(kind_of a.bw)"                   mappability
is "chrom.sizes"     "$(kind_of a.chrom.sizes)"          chrom_sizes
is "masked.bed"      "$(kind_of a.masked.bed)"           mask_bed
is "t2g.tsv"         "$(kind_of a.t2g.tsv)"              t2g
is "dropped_genes"   "$(kind_of a.dropped_genes.tsv)"    dropped_genes
is "metadata.json"   "$(kind_of a.metadata.json)"        provenance
is "md5"             "$(kind_of a.fa.md5)"               checksum_sidecar
is "genuinely other" "$(kind_of a.sam)"                  other

echo
echo "end to end on the test fixtures"
out=$("$BM" -d tests/fixtures --no-sha256 -o - 2>/dev/null)
is "header starts with filename" "$(echo "$out" | head -1 | cut -f1)" "filename"
is "no sha256 column with --no-sha256" "$(echo "$out" | head -1 | grep -c sha256)" "0"
is "one row per fixture file" "$(echo "$out" | tail -n +2 | wc -l)" "$(find tests/fixtures -type f | grep -vcE '\.(fai|log)$')"
# md5 column must be the md5 of the bytes on disk -- the whole point of the file.
f=tests/fixtures/mini/mini.gtf
is "md5 matches md5sum" \
  "$(echo "$out" | awk -F'\t' -v p="$f" '$NF==p {print $3}')" \
  "$(md5sum "$f" | awk '{print $1}')"
is "bytes matches stat" \
  "$(echo "$out" | awk -F'\t' -v p="$f" '$NF==p {print $2}')" \
  "$(stat -c %s "$f")"

echo
echo "pruning build scratch"
mkdir -p "$WORK/t/GRCh38.XX_genmap_index" "$WORK/t"
head -c 2048 /dev/urandom > "$WORK/t/GRCh38.XX_genmap_index/index.lf.drv"
printf 'X\t100\n' > "$WORK/t/GRCh38.XX.chrom.sizes"
n_default=$("$BM" -d "$WORK/t" --no-sha256 -o - 2>/dev/null | tail -n +2 | wc -l)
n_all=$("$BM" -d "$WORK/t" --no-sha256 --include-intermediates -o - 2>/dev/null | tail -n +2 | wc -l)
is "index dir pruned by default"        "$n_default" "1"
is "--include-intermediates keeps it"   "$n_all"     "2"
is "and it would have been kind=other" \
  "$("$BM" -d "$WORK/t" --no-sha256 --include-intermediates -o - 2>/dev/null \
     | awk -F'\t' 'NR>1 && $4=="other"' | wc -l)" "1"

echo
if [ "$FAIL" -eq 0 ]; then printf '\033[32mpassed %d, failed 0\033[0m\n' "$PASS"; else
  printf '\033[31mpassed %d, failed %d\033[0m\n' "$PASS" "$FAIL"; fi
[ "$FAIL" -eq 0 ]
