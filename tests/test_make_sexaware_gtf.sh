#!/usr/bin/env bash
# Behavioural suite for scripts/make_sexaware_gtf.py.
#
# Plain bash asserts, no pytest: the repo's other scripts are bash, and this
# must be runnable in the GenMap env without adding a test dependency. Uses only
# the hand-written tests/fixtures/mini assembly, so it needs no download and
# finishes in about a second.
#
#   bash tests/test_make_sexaware_gtf.sh
#
# Every subtle behaviour of the filter is covered here, which is what makes it
# safe to change: the GRCh37 shared-gene_id trap in particular is a silent
# corruption that no amount of eyeballing a 1.4 GB GTF would catch.
set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
FIX="$HERE/fixtures/mini"
GTF="$FIX/mini.gtf"
MK="$ROOT/scripts/make_sexaware_gtf.py"
PY=${PYTHON:-python3}

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
cp "$FIX"/MINI.* "$WORK/"

PASS=0; FAIL=0
ok()   { PASS=$((PASS+1)); printf '  \033[32mok\033[0m   %s\n' "$1"; }
bad()  { FAIL=$((FAIL+1)); printf '  \033[31mFAIL\033[0m %s\n' "$1"; [ $# -gt 1 ] && printf '       %s\n' "$2"; }
is()   { if [ "$2" = "$3" ]; then ok "$1"; else bad "$1" "expected '$3', got '$2'"; fi; }
isnt() { if [ "$2" != "$3" ]; then ok "$1"; else bad "$1" "expected anything but '$3'"; fi; }
gt()   { if [ "$2" -gt "$3" ] 2>/dev/null; then ok "$1"; else bad "$1" "expected > $3, got '$2'"; fi; }
empty(){ if [ ! -s "$2" ]; then ok "$1"; else bad "$1" "$(head -3 "$2")"; fi; }

# Run the filter into a clean subdirectory, echoing its exit code.
run() {
  local out=$1; shift
  rm -rf "$WORK/$out"; mkdir -p "$WORK/$out"; cp "$FIX"/MINI.* "$WORK/$out/"
  $PY "$MK" --gtf "$GTF" --x-contig X --y-contig Y --prefix MINI \
    --annotation-source mini --annotation-release 1 \
    --outdir "$WORK/$out" --force --quiet "$@" >"$WORK/$out.log" 2>&1
  echo $?
}
meta() { $PY -c "import json,sys;print(json.load(open(sys.argv[1]))[sys.argv[2]])" "$1" "$2"; }
feat() { grep -vc '^#' "$1"; }

echo "== XX: the whole Y contig is removed =="
rc=$(run xx --variant XX --on-unknown-seqname warn --on-out-of-range warn)
is  "exit 0" "$rc" 0
XX=$WORK/xx/MINI.XX.mini-1.gtf.gz
is  "1. zero features remain on Y" "$(zcat "$XX" | grep -v '^#' | awk -F'\t' '$1=="Y"' | wc -l)" 0
# 2. Equivalence in one command: the output IS the input minus Y, byte for byte.
diff <(grep -v '^#' "$GTF" | awk -F'\t' '$1!="Y"') \
     <(zcat "$XX" | grep -v '^#') > "$WORK/xx.diff" 2>&1
empty "2. output == input minus Y, exactly" "$WORK/xx.diff"
is  "   features 40 -> 20" "$(feat <(zcat "$XX"))" 20
is  "   upstream #!genome-build header preserved" \
    "$(zcat "$XX" | grep -c '^#!genome-build MINIv1')" 1
is  "   20 Y features recorded as dropped" \
    "$(meta "$WORK/xx/MINI.XX.mini-1.metadata.json" features_dropped)" 20
is  "   dropped_genes.tsv lists the 5 Y genes" \
    "$(tail -n +2 "$WORK/xx/MINI.XX.mini-1.dropped_genes.tsv" | wc -l)" 5
is  "   reason is whole_contig_masked" \
    "$(awk -F'\t' 'NR>1 {print $11}' "$WORK/xx/MINI.XX.mini-1.dropped_genes.tsv" | sort -u)" \
    "whole_contig_masked"
is  "   the un-named Y gene still appears, with an empty gene_name" \
    "$(awk -F'\t' '$1=="ENSG_NONAME" {print $2"|"$3}' "$WORK/xx/MINI.XX.mini-1.dropped_genes.tsv")" "|lncRNA"

echo "== XY-noYPAR: only the Y-side PAR goes =="
rc=$(run np --variant XY-noYPAR --on-unknown-seqname warn --on-out-of-range warn)
is  "exit 0 (partial overlap only warns by default)" "$rc" 0
NP=$WORK/np/MINI.XY.noYPAR.mini-1.gtf.gz
dropped() { cut -f1 "$WORK/np/MINI.XY.noYPAR.mini-1.dropped_genes.tsv" | tail -n +2 | sort -u | tr '\n' ' '; }
is  "3. PAR1/PAR2 genes dropped, NPY gene kept" \
    "$(zcat "$NP" | grep -v '^#' | awk -F'\t' '$1=="Y" && $3=="gene" {print $9}' | sed -n 's/^gene_id "\([^"]*\)".*/\1/p' | sort | tr '\n' ' ')" \
    "ENSG_NPY1 "
is  "   drop set is exactly the four PAR/PAB genes" "$(dropped)" \
    "ENSG_NONAME ENSG_PAB ENSG_PAR2G ENSG_SHARED "
is  "4. the PAB-straddling gene is flagged partial" \
    "$(awk -F'\t' '$1=="ENSG_PAB" {print $11}' "$WORK/np/MINI.XY.noYPAR.mini-1.dropped_genes.tsv")" \
    "y_par_overlap_partial"
is  "   genes_partially_overlapping == 1" \
    "$(meta "$WORK/np/MINI.XY.noYPAR.mini-1.metadata.json" genes_partially_overlapping)" 1
is  "   zero survivors overlap the mask" \
    "$(meta "$WORK/np/MINI.XY.noYPAR.mini-1.metadata.json" features_on_masked_sequence_kept)" 0

# 5. THE GRCh37 REGRESSION. ENSG_SHARED is on both X and Y; the Y copy must go
#    and the X copy must stay. A drop set keyed on gene_id alone passes every
#    other assertion in this file and fails only this one.
is  "5. the X copy of the shared gene_id survives" \
    "$(zcat "$NP" | grep -v '^#' | awk -F'\t' '$1=="X"' | wc -l)" 4
diff <(awk -F'\t' '$1=="X"' "$GTF") <(zcat "$NP" | grep -v '^#' | awk -F'\t' '$1=="X"') \
    > "$WORK/x.diff" 2>&1
empty "   chrX is bit-for-bit untouched" "$WORK/x.diff"
is  "   and the Y copy is gone" \
    "$(zcat "$NP" | grep -v '^#' | awk -F'\t' '$1=="Y" && /ENSG_SHARED/' | wc -l)" 0

echo "== --on-partial-overlap fail =="
rc=$(run npf --variant XY-noYPAR --on-partial-overlap fail \
         --on-unknown-seqname warn --on-out-of-range warn)
is  "4b. exit 3 when a gene straddles the PAB" "$rc" 3

echo "== --overlap-rule contained =="
rc=$(run cont --variant XY-noYPAR --overlap-rule contained \
         --on-unknown-seqname warn --on-out-of-range warn)
is  "exit 0" "$rc" 0
is  "6. the straddling gene is kept" \
    "$(zcat "$WORK/cont/MINI.XY.noYPAR.mini-1.gtf.gz" | grep -c 'ENSG_PAB')" 4
gt  "   and is reported as sitting on masked sequence" \
    "$(meta "$WORK/cont/MINI.XY.noYPAR.mini-1.metadata.json" features_on_masked_sequence_kept)" 0

echo "== --on-unknown-seqname =="
rc=$(run unkf --variant XY --on-out-of-range warn)
is  "7a. fail (the default) -> exit 2" "$rc" 2
grep -q 'HSCHR1_1_CTG3' "$WORK/unkf.log" && ok "    and names the offending contig" \
  || bad "    and names the offending contig" "$(tail -2 "$WORK/unkf.log")"
rc=$(run unkw --variant XY --on-unknown-seqname warn --on-out-of-range warn)
is  "7b. warn -> exit 0, features kept" "$rc" 0
is  "    features_out == features_in" \
    "$(meta "$WORK/unkw/MINI.XY.mini-1.metadata.json" features_out)" 40
rc=$(run unkd --variant XY --on-unknown-seqname drop --on-out-of-range warn)
is  "7c. drop -> exit 0, the 4 patch features removed" "$rc" 0
is  "    features_out == 36" \
    "$(meta "$WORK/unkd/MINI.XY.mini-1.metadata.json" features_out)" 36
is  "    reason is seqname_not_in_genome" \
    "$(awk -F'\t' 'NR>1 {print $11}' "$WORK/unkd/MINI.XY.mini-1.dropped_genes.tsv" | sort -u)" \
    "seqname_not_in_genome"
# Restricting to the genome is NOT sex-aware filtering. Ensembl r75 ships only
# the all-scaffolds GTF, so a real XY build drops ~5,772 genes on MHC
# haplotypes and patches; counting those as sex-aware drops would tell a user
# their standard XY annotation had been filtered when it had not.
is  "    counted as off-genome, not sex-aware" \
    "$(meta "$WORK/unkd/MINI.XY.mini-1.metadata.json" genes_dropped_off_genome)" 1
is  "    sex-aware drop count is 0 for standard XY" \
    "$(meta "$WORK/unkd/MINI.XY.mini-1.metadata.json" genes_dropped_sexaware)" 0

rc=$(run unks --variant XY-noYPAR --on-unknown-seqname drop --on-out-of-range warn)
is  "7d. XY-noYPAR separates the two populations" "$rc" 0
is  "    off-genome == 1 (the patch gene)" \
    "$(meta "$WORK/unks/MINI.XY.noYPAR.mini-1.metadata.json" genes_dropped_off_genome)" 1
is  "    sex-aware == 4 (PAR1 x2, PAR2, PAB-straddler)" \
    "$(meta "$WORK/unks/MINI.XY.noYPAR.mini-1.metadata.json" genes_dropped_sexaware)" 4
is  "    and the two sum to genes_dropped" \
    "$(meta "$WORK/unks/MINI.XY.noYPAR.mini-1.metadata.json" genes_dropped)" 5

echo "== --on-out-of-range =="
rc=$(run oor --variant XY --on-unknown-seqname warn)
is  "8. a feature past the end of its contig -> exit 2" "$rc" 2
grep -q 'max_end' "$WORK/oor.log" && ok "   and reports the contig and the overrun" \
  || bad "   and reports the contig and the overrun" "$(tail -2 "$WORK/oor.log")"

echo "== XY: standard annotation passes through =="
XY=$WORK/unkw/MINI.XY.mini-1.gtf.gz
is  "9. features_in == features_out" "$(feat <(zcat "$XY"))" 40
diff <(grep -v '^#' "$GTF") <(zcat "$XY" | grep -v '^#') > "$WORK/xy.diff" 2>&1
empty "   every feature line is identical to upstream" "$WORK/xy.diff"
is  "   provenance header records the variant" \
    "$(zcat "$XY" | grep -c '^#!sexaware-variant XY$')" 1
is  "   and binds the GTF to one exact genome build" \
    "$(zcat "$XY" | awk '/^#!sexaware-paired-genome-md5/ {print $2}')" \
    "e5b7e9980fd4b2a6a0a0d1b5c9f0e3a1"

echo "== checksums and reproducibility =="
md5file=$WORK/xx/MINI.XX.mini-1.gtf.gz.md5
is  "10. the sidecar names the file it hashes" \
    "$(awk '{print $2}' "$md5file")" "MINI.XX.mini-1.gtf.gz"
( cd "$WORK/xx" && md5sum -c --status MINI.XX.mini-1.gtf.gz.md5 ) \
  && ok "    md5sum -c verifies" || bad "    md5sum -c verifies"
is  "    and matches metadata.json" \
    "$(awk '{print $1}' "$md5file")" \
    "$(meta "$WORK/xx/MINI.XX.mini-1.metadata.json" output_gtf_md5)"

# 11. Determinism. gzip stamps an mtime by default, so without mtime=0 these two
#     runs differ in bytes and the published md5 could never be reproduced.
sleep 1.1
rc=$(run xx2 --variant XX --on-unknown-seqname warn --on-out-of-range warn)
is  "11. second run exits 0" "$rc" 0
if cmp -s "$XX" "$WORK/xx2/MINI.XX.mini-1.gtf.gz"; then
  ok "    two runs are byte-identical"
else
  bad "    two runs are byte-identical" "gzip header is carrying a timestamp"
fi

echo "== --variant all and --emit-t2g =="
rc=$(run all --variant all --emit-t2g --on-unknown-seqname warn --on-out-of-range warn)
is  "exit 0" "$rc" 0
is  "writes all three variants" \
    "$(ls "$WORK/all"/MINI.*.mini-1.gtf.gz | wc -l)" 3
is  "t2g maps every surviving transcript of XY" \
    "$(tail -n +2 "$WORK/all/MINI.XY.mini-1.t2g.tsv" | wc -l)" 8
# 8 and not 9: the fixture gives the X and Y copies of ENSG_SHARED the same
# transcript_id, so the map collapses them. Real Ensembl releases give the PAR
# copies distinct ids (measured: 0 shared transcript_ids on X/Y in GRCh37 75 and
# 87), but a silent collapse would mis-assign a transcript's reads, so it warns.
# Reports ENST_AUTO1 (two gene_ids -> genuinely ambiguous) and NOT ENST_SHARED
# (one gene_id on two contigs -> collapsing is harmless). Distinguishing the two
# is the whole point; warning about both would train people to ignore it.
is  "reports the AMBIGUOUS transcript_id" \
    "$(meta "$WORK/all/MINI.XY.mini-1.metadata.json" t2g_collisions)" "['ENST_AUTO1']"
is  "t2g for XX omits the Y transcripts" \
    "$(tail -n +2 "$WORK/all/MINI.XX.mini-1.t2g.tsv" | wc -l)" 4
is  "and still contains the X copy of the shared gene_id" \
    "$(awk -F'\t' '$2=="ENSG_SHARED"' "$WORK/all/MINI.XX.mini-1.t2g.tsv" | wc -l)" 1

echo
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ] || exit 1
