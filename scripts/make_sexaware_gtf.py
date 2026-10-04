#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build a sex-chromosome-complement (SCC) informed annotation GTF.

Companion to `mask_genome.py`: that script hard-masks the sequence, this one
filters the annotation so the two agree. The single rule everything follows is

    a sex-aware GTF contains NO feature that overlaps N-masked sequence in the
    corresponding sex-aware FASTA

which is why the preferred filter input is the `<prefix>.<XX|XY>.masked.bed`
that `mask_genome.py` already writes, rather than `config/par/*.bed`: the
annotation is filtered against the mask that was actually applied, so GTF and
FASTA cannot silently diverge, and any `--extra-mask` BED is honoured for free.

Three variants are produced:

  XX          the entire Y contig is removed. Matches `*.XX.fa`, where all of Y
              is N.
  XY          standard upstream annotation, passed through unchanged (plus the
              provenance header). Matches `*.XY.fa`, where the Y PAR is N, so
              the Y-PAR genes are annotated on masked sequence: they will
              produce permanent zero counts and duplicate the X copy's symbol.
              Shipped because many pipelines require the exact upstream GTF.
  XY-noYPAR   as XY, but the features on the Y-side PAR are dropped, which is
              the variant that is actually consistent with `*.XY.fa`.

WARNING -- XX and XY do not share a gene universe. featureCounts/HTSeq on XX
samples emit no chrY rows at all, so naively cbind()-ing an XX and an XY matrix
silently drops genes or fills them with NA. Use the XY gene set as the canonical
index and zero-fill chrY for XX samples; `<out>.dropped_genes.tsv` is exactly
the list you need to do that.

Outputs (under --outdir):
  <prefix>.<variant>.<source>-<release>.gtf.gz           filtered annotation
  <prefix>.<variant>.<source>-<release>.gtf.gz.md5       md5 of those bytes
  <prefix>.<variant>.<source>-<release>.metadata.json    provenance + counts
  <prefix>.<variant>.<source>-<release>.dropped_genes.tsv  what was removed, why
  <prefix>.<variant>.<source>-<release>.t2g.tsv          transcript->gene map

Exit codes: 2 = annotation/genome mismatch (unknown seqname or out-of-range
coordinate), 3 = partial overlap under --on-partial-overlap fail.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mask_genome import (  # noqa: E402  (deliberate: single source of truth)
    contig_matches,
    git_describe,
    is_y_side,
    md5_of_file,
    read_assembly_registry,
    read_bed_intervals,
    strip_chr,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
TOOL_VERSION = "0.2.0"

VARIANTS = ("XX", "XY", "XY-noYPAR")
# Filename token per variant, and which FASTA complement each one pairs with.
VARIANT_TOKEN = {"XX": "XX", "XY": "XY", "XY-noYPAR": "XY.noYPAR"}
VARIANT_COMPLEMENT = {"XX": "XX", "XY": "XY", "XY-noYPAR": "XY"}

# Per-seqname decisions in the hot loop.
KEEP, DROP_ALL, INSPECT = 0, 1, 2

_T0 = time.time()
VERBOSE = True


def log(msg: str) -> None:
    """Print a timestamped progress line to stderr when VERBOSE is on."""
    if VERBOSE:
        print(f"[make_sexaware_gtf {time.strftime('%H:%M:%S')} +{time.time() - _T0:5.1f}s] {msg}",
              file=sys.stderr, flush=True)


# --------------------------------------------------------------------------- #
# coordinates
# --------------------------------------------------------------------------- #
def overlaps_bed(gtf_start: int, gtf_end: int, bed_start: int, bed_end: int) -> bool:
    """True if a 1-based inclusive GTF interval meets a 0-based half-open BED one.

    GTF is 1-based inclusive, BED is 0-based half-open, so the GTF start is
    converted with `gtf_start - 1` and the GTF end needs no adjustment.
    """
    return (gtf_start - 1) < bed_end and gtf_end > bed_start


def contained_in_bed(gtf_start: int, gtf_end: int, bed_start: int, bed_end: int) -> bool:
    """True if a 1-based inclusive GTF interval lies wholly inside a BED one."""
    return (gtf_start - 1) >= bed_start and gtf_end <= bed_end


# --------------------------------------------------------------------------- #
# GTF parsing helpers (bytes; the attribute column may hold non-ASCII)
# --------------------------------------------------------------------------- #
def attr(attrs: bytes, key: bytes) -> bytes | None:
    """Pull one quoted GTF attribute value out of column 9, or None.

    Matches only at the start of the column or directly after '; ', so looking
    up `gene_id` does not accidentally match `ref_gene_id`.
    """
    needle = key + b' "'
    i = 0
    while True:
        i = attrs.find(needle, i)
        if i == -1:
            return None
        if i == 0 or attrs[i - 2:i] == b"; ":
            break
        i += 1
    i += len(needle)
    j = attrs.find(b'"', i)
    return attrs[i:j] if j != -1 else None


def first_attr(attrs: bytes, keys: tuple) -> bytes:
    """First present attribute among `keys` (Ensembl and GENCODE spellings)."""
    for key in keys:
        value = attr(attrs, key)
        if value is not None:
            return value
    return b""


GENE_ID = (b"gene_id",)
TX_ID = (b"transcript_id",)
GENE_NAME = (b"gene_name", b"gene")
GENE_BIOTYPE = (b"gene_biotype", b"gene_type")


def open_gtf_binary(path: Path):
    """Open a plain or gzipped GTF for streaming binary reads."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def open_gtf_out(path: Path):
    """Open a GTF output for binary writes, gzip written deterministically.

    mtime=0 and an empty embedded filename keep the bytes reproducible, so the
    md5 sidecar is a real identity for the content rather than a per-run nonce.
    """
    if str(path).endswith(".gz"):
        raw = open(path, "wb")
        handle = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
        handle.myfileobj = raw
        return handle
    return open(path, "wb")


# --------------------------------------------------------------------------- #
# registries
# --------------------------------------------------------------------------- #
def read_annotation_registry(path: Path) -> list:
    """Parse config/annotations.tsv, skipping comment lines, as a list of rows."""
    with open(path, "r", encoding="utf-8") as handle:
        rows = [line for line in handle if not line.startswith("#") and line.strip()]
    return list(csv.DictReader(rows, delimiter="\t"))


def read_chrom_sizes(path: Path) -> dict:
    """Read a two-column chrom.sizes into {contig_token: length}."""
    sizes = {}
    with open(path, "r", encoding="ascii") as handle:
        for line in handle:
            if not line.strip():
                continue
            name, length = line.split()[:2]
            sizes[name] = int(length)
    return sizes


# --------------------------------------------------------------------------- #
# mask resolution
# --------------------------------------------------------------------------- #
def resolve_mask_intervals(
    variant: str,
    masked_bed: Path | None,
    par_bed: Path | None,
    y_contig: str,
    genome_sizes: dict,
    extra_drop_beds: list,
) -> tuple[dict, list]:
    """Return ({bare_contig: [(start, end, name)]}, provenance_notes).

    Preferred input is the masked.bed that mask_genome.py wrote for the paired
    FASTA, so the annotation is filtered against the mask that was really
    applied. Falling back to a PAR BED reconstructs the same intervals, but
    cannot know about --extra-mask regions.
    """
    mask: dict = {}
    notes: list = []

    def add(bare: str, start: int, end: int, name: str):
        mask.setdefault(bare, []).append((start, end, name))

    if variant == "XY":
        notes.append("no mask: standard upstream annotation")
    elif masked_bed is not None:
        for chrom, start, end, name in read_bed_intervals(masked_bed):
            add(strip_chr(chrom), start, end, name)
        notes.append(f"mask from {masked_bed.name} (the mask applied to the paired FASTA)")
    elif variant == "XX":
        bare_y = strip_chr(y_contig)
        length = next((size for token, size in genome_sizes.items()
                       if contig_matches(token, y_contig)), None)
        if length is None:
            raise SystemExit(f"ERROR: Y contig '{y_contig}' is absent from the chrom.sizes; "
                             "cannot determine the whole-Y interval.")
        add(bare_y, 0, length, "whole_Y_mask")
        notes.append("mask derived from chrom.sizes (whole Y); no masked.bed available")
    else:  # XY-noYPAR
        if par_bed is None:
            raise SystemExit("ERROR: XY-noYPAR needs either --masked-bed (preferred) or --par-bed.")
        n = 0
        for chrom, start, end, name in read_bed_intervals(par_bed):
            if is_y_side(chrom, y_contig):
                add(strip_chr(y_contig), start, end, f"{name}_masked")
                n += 1
        if n == 0:
            raise SystemExit(f"ERROR: {par_bed} has no Y-side interval for y_contig={y_contig}.")
        notes.append(f"mask reconstructed from {par_bed.name} (Y-side PAR rows only); "
                     "--extra-mask regions are NOT represented")

    for bed in extra_drop_beds:
        for chrom, start, end, name in read_bed_intervals(Path(bed)):
            add(strip_chr(chrom), start, end, f"extra:{name}")
        notes.append(f"extra drop intervals from {Path(bed).name}")

    for bare in mask:
        mask[bare].sort()
    return mask, notes


# --------------------------------------------------------------------------- #
# pass 1: scan
# --------------------------------------------------------------------------- #
class GeneRecord:
    """Accumulated evidence about one (seqname, gene_id) on a filtered contig."""

    __slots__ = ("name", "biotype", "strand", "start", "end",
                 "transcripts", "exons", "n_features", "hit", "partial", "reason")

    def __init__(self):
        self.name = ""
        self.biotype = ""
        self.strand = "."
        self.start = None
        self.end = None
        self.transcripts = set()
        self.exons = []
        self.n_features = 0
        self.hit = False
        self.partial = False
        self.reason = ""

    def exonic_bp(self) -> int:
        """Merged exonic span, so overlapping transcript exons are not counted twice."""
        total = 0
        cur_s = cur_e = None
        for start, end in sorted(self.exons):
            if cur_e is None or start > cur_e:
                if cur_e is not None:
                    total += cur_e - cur_s + 1
                cur_s, cur_e = start, end
            else:
                cur_e = max(cur_e, end)
        if cur_e is not None:
            total += cur_e - cur_s + 1
        return total


def scan_gtf(
    gtf: Path,
    mask: dict,
    genome_sizes: dict,
    on_unknown: str,
    overlap_rule: str,
    drop_unit: str,
    y_contig: str,
) -> dict:
    """One pass over the GTF: header, counts, ranges, and the drop set.

    Returns everything the write pass needs, so the provenance header can carry
    final counts and the write pass can be verified against a prediction made
    before a single byte was written.
    """
    header: list = []
    features_in = 0
    per_seq_in: dict = {}
    max_end: dict = {}
    genes: dict = {}          # (raw_seq, gene_id) -> GeneRecord
    tx_hit: set = set()       # (raw_seq, transcript_id)
    tx_of_gene: dict = {}     # (raw_seq, gene_id) -> {transcript_id}
    first_pass_seqnames: dict = {}
    malformed = 0

    # Seqname decisions are discovered lazily, but need the full seqname set to
    # be classified consistently, so classification happens per newly seen name.
    genome_bare = {strip_chr(t): (t, s) for t, s in genome_sizes.items()}
    decision: dict = {}
    unknown_seen: dict = {}

    def classify(raw: bytes) -> int:
        bare = strip_chr(raw.decode("ascii", "replace"))
        if bare not in genome_bare:
            unknown_seen[raw] = 0
            return DROP_ALL if on_unknown == "drop" else KEEP
        intervals = mask.get(bare)
        if not intervals:
            return KEEP
        _token, size = genome_bare[bare]
        covered_start = min(s for s, _e, _n in intervals)
        covered_end = max(e for _s, e, _n in intervals)
        contiguous = all(intervals[i][0] <= intervals[i - 1][1] for i in range(1, len(intervals)))
        if covered_start <= 0 and covered_end >= size and contiguous:
            return DROP_ALL
        return INSPECT

    with open_gtf_binary(gtf) as handle:
        for line in handle:
            if line.startswith(b"#"):
                if features_in == 0:
                    header.append(line)
                continue
            tab = line.find(b"\t")
            if tab == -1:
                malformed += 1
                continue
            raw = line[:tab]
            features_in += 1
            per_seq_in[raw] = per_seq_in.get(raw, 0) + 1

            act = decision.get(raw)
            if act is None:
                act = classify(raw)
                decision[raw] = act
                first_pass_seqnames[raw] = act

            # maxsplit=8 leaves the attribute column unsplit and, more to the
            # point, lets a KEEP line be accounted for without ever touching it.
            fields = line.split(b"\t", 8)
            if len(fields) < 9:
                malformed += 1
                continue
            start, end = int(fields[3]), int(fields[4])
            if end > max_end.get(raw, 0):
                max_end[raw] = end

            if act == KEEP:
                continue

            # Only contigs that are dropped or inspected pay the attribute parse.
            feature, strand, attrs = fields[2], fields[6], fields[8]
            gene_id = first_attr(attrs, GENE_ID)
            key = (raw, gene_id)
            rec = genes.get(key)
            if rec is None:
                rec = genes[key] = GeneRecord()
            rec.n_features += 1
            if rec.start is None or start < rec.start:
                rec.start = start
            if rec.end is None or end > rec.end:
                rec.end = end
            if not rec.name:
                name = first_attr(attrs, GENE_NAME)
                if name:
                    rec.name = name.decode("utf-8", "replace")
            if not rec.biotype:
                biotype = first_attr(attrs, GENE_BIOTYPE)
                if biotype:
                    rec.biotype = biotype.decode("utf-8", "replace")
            if rec.strand == ".":
                rec.strand = strand.decode("ascii", "replace")
            tx = attr(attrs, b"transcript_id")
            if tx:
                rec.transcripts.add(tx)
                tx_of_gene.setdefault(key, set()).add(tx)
            if feature == b"exon":
                rec.exons.append((start, end))

            if act == DROP_ALL:
                if raw in unknown_seen:
                    rec.hit, rec.reason = True, "seqname_not_in_genome"
                else:
                    rec.hit, rec.reason = True, "whole_contig_masked"
                continue

            # INSPECT: test this feature against the mask intervals.
            bare = strip_chr(raw.decode("ascii", "replace"))
            hit = partial = False
            on_y = False
            for m_start, m_end, m_name in mask[bare]:
                if not overlaps_bed(start, end, m_start, m_end):
                    continue
                inside = contained_in_bed(start, end, m_start, m_end)
                if overlap_rule == "contained" and not inside:
                    partial = True
                    continue
                hit = True
                if not inside:
                    partial = True
                if not m_name.startswith("extra:"):
                    on_y = True
                break
            if not hit:
                if partial:
                    rec.partial = True
                continue

            reason_y = "y_par_overlap_partial" if partial else "y_par_overlap"
            reason = reason_y if (on_y and is_y_side(bare, y_contig)) else "extra_mask_overlap"
            rec.hit = True
            if partial:
                rec.partial = True
            if not rec.reason or rec.reason == "y_par_overlap":
                rec.reason = reason
            if tx:
                tx_hit.add((raw, tx))
            elif drop_unit == "feature":
                rec.reason = reason

    # ---- build the drop set --------------------------------------------- #
    # Keyed by (seqname, id) and NEVER by the id alone. Whether a gene that has
    # an X and a Y PAR copy gets one shared gene_id or two distinct ones is a
    # per-release convention, not a guarantee: Ensembl GRCh38 r116 gives the Y
    # copies their own IDs (PPP2R3B on Y is ENSG00000292327, not the X's
    # ENSG00000167393), but nothing in the format requires that, and the GRCh37
    # REST API still reports both locations for one ID. Dropping by bare id
    # would, under the sharing convention, delete the X copy as well -- the
    # exact opposite of what an XY reference is for, and it would publish
    # looking correct. Keying on the pair is right under either convention, so
    # we never have to ask which one a release follows.
    drop_genes: set = set()
    drop_transcripts: set = set()
    if drop_unit == "gene":
        drop_genes = {key for key, rec in genes.items() if rec.hit}
    elif drop_unit == "transcript":
        drop_transcripts = set(tx_hit)
        for key, rec in genes.items():
            if not rec.hit:
                continue
            txs = tx_of_gene.get(key, set())
            if not txs or all((key[0], tx) in drop_transcripts for tx in txs):
                drop_genes.add(key)
    else:  # feature
        pass

    # Two populations get dropped here for completely different reasons, and
    # reporting one number for both is actively misleading: a reader of
    # "[XY] 5,772 gene(s) dropped" would conclude the XY annotation had been
    # sex-aware filtered, when in fact every one of those genes sits on an MHC
    # haplotype or patch scaffold that the primary-assembly FASTA does not
    # contain (Ensembl r75 ships only the all-scaffolds GTF). Keep them apart
    # everywhere they are reported.
    off_genome = {key for key in drop_genes
                  if genes[key].reason == "seqname_not_in_genome"}

    return {
        "header": header,
        "features_in": features_in,
        "per_seq_in": per_seq_in,
        "max_end": max_end,
        "genes": genes,
        "decision": decision,
        "unknown_seqnames": list(unknown_seen),
        "drop_genes": drop_genes,
        "drop_genes_off_genome": off_genome,
        "drop_genes_sexaware": drop_genes - off_genome,
        "drop_transcripts": drop_transcripts,
        "malformed": malformed,
        "seqname_actions": first_pass_seqnames,
    }


def predict_output(scan: dict, drop_unit: str) -> int:
    """How many feature lines the write pass should emit, computed before writing.

    Comparing this with the real count afterwards is a free internal consistency
    check: if the two disagree the filter is not doing what it just claimed.
    """
    if drop_unit == "feature":
        return -1  # per-line decisions cannot be predicted from aggregates
    dropped = sum(rec.n_features for key, rec in scan["genes"].items()
                  if key in scan["drop_genes"])
    return scan["features_in"] - dropped


# --------------------------------------------------------------------------- #
# pass 2: write
# --------------------------------------------------------------------------- #
def write_gtf(
    gtf: Path,
    out_path: Path,
    scan: dict,
    mask: dict,
    header_lines: list,
    drop_unit: str,
    overlap_rule: str,
    x_contig: str,
    y_contig: str,
    emit_t2g: bool,
    t2g_path: Path | None,
) -> dict:
    """Write the filtered GTF and return the post-condition counters."""
    decision = scan["decision"]
    drop_genes = scan["drop_genes"]
    drop_transcripts = scan["drop_transcripts"]

    features_out = 0
    per_seq_out: dict = {}
    y_features_out = 0
    x_features_out = 0
    kept_on_masked = 0
    t2g: dict = {}
    t2g_collisions: set = set()

    with open_gtf_binary(gtf) as handle, open_gtf_out(out_path) as out:
        for line in header_lines:
            out.write(line)

        for line in handle:
            if line.startswith(b"#"):
                continue
            tab = line.find(b"\t")
            if tab == -1:
                continue
            raw = line[:tab]
            act = decision.get(raw, KEEP)
            if act == DROP_ALL:
                continue

            if act == INSPECT:
                fields = line.split(b"\t", 8)
                if len(fields) < 9:
                    continue
                attrs = fields[8]
                gene_id = first_attr(attrs, GENE_ID)
                if (raw, gene_id) in drop_genes:
                    continue
                tx = attr(attrs, b"transcript_id")
                if tx and (raw, tx) in drop_transcripts:
                    continue
                if drop_unit == "feature":
                    start, end = int(fields[3]), int(fields[4])
                    bare = strip_chr(raw.decode("ascii", "replace"))
                    drop_this = False
                    for m_start, m_end, _m_name in mask.get(bare, ()):
                        if not overlaps_bed(start, end, m_start, m_end):
                            continue
                        if overlap_rule == "contained" and not contained_in_bed(
                                start, end, m_start, m_end):
                            continue
                        drop_this = True
                        break
                    if drop_this:
                        continue
                start, end = int(fields[3]), int(fields[4])
                bare = strip_chr(raw.decode("ascii", "replace"))
                if any(overlaps_bed(start, end, s, e) for s, e, _n in mask.get(bare, ())):
                    kept_on_masked += 1
                if emit_t2g:
                    tx_id = attr(attrs, b"transcript_id")
                    if tx_id:
                        prev = t2g.get(tx_id)
                        if prev is None:
                            t2g[tx_id] = (gene_id, first_attr(attrs, GENE_NAME))
                        elif prev[0] != gene_id:
                            t2g_collisions.add(tx_id)
            elif emit_t2g:
                parts = line.split(b"\t", 8)
                if len(parts) < 9:
                    continue
                attrs = parts[8]
                tx_id = attr(attrs, b"transcript_id")
                if tx_id:
                    gene_id = first_attr(attrs, GENE_ID)
                    prev = t2g.get(tx_id)
                    if prev is None:
                        t2g[tx_id] = (gene_id, first_attr(attrs, GENE_NAME))
                    elif prev[0] != gene_id:
                        t2g_collisions.add(tx_id)

            out.write(line)
            features_out += 1
            per_seq_out[raw] = per_seq_out.get(raw, 0) + 1
            bare = strip_chr(raw.decode("ascii", "replace"))
            if contig_matches(bare, y_contig):
                y_features_out += 1
            elif contig_matches(bare, x_contig):
                x_features_out += 1

    if emit_t2g and t2g_path is not None:
        with open(t2g_path, "w", encoding="utf-8") as handle:
            handle.write("transcript_id\tgene_id\tgene_name\n")
            for tx_id, (gene_id, gene_name) in t2g.items():
                handle.write(f"{tx_id.decode('utf-8', 'replace')}\t"
                             f"{gene_id.decode('utf-8', 'replace')}\t"
                             f"{gene_name.decode('utf-8', 'replace')}\n")

    return {
        "features_out": features_out,
        "per_seq_out": per_seq_out,
        "y_features_out": y_features_out,
        "x_features_out": x_features_out,
        "features_on_masked_sequence_kept": kept_on_masked,
        "n_transcripts": len(t2g),
        "t2g_collisions": sorted(t.decode("utf-8", "replace") for t in t2g_collisions),
    }


def write_dropped_genes(path: Path, scan: dict) -> int:
    """Write the per-gene record of what was removed and why."""
    rows = []
    for (raw, gene_id), rec in scan["genes"].items():
        if not rec.hit:
            continue
        rows.append((
            gene_id.decode("utf-8", "replace"),
            rec.name,
            rec.biotype,
            raw.decode("ascii", "replace"),
            rec.start, rec.end, rec.strand,
            len(rec.transcripts), len(rec.exons), rec.exonic_bp(),
            rec.reason,
        ))
    rows.sort(key=lambda r: (r[3], r[4] or 0, r[0]))
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("gene_id\tgene_name\tgene_biotype\tseqname\tstart\tend\tstrand\t"
                     "n_transcripts\tn_exons\texonic_bp\treason\n")
        for row in rows:
            handle.write("\t".join(str(x) for x in row) + "\n")
    return len(rows)


# --------------------------------------------------------------------------- #
# validation
# --------------------------------------------------------------------------- #
def validate_against_genome(scan: dict, genome_sizes: dict, on_unknown: str,
                            on_out_of_range: str, chrom_sizes_path: Path) -> dict:
    """Check the GTF really belongs to this genome before anything is published.

    This is the guard that turns FASTA/GTF pairing from an assumption into a
    checked fact: a release bump that changes the contig set, or a GTF flavour
    that annotates patches, is caught here instead of surfacing as silently
    unmapped features months later.
    """
    genome_bare = {strip_chr(t): (t, s) for t, s in genome_sizes.items()}
    gtf_bare = {strip_chr(r.decode("ascii", "replace")) for r in scan["per_seq_in"]}

    # Naming style first: it produces by far the most confusing failure mode.
    gtf_prefixed = sum(1 for r in scan["per_seq_in"] if r.startswith(b"chr"))
    genome_prefixed = sum(1 for t in genome_sizes if t.startswith("chr"))
    gtf_style = "chr-prefixed" if gtf_prefixed > len(scan["per_seq_in"]) / 2 else "bare"
    genome_style = "chr-prefixed" if genome_prefixed > len(genome_sizes) / 2 else "bare"
    if gtf_style != genome_style:
        log(f"NOTE: GTF seqnames are {gtf_style} and the genome is {genome_style}; "
            "matching after normalising away the 'chr' prefix")

    unknown = sorted(gtf_bare - set(genome_bare))
    if unknown:
        counts = {}
        for raw, n in scan["per_seq_in"].items():
            bare = strip_chr(raw.decode("ascii", "replace"))
            if bare in set(unknown):
                counts[bare] = counts.get(bare, 0) + n
        shown = sorted(counts.items(), key=lambda kv: -kv[1])[:20]
        msg = (f"{len(unknown)} seqname(s) in the GTF are absent from {chrom_sizes_path.name} "
               f"({sum(counts.values()):,} features). Top: "
               + ", ".join(f"{k}({v:,})" for k, v in shown))
        if on_unknown == "fail":
            print(f"[make_sexaware_gtf] ERROR: {msg}\n"
                  "  The GTF and the genome do not match. Most often this means a\n"
                  "  chr_patch_hapl_scaff GTF was paired with a primary_assembly FASTA,\n"
                  "  or an Ensembl release that predates the GTF flavour split (e.g. 75).\n"
                  "  Use the 'primary' flavour in config/annotations.tsv, or pass\n"
                  "  --on-unknown-seqname drop to restrict the annotation to the genome.",
                  file=sys.stderr)
            sys.exit(2)
        log(f"{'dropping' if on_unknown == 'drop' else 'WARNING: keeping'} features on "
            f"unknown seqnames -- {msg}")

    out_of_range = []
    for raw, end in scan["max_end"].items():
        bare = strip_chr(raw.decode("ascii", "replace"))
        if bare in genome_bare and end > genome_bare[bare][1]:
            out_of_range.append((bare, end, genome_bare[bare][1]))
    if out_of_range:
        detail = ", ".join(f"{c}: max_end={e:,} > length={L:,}" for c, e, L in out_of_range[:20])
        if on_out_of_range == "fail":
            print(f"[make_sexaware_gtf] ERROR: {len(out_of_range)} contig(s) carry features past "
                  f"the end of the sequence -- {detail}\n"
                  "  This is unambiguous corruption: the GTF belongs to a different assembly\n"
                  "  version than the FASTA. Re-pin both in config/sources.tsv and\n"
                  "  config/annotations.tsv.", file=sys.stderr)
            sys.exit(2)
        log(f"WARNING: {len(out_of_range)} contig(s) with out-of-range features -- {detail}")

    missing = sorted(set(genome_bare) - gtf_bare)
    return {
        "seqnames_in_gtf": len(gtf_bare),
        "seqnames_in_genome": len(genome_bare),
        "seqnames_in_gtf_not_in_genome": unknown,
        "n_seqnames_in_genome_not_in_gtf": len(missing),
        "out_of_range_contigs": [{"seqname": c, "max_end": e, "length": L}
                                 for c, e, L in out_of_range],
    }


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    src = parser.add_argument_group("annotation source")
    src.add_argument("--assembly", help="Registered assembly key (see config/assemblies.tsv).")
    src.add_argument("--gtf", required=True, help="Source GTF (.gtf/.gtf.gz).")
    src.add_argument("--annotations", default=str(REPO_ROOT / "config" / "annotations.tsv"))
    src.add_argument("--registry", default=str(REPO_ROOT / "config" / "assemblies.tsv"))
    src.add_argument("--annotation-source", default=None,
                     help="Provider label. Default: inferred from config/annotations.tsv.")
    src.add_argument("--annotation-release", default=None,
                     help="Release label. Default: inferred from config/annotations.tsv.")

    filt = parser.add_argument_group("filtering")
    filt.add_argument("--variant", required=True, choices=list(VARIANTS) + ["all"],
                      help="Which annotation variant(s) to write. 'all' writes the three.")
    filt.add_argument("--masked-bed", default=None,
                      help="PREFERRED filter input: the masked.bed mask_genome.py wrote for the "
                           "paired FASTA. Default: resolved from --outdir/--prefix.")
    filt.add_argument("--par-bed", default=None,
                      help="Fallback when no masked.bed exists. Cannot see --extra-mask regions.")
    filt.add_argument("--x-contig", default=None, help="X contig name. Default from registry.")
    filt.add_argument("--y-contig", default=None, help="Y contig name. Default from registry.")
    filt.add_argument("--overlap-rule", choices=["any", "contained"], default="any",
                      help="'any': a feature touching the mask is a hit (default). "
                           "'contained': only features wholly inside it.")
    filt.add_argument("--drop-unit", choices=["gene", "transcript", "feature"], default="gene",
                      help="Granularity of removal. 'gene' (default) is the only one that "
                           "cannot leave a transcript model missing exons.")
    filt.add_argument("--on-partial-overlap", choices=["warn", "fail"], default="warn",
                      help="What to do when a gene straddles a mask boundary (the PAB).")
    filt.add_argument("--extra-drop-bed", action="append", default=[],
                      help="Extra BED of regions whose features to drop (repeatable).")
    filt.add_argument("--allow-approximate-mouse-par", action="store_true",
                      help="Permit XY-noYPAR filtering with an 'approximate' (mouse) PAR.")

    val = parser.add_argument_group("validation")
    val.add_argument("--chrom-sizes", default=None,
                     help="chrom.sizes of the paired genome. Default: resolved from --outdir.")
    val.add_argument("--on-unknown-seqname", choices=["fail", "warn", "drop"], default="fail",
                     help="Features on contigs absent from the genome. Default fail.")
    val.add_argument("--on-out-of-range", choices=["fail", "warn"], default="fail",
                     help="Features past the end of their contig. Default fail.")

    out = parser.add_argument_group("output")
    out.add_argument("--outdir", required=True)
    out.add_argument("--prefix", default=None, help="Output basename. Default: <assembly>.")
    out.add_argument("--no-gzip", action="store_true", help="Write plain .gtf (a human GTF is ~1.4 GB).")
    out.add_argument("--emit-t2g", action="store_true",
                     help="Also write a transcript->gene table for salmon/kallisto + tximport.")
    out.add_argument("--force", action="store_true")
    out.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    global VERBOSE
    VERBOSE = not args.quiet

    gtf_path = Path(args.gtf)
    if not gtf_path.exists():
        parser.error(f"--gtf {gtf_path} does not exist")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- registry lookups ------------------------------------------------ #
    x_contig = args.x_contig or "chrX"
    y_contig = args.y_contig or "chrY"
    par_status = "verified"
    registry_par_bed = None
    prefix = args.prefix
    if args.assembly:
        registry = read_assembly_registry(Path(args.registry))
        if args.assembly not in registry:
            parser.error(f"Unknown assembly '{args.assembly}'. Known: {', '.join(registry)}")
        rec = registry[args.assembly]
        x_contig = args.x_contig or rec["x_contig"]
        y_contig = args.y_contig or rec["y_contig"]
        par_status = rec["par_status"]
        if rec["par_bed"] != "NA":
            registry_par_bed = REPO_ROOT / rec["par_bed"]
        prefix = prefix or args.assembly
    else:
        if args.x_contig is None or args.y_contig is None:
            parser.error("Custom mode (no --assembly) requires --x-contig and --y-contig.")
        prefix = prefix or gtf_path.name.split(".gtf")[0]

    # Match the GTF back to its registry row, so source/release/restrict policy
    # come from data rather than from a filename guess.
    ann_source, ann_release, restrict = args.annotation_source, args.annotation_release, False
    ann_flavour = None
    ann_rows = read_annotation_registry(Path(args.annotations)) if Path(args.annotations).exists() else []
    for row in ann_rows:
        if Path(row["url"]).name == gtf_path.name:
            ann_source = ann_source or row["source"]
            ann_release = ann_release or row["release"]
            ann_flavour = row.get("flavour")
            restrict = row.get("restrict_to_genome", "0") == "1"
            break
    if ann_release is None:
        # Ensembl names files <Species>.<ASM>.<REL>.gtf.gz
        parts = gtf_path.name.split(".")
        ann_release = next((p for p in parts if p.isdigit()), "NA")
    ann_source = ann_source or "unknown"
    on_unknown = args.on_unknown_seqname
    if restrict and on_unknown == "fail":
        log(f"config/annotations.tsv marks {gtf_path.name} restrict_to_genome=1; "
            "restricting the annotation to the genome contigs")
        on_unknown = "drop"

    variants = list(VARIANTS) if args.variant == "all" else [args.variant]
    if "XY-noYPAR" in variants and par_status == "approximate" and not args.allow_approximate_mouse_par:
        parser.error(
            f"PAR for {args.assembly} is marked APPROXIMATE (mouse). Refusing XY-noYPAR filtering.\n"
            "The Y-PAR drop set would come from a DERIVED, unverified interval.\n"
            "Supply a curated --par-bed/--masked-bed, or pass --allow-approximate-mouse-par."
        )

    log(f"source GTF  : {gtf_path.name}  [{ann_source}-{ann_release}]")
    log(f"assembly    : {args.assembly or '(custom)'}  |  X={x_contig}  Y={y_contig}  PAR={par_status}")
    log(f"variants    : {', '.join(variants)}")
    log("computing source GTF md5 ...")
    gtf_md5 = md5_of_file(gtf_path)

    exit_code = 0
    for variant in variants:
        complement = VARIANT_COMPLEMENT[variant]
        token = VARIANT_TOKEN[variant]
        stem = f"{prefix}.{token}.{ann_source}-{ann_release}"
        suffix = "gtf" if args.no_gzip else "gtf.gz"
        out_gtf = outdir / f"{stem}.{suffix}"
        out_md5 = outdir / f"{stem}.{suffix}.md5"
        out_meta = outdir / f"{stem}.metadata.json"
        out_dropped = outdir / f"{stem}.dropped_genes.tsv"
        out_t2g = outdir / f"{stem}.t2g.tsv" if args.emit_t2g else None
        if out_gtf.exists() and not args.force:
            parser.error(f"{out_gtf} exists; pass --force to overwrite.")

        # chrom.sizes of the paired genome: the annotation cannot be validated
        # without the genome it is meant to match.
        sizes_path = Path(args.chrom_sizes) if args.chrom_sizes else \
            outdir / f"{prefix}.{complement}.chrom.sizes"
        if not sizes_path.exists():
            parser.error(
                f"{sizes_path} not found. A sex-aware GTF is defined relative to its sex-aware\n"
                f"genome, so build it first:\n"
                f"  python scripts/mask_genome.py --assembly {args.assembly or '<KEY>'} "
                f"--fasta <source.fa.gz> --complement both --outdir {outdir} --gzip\n"
                f"or point --chrom-sizes at the matching file."
            )
        genome_sizes = read_chrom_sizes(sizes_path)

        masked_bed = Path(args.masked_bed) if args.masked_bed else None
        if masked_bed is None and variant != "XY":
            candidate = outdir / f"{prefix}.{complement}.masked.bed"
            masked_bed = candidate if candidate.exists() else None
        par_bed = Path(args.par_bed) if args.par_bed else registry_par_bed

        log(f"=== {variant}: {out_gtf.name} ===")
        mask, mask_notes = resolve_mask_intervals(
            variant=variant, masked_bed=masked_bed if variant != "XY" else None,
            par_bed=par_bed, y_contig=y_contig, genome_sizes=genome_sizes,
            extra_drop_beds=args.extra_drop_bed,
        )
        for note in mask_notes:
            log(f"  {note}")

        log("  pass 1/2: scanning annotation ...")
        scan = scan_gtf(gtf_path, mask, genome_sizes, on_unknown,
                        args.overlap_rule, args.drop_unit, y_contig)
        log(f"  {scan['features_in']:,} input features over "
            f"{len(scan['per_seq_in'])} seqname(s)")
        if scan["malformed"]:
            log(f"  WARNING: skipped {scan['malformed']} malformed line(s)")

        validation = validate_against_genome(scan, genome_sizes, on_unknown,
                                             args.on_out_of_range, sizes_path)

        n_partial = sum(1 for rec in scan["genes"].values() if rec.partial)
        if n_partial:
            names = [r.name or "?" for r in scan["genes"].values() if r.partial][:10]
            log(f"  WARNING: {n_partial} gene(s) straddle a mask boundary "
                f"(the PAB): {', '.join(names)}")
            if args.on_partial_overlap == "fail":
                print("[make_sexaware_gtf] ERROR: partial mask overlap with "
                      "--on-partial-overlap fail", file=sys.stderr)
                sys.exit(3)

        predicted = predict_output(scan, args.drop_unit)
        log(f"  drop set: {len(scan['drop_genes'])} gene(s), "
            f"{len(scan['drop_transcripts'])} transcript(s)")

        log("  pass 2/2: writing ...")
        result = write_gtf(gtf_path, out_gtf, scan, mask, scan["header"] + build_header(
            variant=variant, assembly=args.assembly, prefix=prefix, complement=complement,
            gtf_path=gtf_path, gtf_md5=gtf_md5, ann_source=ann_source, ann_release=ann_release,
            masked_bed=masked_bed if variant != "XY" else None, outdir=outdir,
            overlap_rule=args.overlap_rule, drop_unit=args.drop_unit,
            features_in=scan["features_in"], features_out=predicted,
            genes_dropped=len(scan["drop_genes_sexaware"]),
            genes_dropped_off_genome=len(scan["drop_genes_off_genome"]),
            mask_notes=mask_notes,
            on_unknown=on_unknown,
        ), args.drop_unit, args.overlap_rule, x_contig, y_contig,
            args.emit_t2g, out_t2g)

        n_dropped_rows = write_dropped_genes(out_dropped, scan)

        if result["t2g_collisions"]:
            log(f"  WARNING: {len(result['t2g_collisions'])} transcript_id(s) map to more than "
                f"one gene_id, so the t2g table collapsed them: "
                f"{', '.join(result['t2g_collisions'][:5])}")

        # Measured on GRCh37 r75 and r87: the Y PAR is all-N in Ensembl's source
        # FASTA, so the genebuild annotates nothing there (the first gene on Y
        # starts at 2,652,790, past the 2,649,520 PAR1 end) and XY-noYPAR comes
        # out identical to XY. Say so rather than shipping a look-alike file.
        # Only sex-aware drops make XY-noYPAR differ from XY in substance; an
        # off-genome restriction applies identically to both variants.
        identical_to_xy = variant == "XY-noYPAR" and len(scan["drop_genes_sexaware"]) == 0
        if identical_to_xy:
            log("  NOTE: nothing was dropped, so this XY-noYPAR annotation is identical to the "
                "standard XY one. The Y PAR of this assembly carries no annotation (its source "
                "FASTA already has that region as N), so there is nothing for the filter to do.")

        # ---- post-conditions --------------------------------------------- #
        problems = []
        if predicted != -1 and result["features_out"] != predicted:
            problems.append(f"predicted {predicted:,} output features, wrote "
                            f"{result['features_out']:,}")
        if variant == "XX" and result["y_features_out"] != 0:
            problems.append(f"XX still carries {result['y_features_out']:,} features on {y_contig}")
        if args.overlap_rule == "any" and result["features_on_masked_sequence_kept"] != 0:
            problems.append(f"{result['features_on_masked_sequence_kept']:,} surviving features "
                            "still overlap the mask")
        x_in = sum(n for raw, n in scan["per_seq_in"].items()
                   if contig_matches(strip_chr(raw.decode("ascii", "replace")), x_contig))
        if result["x_features_out"] != x_in:
            problems.append(f"{x_contig} features changed: {x_in:,} in -> "
                            f"{result['x_features_out']:,} out. chrX must never lose a feature: "
                            "the usual cause is a drop set keyed on a bare gene_id where the X and "
                            "Y PAR copies share one. Key it on (seqname, gene_id).")
        if problems:
            for p in problems:
                print(f"[make_sexaware_gtf] POST-CONDITION FAILED: {p}", file=sys.stderr)
            sys.exit(4)

        log("  computing md5 ...")
        checksum = md5_of_file(out_gtf)
        out_md5.write_text(f"{checksum}  {out_gtf.name}\n")

        metadata = {
            "tool": "make_sexaware_gtf.py",
            "tool_version": TOOL_VERSION,
            "repo_revision": git_describe(),
            "assembly": args.assembly,
            "variant": variant,
            "paired_complement": complement,
            "paired_genome_fasta": paired_fasta_name(outdir, prefix, complement),
            "paired_genome_fasta_md5": paired_fasta_md5(outdir, prefix, complement),
            "chrom_sizes": sizes_path.name,
            "annotation_source": ann_source,
            "annotation_release": ann_release,
            "annotation_flavour": ann_flavour or "unknown",
            "source_gtf": gtf_path.name,
            "source_gtf_md5": gtf_md5,
            "mask_bed": masked_bed.name if (masked_bed and variant != "XY") else None,
            "mask_bed_md5": md5_of_file(masked_bed) if (masked_bed and variant != "XY") else None,
            "par_bed": str(par_bed) if par_bed else None,
            "par_status": par_status,
            "mask_notes": mask_notes,
            "extra_drop_beds": list(args.extra_drop_bed),
            "x_contig": x_contig,
            "y_contig": y_contig,
            "overlap_rule": args.overlap_rule,
            "drop_unit": args.drop_unit,
            "on_unknown_seqname": on_unknown,
            "features_in": scan["features_in"],
            "features_out": result["features_out"],
            "features_dropped": scan["features_in"] - result["features_out"],
            "genes_dropped": len(scan["drop_genes"]),
            "genes_dropped_sexaware": len(scan["drop_genes_sexaware"]),
            "genes_dropped_off_genome": len(scan["drop_genes_off_genome"]),
            "restrict_to_genome": bool(restrict),
            "genes_in_dropped_table": n_dropped_rows,
            "genes_partially_overlapping": n_partial,
            "transcripts_out": result["n_transcripts"] if args.emit_t2g else None,
            "features_on_masked_sequence_kept": result["features_on_masked_sequence_kept"],
            "t2g_collisions": result["t2g_collisions"],
            "identical_to_standard_xy": identical_to_xy,
            "features_on_x_contig_in": x_in,
            "features_on_x_contig_out": result["x_features_out"],
            "features_on_y_contig_out": result["y_features_out"],
            "validation": validation,
            "output_gtf": out_gtf.name,
            "output_gtf_bytes": out_gtf.stat().st_size,
            "output_gtf_md5": checksum,
            "method": (
                "A sex-aware GTF contains no feature overlapping N-masked sequence in the "
                "corresponding sex-aware FASTA. XX: the whole Y contig is removed. "
                "XY: upstream annotation passed through. XY-noYPAR: features on the Y-side PAR "
                "are removed. The drop set is keyed on (seqname, gene_id), never on gene_id "
                "alone, so that a release which gives an X/Y PAR gene pair one shared id cannot "
                "lose the X copy."
            ),
        }
        out_meta.write_text(json.dumps(metadata, indent=2) + "\n")

        dropped_note = f"{len(scan['drop_genes_sexaware'])} gene(s) dropped"
        if scan["drop_genes_off_genome"]:
            dropped_note += (f" + {len(scan['drop_genes_off_genome'])} off-genome"
                             " (not sex-aware filtering)")
        print(f"[{variant}] {out_gtf}  ({scan['features_in']:,} -> "
              f"{result['features_out']:,} features; {dropped_note}; "
              f"md5 {checksum})", file=sys.stderr)

    sys.exit(exit_code)


def paired_fasta_name(outdir: Path, prefix: str, complement: str) -> str | None:
    """Basename of the sex-aware FASTA this GTF is meant to be used with."""
    for suffix in ("fa.gz", "fa"):
        candidate = outdir / f"{prefix}.{complement}.{suffix}"
        if candidate.exists():
            return candidate.name
    return f"{prefix}.{complement}.fa.gz"


def paired_fasta_md5(outdir: Path, prefix: str, complement: str) -> str | None:
    """md5 of the paired FASTA, read from the sidecar mask_genome.py wrote.

    Recording it binds this GTF to one exact genome build, so a user can prove
    the pair they hold is the pair that was tested together.
    """
    sidecar = outdir / f"{prefix}.{complement}.fa.md5"
    if sidecar.exists():
        text = sidecar.read_text().split()
        return text[0] if text else None
    return None


def build_header(**kw) -> list:
    """Provenance comment lines, prefixed '#!' so any `grep -v '^#'` drops them.

    Deliberately carries NO build timestamp: the output stays byte-reproducible,
    which is what makes the md5 sidecar a real identity. The build time lives in
    the metadata JSON instead.
    """
    masked_bed = kw["masked_bed"]
    lines = [
        f"#!sexaware-variant {kw['variant']}",
        f"#!sexaware-assembly {kw['assembly'] or '(custom)'}",
        f"#!sexaware-paired-genome {paired_fasta_name(kw['outdir'], kw['prefix'], kw['complement'])}",
        f"#!sexaware-paired-genome-md5 "
        f"{paired_fasta_md5(kw['outdir'], kw['prefix'], kw['complement']) or 'NA'}",
        f"#!sexaware-annotation-source {kw['ann_source']}",
        f"#!sexaware-annotation-release {kw['ann_release']}",
        f"#!sexaware-source-gtf {kw['gtf_path'].name}",
        f"#!sexaware-source-gtf-md5 {kw['gtf_md5']}",
        f"#!sexaware-mask-bed {masked_bed.name if masked_bed else 'none'}",
        f"#!sexaware-mask-bed-md5 {md5_of_file(masked_bed) if masked_bed else 'NA'}",
        f"#!sexaware-overlap-rule {kw['overlap_rule']}",
        f"#!sexaware-drop-unit {kw['drop_unit']}",
        f"#!sexaware-on-unknown-seqname {kw['on_unknown']}",
        f"#!sexaware-features-in {kw['features_in']}",
        f"#!sexaware-features-out {kw['features_out']}",
        f"#!sexaware-genes-dropped {kw['genes_dropped']}",
        # Only emitted when non-zero, so the common case stays uncluttered and
        # the line's presence is itself the signal that restriction happened.
        *([f"#!sexaware-genes-dropped-off-genome {kw['genes_dropped_off_genome']}"]
          if kw.get("genes_dropped_off_genome") else []),
        f"#!sexaware-tool make_sexaware_gtf.py {TOOL_VERSION}",
        f"#!sexaware-repo {git_describe()}",
    ]
    for note in kw["mask_notes"]:
        lines.append(f"#!sexaware-note {note}")
    return [line.encode("utf-8") + b"\n" for line in lines]


if __name__ == "__main__":
    main()
