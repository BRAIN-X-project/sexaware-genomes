#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Infer a sample's sex-chromosome complement before choosing XX vs XY reference.

Data product #4 (the gatekeeper). Estimates chrX/autosome and chrY/autosome
normalised coverage from a BAM's index statistics and reports a call:

  XX   : chrY/autosome ~ 0            (no Y)
  XY   : chrY/autosome ~ 0.5 diploid-equivalent, chrX/autosome ~ 0.5
  other: ratios inconsistent with XX/XY (flag: XXY, XYY, X0/Turner, mosaicism,
         contamination, sample swap) -> inspect before aligning.

This is a fast idxstats-based screen (reads-per-base normalised by mappable
length is more robust; pass --mappable-sizes to use mappable bp per contig from
your callable BEDs instead of raw contig length). For production QC also inspect
X-heterozygosity (a true XX has heterozygous SNPs across chrX non-PAR; a true XY
does not), e.g. with XYalign or bcftools.

Requires samtools on PATH (uses `samtools idxstats`).
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time

_T0 = time.time()
VERBOSE = True


def log(msg: str) -> None:
    """Timestamped progress line to stderr when VERBOSE is on (result goes to stdout)."""
    if VERBOSE:
        print(f"[infer_sex {time.strftime('%H:%M:%S')} +{time.time() - _T0:5.1f}s] {msg}",
              file=sys.stderr, flush=True)


def strip(name):
    return name[3:] if name.startswith("chr") else name


def idxstats(bam):
    out = subprocess.check_output(["samtools", "idxstats", bam], text=True)
    rows = {}
    for line in out.splitlines():
        name, length, mapped, _unmapped = line.split("\t")
        rows[name] = (int(length), int(mapped))
    return rows


def load_mappable_sizes(path):
    sizes = {}
    for line in open(path):
        if line.strip():
            name, val = line.split()[:2]
            sizes[name] = int(val)
    return sizes


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam", required=True)
    p.add_argument("--x-contig", default="chrX")
    p.add_argument("--y-contig", default="chrY")
    p.add_argument("--mappable-sizes", default=None,
                   help="Optional 2-col file (contig, mappable_bp) to normalise by "
                        "mappable length instead of raw contig length.")
    p.add_argument("--xx-max-y", type=float, default=0.05,
                   help="Max chrY/autosome ratio still called XX. Default 0.05.")
    p.add_argument("--xy-min-y", type=float, default=0.15,
                   help="Min chrY/autosome ratio to call a Y present. Default 0.15.")
    p.add_argument("--quiet", action="store_true", help="Suppress the progress log (result still on stdout).")
    args = p.parse_args()

    global VERBOSE
    VERBOSE = not args.quiet

    log(f"running `samtools idxstats {args.bam}` ...")
    stats = idxstats(args.bam)
    norm = load_mappable_sizes(args.mappable_sizes) if args.mappable_sizes else None
    log(f"idxstats: {len(stats)} contig(s)"
        + (f"; normalising by mappable bp from {args.mappable_sizes}" if norm else "; normalising by contig length"))

    def depth(name):
        length, mapped = stats[name]
        denom = norm.get(name, length) if norm else length
        return mapped / denom if denom else 0.0

    auto = [n for n in stats if strip(n).isdigit()]
    if not auto:
        sys.exit("No autosomes found in idxstats; check contig naming.")
    auto_depth = sum(depth(n) for n in auto) / len(auto)
    log(f"mean autosomal depth over {len(auto)} contig(s): {auto_depth:.4g}")

    x_name = next((n for n in stats if strip(n) == strip(args.x_contig)), None)
    y_name = next((n for n in stats if strip(n) == strip(args.y_contig)), None)
    x_ratio = (depth(x_name) / auto_depth) if (x_name and auto_depth) else float("nan")
    y_ratio = (depth(y_name) / auto_depth) if (y_name and auto_depth) else 0.0
    log(f"X contig {x_name}: ratio {x_ratio:.3f}  |  Y contig {y_name}: ratio {y_ratio:.3f}")

    if y_ratio <= args.xx_max_y and x_ratio >= 0.8:
        call, ref = "XX", "use the XX (whole-Y-masked) reference"
    elif y_ratio >= args.xy_min_y and 0.3 <= x_ratio <= 0.75:
        call, ref = "XY", "use the XY (Y-PAR-masked) reference"
    else:
        call, ref = "OTHER/CHECK", "ambiguous — inspect (XXY/XYY/X0/mosaic/swap/contam) before aligning"

    print(f"sample_bam\t{args.bam}")
    print(f"chrX/autosome\t{x_ratio:.3f}")
    print(f"chrY/autosome\t{y_ratio:.3f}")
    print(f"call\t{call}")
    print(f"recommendation\t{ref}")


if __name__ == "__main__":
    main()
