#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Emit ploidy / PAR helper files for sex-aware variant calling.

Data product #2. Given an assembly's PAR BED and its chrom.sizes, write the
region files callers need so that, in an XY sample, the non-PAR X and non-PAR Y
are treated as HAPLOID while PAR + autosomes stay DIPLOID; in an XX sample the X
is diploid and Y is absent.

Outputs (under --outdir, prefixed by --label):
  <label>.PAR.bed                 PAR intervals (chrX & chrY)
  <label>.chrX_nonPAR.bed         diploid-in-XX / haploid-in-XY region of X
  <label>.chrY_nonPAR.bed         haploid-in-XY region of Y (empty concept in XX)
  <label>.XY.haploid_regions.bed  chrX_nonPAR + chrY_nonPAR (ploidy 1 in XY)
  <label>.XY.gatk_ploidy.md       ready-to-adapt GATK / DeepVariant snippets

This produces region files, not a full pipeline; you still choose the caller.
"""
from __future__ import annotations

import argparse
from pathlib import Path


def read_sizes(path):
    sizes = {}
    for line in Path(path).read_text().splitlines():
        if line.strip():
            name, length = line.split()[:2]
            sizes[name] = int(length)
    return sizes


def read_bed(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.split("\t")
        rows.append((f[0], int(f[1]), int(f[2])))
    return rows


def match(contig, target):
    strip = lambda s: s[3:] if s.startswith("chr") else s
    return strip(contig) == strip(target)


def complement(contig, length, par_rows):
    """Return the non-PAR intervals of one contig as (start, end) list."""
    pars = sorted((s, e) for c, s, e in par_rows if match(c, contig))
    out, cursor = [], 0
    for s, e in pars:
        if s > cursor:
            out.append((cursor, s))
        cursor = max(cursor, e)
    if cursor < length:
        out.append((cursor, length))
    return out


def write_bed(path, contig, intervals, name):
    with open(path, "w") as h:
        for s, e in intervals:
            h.write(f"{contig}\t{s}\t{e}\t{name}\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--par-bed", required=True)
    p.add_argument("--chrom-sizes", required=True)
    p.add_argument("--label", required=True, help="e.g. GRCh38")
    p.add_argument("--x-contig", default="chrX")
    p.add_argument("--y-contig", default="chrY")
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    sizes = read_sizes(args.chrom_sizes)
    par = read_bed(args.par_bed)

    x_len = next(v for k, v in sizes.items() if match(k, args.x_contig))
    y_len = next((v for k, v in sizes.items() if match(k, args.y_contig)), None)
    x_name = next(k for k in sizes if match(k, args.x_contig))
    y_name = next((k for k in sizes if match(k, args.y_contig)), args.y_contig)

    # PAR passthrough (both chromosomes).
    with open(outdir / f"{args.label}.PAR.bed", "w") as h:
        for c, s, e in par:
            h.write(f"{c}\t{s}\t{e}\tPAR\n")

    x_nonpar = complement(args.x_contig, x_len, par)
    write_bed(outdir / f"{args.label}.chrX_nonPAR.bed", x_name, x_nonpar, "X_nonPAR")

    y_nonpar = complement(args.y_contig, y_len, par) if y_len else []
    if y_len:
        write_bed(outdir / f"{args.label}.chrY_nonPAR.bed", y_name, y_nonpar, "Y_nonPAR")

    with open(outdir / f"{args.label}.XY.haploid_regions.bed", "w") as h:
        for s, e in x_nonpar:
            h.write(f"{x_name}\t{s}\t{e}\tX_nonPAR_haploid\n")
        for s, e in y_nonpar:
            h.write(f"{y_name}\t{s}\t{e}\tY_nonPAR_haploid\n")

    (outdir / f"{args.label}.XY.gatk_ploidy.md").write_text(f"""# Sex-aware ploidy — {args.label}

## Concept
| sample | region | ploidy |
|--------|--------|--------|
| XX | autosomes, chrX (whole) | 2 |
| XX | chrY | absent (use the XX-masked reference) |
| XY | autosomes, PAR (chrX & chrY) | 2 |
| XY | chrX non-PAR, chrY non-PAR | 1 |

## GATK HaplotypeCaller (XY)
Call the diploid and haploid partitions separately, then merge:

```bash
# diploid: autosomes + PAR
gatk HaplotypeCaller -R {args.label}.XY.fa -I sample.bam \\
  -L autosomes.interval_list -L {args.label}.PAR.bed \\
  --sample-ploidy 2 -O sample.diploid.g.vcf.gz

# haploid: non-PAR X and Y
gatk HaplotypeCaller -R {args.label}.XY.fa -I sample.bam \\
  -L {args.label}.XY.haploid_regions.bed \\
  --sample-ploidy 1 -O sample.haploid.g.vcf.gz
```

## DeepVariant (XY)
```bash
run_deepvariant --model_type=WGS \\
  --ref={args.label}.XY.fa --reads=sample.bam \\
  --haploid_contigs="{y_name}" \\
  --par_regions_bed={args.label}.PAR.bed \\
  --output_vcf=sample.vcf.gz
```
Note: DeepVariant treats the given contigs as haploid except inside
`--par_regions_bed`. For X, restrict/haploid-ise via regions as above.
""")
    print(f"Wrote ploidy/PAR files for {args.label} to {outdir}")


if __name__ == "__main__":
    main()
