#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Threshold a mappability bigWig into a callable/mappable BED.

Data product #1: for each assembly x complement (XX/XY) x k, turn the continuous
mappability bigWig into a BED of positions that pass a score threshold, already
consistent with the sex-aware masking (masked-out Y / Y-PAR read as 0 and are
therefore excluded automatically).

Typical use:
  # strict single-read (uniquely mappable) callable regions at k=100
  build_callable_beds.py --bigwig GRCh38.XY.k100.single_read.bw \\
      --threshold 1.0 --min-length 1 --out GRCh38.XY.k100.callable.bed

  # relaxed multi-read >= 0.9
  build_callable_beds.py --bigwig GRCh38.XY.k100.multi_read.bw \\
      --threshold 0.9 --out GRCh38.XY.k100.mappable_ge0.9.bed

Requires pyBigWig.
"""
from __future__ import annotations

import argparse
import numpy as np


def callable_intervals(values, threshold, min_length):
    """Yield (start, end) runs where values >= threshold and run length >= min_length."""
    passing = np.nan_to_num(values, nan=0.0) >= threshold
    if not passing.any():
        return
    padded = np.concatenate(([False], passing, [False]))
    edges = np.diff(padded.astype(np.int8))
    starts = np.nonzero(edges == 1)[0]
    ends = np.nonzero(edges == -1)[0]
    for start, end in zip(starts, ends):
        if end - start >= min_length:
            yield int(start), int(end)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bigwig", required=True, help="Mappability bigWig (single_read or multi_read).")
    parser.add_argument("--out", required=True, help="Output BED path.")
    parser.add_argument("--threshold", type=float, default=1.0,
                        help="Minimum score to be callable (1.0 = uniquely mappable). Default 1.0.")
    parser.add_argument("--min-length", type=int, default=1, help="Drop callable runs shorter than this.")
    parser.add_argument("--chunk", type=int, default=10_000_000, help="Bases fetched per bigWig read.")
    parser.add_argument("--name", default=None, help="Optional 4th-column label for each interval.")
    args = parser.parse_args()

    import pyBigWig

    bw = pyBigWig.open(args.bigwig)
    label = args.name
    total = 0
    with open(args.out, "w") as out:
        for chrom, length in bw.chroms().items():
            for offset in range(0, length, args.chunk):
                end = min(offset + args.chunk, length)
                values = bw.values(chrom, offset, end, numpy=True)
                for lo, hi in callable_intervals(values, args.threshold, args.min_length):
                    a, b = lo + offset, hi + offset
                    total += b - a
                    if label:
                        out.write(f"{chrom}\t{a}\t{b}\t{label}\n")
                    else:
                        out.write(f"{chrom}\t{a}\t{b}\n")
    bw.close()
    print(f"Wrote {args.out}: {total:,} callable bp (threshold {args.threshold}).")


if __name__ == "__main__":
    main()
