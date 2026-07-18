#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Convert a GenMap per-k-mer bedGraph into Umap-style single/multi-read tracks.

GenMap emits, for each position i that starts a valid k-mer, the value
1/frequency (== 1.0 when the k-mer is unique). This script produces:

  single_read(i) = 1.0 if the k-mer starting at i is unique (freq == 1), else 0.0
                   (Umap "single-read mappability": is this exact k-mer unique?)
  multi_read(i)  = fraction of the up-to-k k-mers overlapping position i whose
                   starting k-mer is unique — i.e. a length-k trailing mean of the
                   single_read indicator (Umap "multi-read mappability").

The multi-read value at base i averages the single-read indicator over the k
k-mer start positions [i-k+1 .. i] that cover base i (clipped at contig starts).
Both outputs are written as bedGraph and can be fed to bedGraphToBigWig.

Runs one contig at a time to keep memory bounded; uses a cumulative-sum sliding
window so it is O(genome length).
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np

_T0 = time.time()
VERBOSE = True


def log(msg: str) -> None:
    """Timestamped progress line to stderr when VERBOSE is on."""
    if VERBOSE:
        print(f"[genmap2umap {time.strftime('%H:%M:%S')} +{time.time() - _T0:5.1f}s] {msg}",
              file=sys.stderr, flush=True)


def read_chrom_sizes(path):
    sizes = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            name, length = line.split()[:2]
            sizes[name] = int(length)
    return sizes


def iter_genmap_bedgraph(path):
    """Yield (chrom, start, end, value) from a GenMap bedGraph, in file order."""
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("track", "#")):
                continue
            chrom, start, end, value = line.rstrip("\n").split("\t")[:4]
            yield chrom, int(start), int(end), float(value)


def unique_indicator_array(chrom_length, intervals):
    """Build a per-position single-read indicator for one contig.

    `intervals` is a list of (start, end, value) covering k-mer START positions.
    Positions with value == 1.0 are unique (indicator 1); positions never emitted
    by GenMap (inside Ns or without a valid k-mer) stay 0.
    """
    indicator = np.zeros(chrom_length, dtype=np.float32)
    for start, end, value in intervals:
        if value >= 1.0:  # unique k-mer
            indicator[start:end] = 1.0
    return indicator


def multi_read_from_single(single, kmer):
    """Trailing length-k mean of the single-read indicator (Umap multi-read).

    multi[i] = mean(single[max(0, i-k+1) .. i]) using a cumulative sum; the
    denominator is the actual (clipped) window length so contig starts are not
    penalised for missing upstream positions.
    """
    n = single.shape[0]
    if n == 0:
        return single
    cumulative = np.zeros(n + 1, dtype=np.float64)
    np.cumsum(single, dtype=np.float64, out=cumulative[1:])
    idx = np.arange(n)
    lo = np.maximum(0, idx - kmer + 1)
    window_sum = cumulative[idx + 1] - cumulative[lo]
    window_len = (idx + 1 - lo).astype(np.float64)
    return (window_sum / window_len).astype(np.float32)


def write_bedgraph_runs(handle, chrom, values):
    """Write a per-base array as run-length-encoded bedGraph (0 runs skipped)."""
    n = values.shape[0]
    if n == 0:
        return
    change = np.nonzero(np.diff(values))[0] + 1
    starts = np.concatenate(([0], change))
    ends = np.concatenate((change, [n]))
    for start, end in zip(starts, ends):
        value = values[start]
        if value == 0.0:
            continue
        handle.write(f"{chrom}\t{start}\t{end}\t{value:.6g}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genmap-bedgraph", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--kmer", type=int, required=True)
    parser.add_argument("--single-out", required=True)
    parser.add_argument("--multi-out", required=True)
    parser.add_argument("--quiet", action="store_true", help="Suppress per-contig progress log.")
    args = parser.parse_args()

    global VERBOSE
    VERBOSE = not args.quiet

    sizes = read_chrom_sizes(args.chrom_sizes)
    log(f"k={args.kmer}  |  {len(sizes)} contig(s) from {args.chrom_sizes}")

    # Group GenMap intervals by contig (input is sorted per contig by GenMap).
    log(f"reading GenMap bedGraph {args.genmap_bedgraph} ...")
    per_chrom = {}
    for chrom, start, end, value in iter_genmap_bedgraph(args.genmap_bedgraph):
        per_chrom.setdefault(chrom, []).append((start, end, value))
    log(f"loaded intervals for {len(per_chrom)} contig(s); deriving tracks")

    with open(args.single_out, "w") as single_handle, open(args.multi_out, "w") as multi_handle:
        for chrom, length in sizes.items():
            intervals = per_chrom.get(chrom, [])
            single = unique_indicator_array(length, intervals)
            multi = multi_read_from_single(single, args.kmer)
            write_bedgraph_runs(single_handle, chrom, single)
            write_bedgraph_runs(multi_handle, chrom, multi)
            unique_bp = int(single.sum())
            log(f"  {chrom}: {length:,} bp, {unique_bp:,} uniquely-mappable "
                f"({100.0 * unique_bp / length if length else 0:.1f}%)")
    log(f"wrote {args.single_out} and {args.multi_out}")


if __name__ == "__main__":
    main()
