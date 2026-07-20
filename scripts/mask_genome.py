#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build a sex-chromosome-complement (SCC) informed reference genome.

This is a self-contained, organism-agnostic generalisation of the original
`1A_map_hg38_ref.py`. It hard-masks (with Ns) the regions that create X<->Y
cross-mapping, following the consensus of:

  * Olney et al. 2020, Biology of Sex Differences (10.1186/s13293-020-00312-9)
  * "Best practices for improving alignment and variant calling on human sex
    chromosomes" (PMC12190741, AJHG 2025)

Two reference flavours are produced, one per sample sex-chromosome complement:

  XX  (no Y)  -> hard-mask the ENTIRE Y contig.
  XY  (has Y) -> hard-mask ONLY the PAR intervals ON THE Y contig, so PAR reads
                 map unambiguously to the intact X copy. The X is left intact.

You either:
  (a) pass --assembly <KEY> to use a registered assembly + its verified PAR BED
      (config/assemblies.tsv, config/par/*.bed), or
  (b) pass --fasta <FILE> together with the contig names and --par-bed for a
      fully custom genome.

Optional extra masking (--extra-mask BED, repeatable) lets you additionally
hard-mask stricter "best-practices" regions (XTR, ampliconic/palindromic Y,
gametolog genes) on top of the standard masking. By default nothing extra is
masked; the extra BEDs are shipped in config/extra_mask/ so downstream tools can
also use them purely as blacklists.

Outputs (under --outdir):
  <prefix>.<XX|XY>.fa[.gz]      masked FASTA
  <prefix>.<XX|XY>.chrom.sizes  two-column contig sizes
  <prefix>.<XX|XY>.masked.bed   every interval that was set to N
  <prefix>.<XX|XY>.metadata.json provenance + parameters
  <prefix>.<XX|XY>.fa.md5       checksum for verification
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# --------------------------------------------------------------------------- #
# verbose logging (timestamped, to stderr; silence with --quiet)
# --------------------------------------------------------------------------- #
_T0 = time.time()
VERBOSE = True


def log(msg: str) -> None:
    """Print a timestamped progress line to stderr when VERBOSE is on."""
    if VERBOSE:
        print(f"[mask_genome {time.strftime('%H:%M:%S')} +{time.time() - _T0:5.1f}s] {msg}",
              file=sys.stderr, flush=True)


# --------------------------------------------------------------------------- #
# I/O helpers
# --------------------------------------------------------------------------- #
def open_maybe_gzip_text(path: Path):
    """Open a plain or gzip-compressed file for streaming text reads."""
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="ascii", newline="")
    return open(path, "r", encoding="ascii", newline="")


def open_out_text(path: Path):
    """Open a plain or gzip-compressed output for text writes."""
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "wb"), encoding="ascii", newline="")
    return open(path, "w", encoding="ascii", newline="")


def read_assembly_registry(path: Path) -> dict:
    """Parse config/assemblies.tsv, skipping comment lines, keyed by assembly."""
    registry = {}
    with open(path, "r", encoding="utf-8") as handle:
        rows = [line for line in handle if not line.startswith("#") and line.strip()]
    reader = csv.DictReader(rows, delimiter="\t")
    for row in reader:
        registry[row["assembly"]] = row
    return registry


def read_bed_intervals(path: Path) -> list:
    """Read a BED file (>=3 cols) into (chrom, start, end, name) tuples.

    Lines starting with '#' or 'track'/'browser' are ignored. Coordinates are
    0-based half-open, as in the BED spec.
    """
    intervals = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            name = fields[3] if len(fields) > 3 else f"{chrom}:{start}-{end}"
            if end <= start:
                raise ValueError(f"Invalid interval in {path}: {line.strip()}")
            intervals.append((chrom, start, end, name))
    return intervals


# --------------------------------------------------------------------------- #
# chromosome-name matching
# --------------------------------------------------------------------------- #
def strip_chr(name: str) -> str:
    """Normalise a contig token to its bare form (drop leading 'chr')."""
    token = name.split()[0]
    return token[3:] if token.startswith("chr") else token


def contig_matches(fasta_name: str, target: str) -> bool:
    """True if a FASTA header token refers to the same contig as `target`."""
    return strip_chr(fasta_name) == strip_chr(target)


# --------------------------------------------------------------------------- #
# masking core
# --------------------------------------------------------------------------- #
def build_mask_plan(
    par_intervals: list,
    extra_intervals: list,
    x_contig: str,
    y_contig: str,
    complement: str,
    y_length_hint: int | None,
) -> dict:
    """Return {bare_contig: [(start, end, name), ...]} of regions to N-mask.

    For XX the whole Y is scheduled (using y_length_hint if known, otherwise the
    length is filled in while streaming). For XY only the Y-side PAR intervals
    are used. Extra intervals are always added on top when supplied.
    """
    plan: dict = {}

    def add(bare: str, start: int, end: int, name: str):
        plan.setdefault(bare, []).append((start, end, name))

    if complement == "XX":
        # Whole-Y mask; end filled at stream time if hint is None.
        add(strip_chr(y_contig), 0, y_length_hint if y_length_hint else -1, "whole_Y_mask")
    elif complement == "XY":
        # Route the Y-side PAR intervals onto the actual Y contig, regardless of
        # whether the source FASTA names it chrY / Y / an accession (e.g. NCBI
        # T2T NC_060948.1). A PAR row is Y-side if its canonical BED label is
        # 'Y', or if it already matches the real y_contig name.
        for chrom, start, end, name in par_intervals:
            if strip_chr(chrom).upper() == "Y" or contig_matches(chrom, y_contig):
                add(strip_chr(y_contig), start, end, f"{name}_masked")
    else:
        raise ValueError("complement must be 'XX' or 'XY'")

    for chrom, start, end, name in extra_intervals:
        add(strip_chr(chrom), start, end, f"extra:{name}")

    return plan


def mask_fasta_stream(
    source_fasta: Path,
    out_fasta: Path,
    mask_plan: dict,
    allowlist: set | None,
) -> tuple[dict, list]:
    """Stream a FASTA, N-mask scheduled intervals, write the result.

    Returns (lengths, applied_intervals). `lengths` maps output contig name to
    length; `applied_intervals` is the list of (contig, start, end, name) that
    were actually written as N (whole-Y masks get their real end resolved here).

    Every contig is written by default (the reference must contain the whole
    genome); only scheduled intervals are N-masked. If `allowlist` is given,
    contigs whose bare name is not in it are dropped.
    """
    lengths: dict = {}
    applied: list = []
    current_bare = None
    current_name = None
    current_pos = 0  # 0-based position of the next base to be written
    current_masks: list = []
    write_current = False

    with open_maybe_gzip_text(source_fasta) as handle, open_out_text(out_fasta) as out:

        def finalize_contig():
            # Resolve any open-ended whole-contig mask (end == -1) to full length.
            if current_bare is None:
                return
            lengths[current_name] = current_pos
            masked_bp = 0
            for start, end, name in current_masks:
                resolved_end = current_pos if end == -1 else min(end, current_pos)
                if resolved_end > start:
                    applied.append((current_name, start, resolved_end, name))
                    masked_bp += resolved_end - start
            if current_masks:
                log(f"  contig {current_name}: {current_pos:,} bp written, "
                    f"{masked_bp:,} bp N-masked ({len(current_masks)} interval(s))")
            else:
                log(f"  contig {current_name}: {current_pos:,} bp written, unmasked")

        for line in handle:
            if line.startswith(">"):
                finalize_contig()
                raw = line[1:].strip()
                token = raw.split()[0]
                bare = strip_chr(token)
                if allowlist is not None and bare not in allowlist:
                    current_bare = None
                    write_current = False
                    continue
                masks = mask_plan.get(bare)
                current_bare = bare
                current_name = token
                current_pos = 0
                current_masks = sorted(masks) if masks else []
                write_current = True
                out.write(f">{token}\n")
                continue

            if not write_current or current_bare is None:
                continue

            seq = line.strip()
            if not seq:
                continue
            if current_masks:
                seq = _apply_line_masks(seq, current_pos, current_masks)
            out.write(seq + "\n")
            current_pos += len(seq)

        finalize_contig()

    return lengths, applied


def _apply_line_masks(seq: str, line_start: int, masks: list) -> str:
    """N-mask the parts of one sequence line covered by any mask interval."""
    line_end = line_start + len(seq)  # 0-based half-open
    chars = None
    for start, end, _name in masks:
        eff_end = line_end if end == -1 else end
        overlap_start = max(line_start, start)
        overlap_end = min(line_end, eff_end)
        if overlap_start < overlap_end:
            if chars is None:
                chars = list(seq)
            local_start = overlap_start - line_start
            local_end = overlap_end - line_start
            chars[local_start:local_end] = "N" * (local_end - local_start)
    return "".join(chars) if chars is not None else seq


# --------------------------------------------------------------------------- #
# outputs
# --------------------------------------------------------------------------- #
def write_chrom_sizes(path: Path, lengths: dict) -> None:
    with open(path, "w", encoding="ascii") as handle:
        for name, length in lengths.items():
            handle.write(f"{name}\t{length}\n")


def write_masked_bed(path: Path, applied: list) -> None:
    with open(path, "w", encoding="ascii") as handle:
        for chrom, start, end, name in sorted(applied):
            handle.write(f"{chrom}\t{start}\t{end}\t{name}\n")


def md5_of_file(path: Path, chunk: int = 1 << 20) -> str:
    digest = hashlib.md5()
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    src = parser.add_argument_group("genome source")
    src.add_argument("--assembly", help="Registered assembly key (see config/assemblies.tsv).")
    src.add_argument("--fasta", help="Source FASTA (.fa/.fa.gz). Required in custom mode.")
    src.add_argument("--par-bed", help="PAR BED. Overrides the registry BED; required in custom mode for XY.")
    src.add_argument("--x-contig", default=None, help="X contig name (custom mode). Default from registry.")
    src.add_argument("--y-contig", default=None, help="Y contig name (custom mode). Default from registry.")
    src.add_argument("--registry", default=str(REPO_ROOT / "config" / "assemblies.tsv"))

    mask = parser.add_argument_group("masking")
    mask.add_argument("--complement", required=True, choices=["XX", "XY", "both"],
                      help="Sample sex-chromosome complement. 'both' writes both references.")
    mask.add_argument("--extra-mask", action="append", default=[],
                      help="Extra BED of regions to also hard-mask (repeatable): XTR, ampliconic Y, gametologs.")
    mask.add_argument("--allow-approximate-mouse-par", action="store_true",
                      help="Permit XY masking with an 'approximate' (mouse) PAR BED. Off by default.")

    out = parser.add_argument_group("output")
    out.add_argument("--outdir", required=True)
    out.add_argument("--prefix", default=None, help="Output basename. Default: <assembly> or FASTA stem.")
    out.add_argument("--contigs", default=None,
                     help="Optional comma-separated allowlist of contigs to KEEP (bare names, "
                          "e.g. '1,2,...,22,X,Y'). Default keeps every contig in the source FASTA.")
    out.add_argument("--gzip", action="store_true", help="Write the masked FASTA gzip-compressed.")
    out.add_argument("--force", action="store_true")
    out.add_argument("--quiet", action="store_true", help="Suppress the verbose per-contig progress log.")
    args = parser.parse_args()

    global VERBOSE
    VERBOSE = not args.quiet

    # Resolve source + PAR from registry or custom flags.
    par_status = "verified"
    x_contig = args.x_contig or "chrX"
    y_contig = args.y_contig or "chrY"
    par_bed_path = Path(args.par_bed) if args.par_bed else None
    fasta_path = Path(args.fasta) if args.fasta else None
    default_prefix = args.prefix

    if args.assembly:
        registry = read_assembly_registry(Path(args.registry))
        if args.assembly not in registry:
            parser.error(f"Unknown assembly '{args.assembly}'. Known: {', '.join(registry)}")
        rec = registry[args.assembly]
        x_contig = args.x_contig or rec["x_contig"]
        y_contig = args.y_contig or rec["y_contig"]
        par_status = rec["par_status"]
        if par_bed_path is None and rec["par_bed"] != "NA":
            par_bed_path = REPO_ROOT / rec["par_bed"]
        if fasta_path is None:
            parser.error("--assembly still needs --fasta pointing to the downloaded source FASTA.")
        default_prefix = default_prefix or args.assembly
    else:
        if fasta_path is None:
            parser.error("Custom mode requires --fasta.")
        default_prefix = default_prefix or fasta_path.name.split(".fa")[0]

    allowlist = None
    if args.contigs:
        allowlist = {strip_chr(c.strip()) for c in args.contigs.split(",") if c.strip()}

    complements = ["XX", "XY"] if args.complement == "both" else [args.complement]

    if "XY" in complements:
        if par_bed_path is None:
            parser.error("XY masking requires a PAR BED (--par-bed or a registered assembly).")
        if par_status == "approximate" and not args.allow_approximate_mouse_par:
            parser.error(
                f"PAR BED {par_bed_path} is marked APPROXIMATE (mouse). Refusing XY masking.\n"
                "Supply a curated --par-bed, or pass --allow-approximate-mouse-par to proceed anyway."
            )

    log(f"source FASTA : {fasta_path}")
    log(f"assembly     : {args.assembly or '(custom)'}  |  X={x_contig}  Y={y_contig}")
    log(f"PAR BED      : {par_bed_path or '(none)'}  [{par_status}]")
    log(f"complements  : {', '.join(complements)}")

    par_intervals = read_bed_intervals(par_bed_path) if par_bed_path else []
    if par_intervals:
        log(f"read {len(par_intervals)} PAR interval(s)")
    extra_intervals = []
    for extra in args.extra_mask:
        rows = read_bed_intervals(Path(extra))
        extra_intervals.extend(rows)
        log(f"read {len(rows)} extra-mask interval(s) from {extra}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for complement in complements:
        suffix = "fa.gz" if args.gzip else "fa"
        out_fasta = outdir / f"{default_prefix}.{complement}.{suffix}"
        out_sizes = outdir / f"{default_prefix}.{complement}.chrom.sizes"
        out_bed = outdir / f"{default_prefix}.{complement}.masked.bed"
        out_meta = outdir / f"{default_prefix}.{complement}.metadata.json"
        out_md5 = outdir / f"{default_prefix}.{complement}.fa.md5"

        if out_fasta.exists() and not args.force:
            parser.error(f"{out_fasta} exists; pass --force to overwrite.")

        log(f"=== {complement}: streaming + masking -> {out_fasta.name} ===")
        mask_plan = build_mask_plan(
            par_intervals=par_intervals,
            extra_intervals=extra_intervals,
            x_contig=x_contig,
            y_contig=y_contig,
            complement=complement,
            y_length_hint=None,
        )
        lengths, applied = mask_fasta_stream(
            source_fasta=fasta_path,
            out_fasta=out_fasta,
            mask_plan=mask_plan,
            allowlist=allowlist,
        )
        write_chrom_sizes(out_sizes, lengths)
        log(f"  wrote {out_sizes.name} ({len(lengths)} contigs)")
        write_masked_bed(out_bed, applied)
        log(f"  wrote {out_bed.name} ({len(applied)} masked interval(s))")
        log("  computing md5 checksum ...")
        checksum = md5_of_file(out_fasta)
        out_md5.write_text(f"{checksum}  {out_fasta.name}\n")

        metadata = {
            "assembly": args.assembly,
            "complement": complement,
            "source_fasta": str(fasta_path),
            "par_bed": str(par_bed_path) if par_bed_path else None,
            "par_status": par_status,
            "extra_mask_beds": list(args.extra_mask),
            "x_contig": x_contig,
            "y_contig": y_contig,
            "n_masked_intervals": len(applied),
            "masked_bp": sum(end - start for _c, start, end, _n in applied),
            "output_fasta": out_fasta.name,
            "output_fasta_md5": checksum,
            "method": (
                "XX: whole-Y hard-mask. XY: Y-PAR hard-mask (X intact). "
                "Standard follows Olney 2020 & AJHG 2025 best practices."
            ),
        }
        out_meta.write_text(json.dumps(metadata, indent=2) + "\n")

        print(f"[{complement}] {out_fasta}  (masked {metadata['masked_bp']:,} bp "
              f"in {len(applied)} interval(s); md5 {checksum})", file=sys.stderr)


if __name__ == "__main__":
    main()
