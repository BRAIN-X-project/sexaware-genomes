#!/usr/bin/env python
"""Build a unified chromosome/region territory table for the primary human assembly.

The wide output has ``chrom``, ``region`` and ``region_group``. Autosomes are
reported as whole-chromosome regions, while X/Y are kept both as complete
chromosome rows and, when supplied, as rows from the ``xy_regions`` partition.
The table quantifies how much sequence is usable at
several increasingly strict definitions:

- ``raw_bp``: total chromosome length.
- ``nongap_bp``: length minus assembly gaps.
- ``structural_callable_bp``: length minus gaps and structural cytobands
  (centromere, variable heterochromatin, acrocentric stalk).
- mappable base pairs at Umap k100 score thresholds (0.50/0.70/0.90).
- the same, additionally excluding an ENCODE-style blacklist.

Inputs are standard UCSC-style tab files plus a Umap mappability bigWig.
When a Y-PAR-masked Umap bigWig is used, PAR1Y/PAR2Y mappability-derived
values can be filled from PAR1X/PAR2X while preserving Y raw/nongap/structural
lengths.
"""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

PRIMARY_CHROMS = [str(index) for index in range(1, 23)] + ["X", "Y"]
DEFAULT_STRUCTURAL_EXCLUDE_STAINS = ("acen", "gvar", "stalk")
PAR_POLICY_CHOICES = (
    "reference_hg38_duplicated",
    "y_par_masked",
    "par_collapsed",
    "empirical_bam_callable",
)

TERRITORY_BP_COLUMNS = [
    "raw_bp",
    "nongap_bp",
    "structural_callable_bp",
    "mappable_k100_t50_bp",
    "mappable_k100_t70_bp",
    "mappable_k100_t90_bp",
    "effective_mappable_k100_bp",
    "callable_k100_t90_no_blacklist_bp",
    "effective_callable_k100_no_blacklist_bp",
]
MAPPABILITY_BP_COLUMNS = [
    "mappable_k100_t50_bp",
    "mappable_k100_t70_bp",
    "mappable_k100_t90_bp",
    "effective_mappable_k100_bp",
    "callable_k100_t90_no_blacklist_bp",
    "effective_callable_k100_no_blacklist_bp",
]
UNIFIED_ID_COLUMNS = ["chrom", "region", "region_group"]


def cast_numeric_territory_columns_to_int64(
    territory: pd.DataFrame,
    columns: list[str] | tuple[str, ...] | None = None,
) -> pd.DataFrame:
    """Return ``territory`` with all territory metric columns as int64.

    Umap ``effective_*`` values are score-weighted sums and can arrive as
    floats. For this workflow they are written as base-pair-like denominators,
    so they are rounded to the nearest integer before casting to ``int64``.
    """
    out = territory.copy()
    if columns is None:
        columns = [column for column in TERRITORY_BP_COLUMNS if column in out.columns]
    for column in columns:
        if column not in out.columns:
            continue
        values = pd.to_numeric(out[column], errors="coerce").fillna(0)
        out[column] = np.rint(values).astype("int64")
    return out


def _parse_comma_list(value: str | None) -> tuple[str, ...]:
    """Parse a comma-separated CLI option into a tuple of non-empty strings."""
    if value is None:
        return DEFAULT_STRUCTURAL_EXCLUDE_STAINS
    items = tuple(item.strip() for item in value.split(",") if item.strip())
    return items


def _threshold_label(threshold: float) -> str:
    """Return the integer-percent column suffix for a score threshold.

    ``round`` is used instead of ``int`` because ``int(0.58 * 100)`` truncates
    to ``57`` due to floating-point representation.
    """
    return f"t{round(threshold * 100):02d}_bp"


def normalize_chrom_name(chrom: pd.Series) -> pd.Series:
    """Convert UCSC ``chr1``/``chrX`` names to Ensembl-like ``1``/``X`` names.

    Parameters
    ----------
    chrom : pandas.Series
        Chromosome names in either UCSC or Ensembl style.

    Returns
    -------
    pandas.Series
        Names with a leading ``chr`` prefix removed.
    """
    return chrom.astype(str).str.replace(r"^chr", "", regex=True)


def read_chrom_info(path: str | Path) -> pd.DataFrame:
    """Read a UCSC ``chromInfo``/``chrom.sizes`` file for primary chromosomes.

    Parameters
    ----------
    path : str or pathlib.Path
        Tab-separated file whose first two columns are chromosome name and
        length in base pairs.

    Returns
    -------
    pandas.DataFrame
        Columns ``chrom`` (ordered categorical over :data:`PRIMARY_CHROMS`)
        and ``raw_bp`` (int64), restricted to the primary chromosomes and
        sorted in canonical order.
    """
    chrom_info = pd.read_csv(
        path,
        sep="\t",
        header=None,
        usecols=[0, 1],
        names=["chrom", "raw_bp"],
    )
    chrom_info["chrom"] = normalize_chrom_name(chrom_info["chrom"])
    chrom_info = chrom_info[chrom_info["chrom"].isin(PRIMARY_CHROMS)].copy()
    chrom_info["raw_bp"] = chrom_info["raw_bp"].astype("int64")
    chrom_info["chrom"] = pd.Categorical(
        chrom_info["chrom"],
        categories=PRIMARY_CHROMS,
        ordered=True,
    )
    return chrom_info.sort_values("chrom").reset_index(drop=True)


def read_gap(path: str | Path) -> pd.DataFrame:
    """Read a UCSC ``gap`` table, keeping only the interval and gap type.

    Parameters
    ----------
    path : str or pathlib.Path
        Tab-separated UCSC gap table with columns
        ``bin, chrom, chromStart, chromEnd, ix, n, size, type, bridge``.

    Returns
    -------
    pandas.DataFrame
        Columns ``chrom``, ``start``, ``end``, ``type`` with normalized
        chromosome names. Coordinates are 0-based half-open.
    """
    gap_columns = [
        "bin",
        "chrom",
        "start",
        "end",
        "ix",
        "n",
        "size",
        "type",
        "bridge",
    ]
    gap = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=gap_columns,
        usecols=["chrom", "start", "end", "type"],
    )
    gap["chrom"] = normalize_chrom_name(gap["chrom"])
    return gap[["chrom", "start", "end", "type"]]


def read_cytoband(path: str | Path) -> pd.DataFrame:
    """Read a UCSC ``cytoBand`` table.

    Parameters
    ----------
    path : str or pathlib.Path
        Tab-separated cytoband table with columns
        ``chrom, chromStart, chromEnd, name, gieStain``.

    Returns
    -------
    pandas.DataFrame
        Columns ``chrom``, ``start``, ``end``, ``name``, ``gieStain`` with
        normalized chromosome names. Coordinates are 0-based half-open.
    """
    cytoband = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name", "gieStain"],
    )
    cytoband["chrom"] = normalize_chrom_name(cytoband["chrom"])
    return cytoband[["chrom", "start", "end", "name", "gieStain"]]


def read_bed3(path: str | Path) -> pd.DataFrame:
    """Read the first three columns of a BED file as an interval table.

    Parameters
    ----------
    path : str or pathlib.Path
        BED file; ``#``-prefixed lines are treated as comments.

    Returns
    -------
    pandas.DataFrame
        Columns ``chrom``, ``start``, ``end`` with normalized chromosome
        names. Coordinates are 0-based half-open.
    """
    bed = pd.read_csv(
        path,
        sep="\t",
        header=None,
        comment="#",
        usecols=[0, 1, 2],
        names=["chrom", "start", "end"],
    )
    bed["chrom"] = normalize_chrom_name(bed["chrom"])
    return bed[["chrom", "start", "end"]]


def clip_intervals_to_chroms(
    intervals: pd.DataFrame,
    chrom_info: pd.DataFrame,
) -> pd.DataFrame:
    """Restrict intervals to known chromosomes and their bounds.

    Intervals on chromosomes absent from ``chrom_info`` are dropped. Remaining
    intervals are clamped to ``[0, chromosome_length)`` and any that become
    empty are removed.

    Parameters
    ----------
    intervals : pandas.DataFrame
        Must contain ``chrom``, ``start``, ``end`` (0-based half-open).
    chrom_info : pandas.DataFrame
        Output of :func:`read_chrom_info`, providing chromosome lengths.

    Returns
    -------
    pandas.DataFrame
        Columns ``chrom``, ``start``, ``end`` with non-empty, in-bounds
        intervals.
    """
    chrom_lengths = chrom_info.set_index("chrom")["raw_bp"].astype("int64").to_dict()
    clipped = intervals.copy()
    clipped = clipped[clipped["chrom"].isin(chrom_lengths)].copy()
    if clipped.empty:
        return pd.DataFrame(columns=["chrom", "start", "end"])
    clipped["start"] = clipped["start"].astype("int64").clip(lower=0)
    clipped["end"] = clipped["end"].astype("int64")
    length_per_row = clipped["chrom"].map(chrom_lengths).astype("int64")
    clipped["end"] = np.minimum(clipped["end"], length_per_row)
    clipped = clipped[clipped["end"] > clipped["start"]].copy()
    return clipped[["chrom", "start", "end"]]


def merge_intervals(intervals: pd.DataFrame) -> pd.DataFrame:
    """Merge overlapping or adjacent BED-like intervals per chromosome.

    Uses the standard sort-and-sweep approach, vectorized: after sorting by
    ``(chrom, start, end)``, a new merged block starts whenever the row's
    ``start`` lies beyond the running maximum ``end`` of the preceding rows on
    the same chromosome, or the chromosome changes. Adjacent intervals that
    touch at a boundary are merged (matching half-open semantics).

    Parameters
    ----------
    intervals : pandas.DataFrame
        Must contain ``chrom``, ``start``, ``end`` (0-based half-open).

    Returns
    -------
    pandas.DataFrame
        Columns ``chrom``, ``start``, ``end`` with disjoint, sorted intervals.
    """
    if intervals.empty:
        return pd.DataFrame(columns=["chrom", "start", "end"])
    sorted_intervals = intervals[["chrom", "start", "end"]].copy()
    sorted_intervals["start"] = sorted_intervals["start"].astype("int64")
    sorted_intervals["end"] = sorted_intervals["end"].astype("int64")
    sorted_intervals = sorted_intervals.sort_values(["chrom", "start", "end"]).reset_index(
        drop=True
    )

    # Running maximum end of the rows seen so far within each chromosome,
    # offset by one so it refers to the *preceding* rows only.
    running_max_end = sorted_intervals.groupby("chrom", observed=True)["end"].cummax().shift()
    same_chrom_as_previous = sorted_intervals["chrom"].eq(sorted_intervals["chrom"].shift())
    starts_new_block = (~same_chrom_as_previous) | (sorted_intervals["start"] > running_max_end)
    block_id = starts_new_block.cumsum()

    merged = sorted_intervals.groupby(block_id, sort=False).agg(
        chrom=("chrom", "first"),
        start=("start", "min"),
        end=("end", "max"),
    )
    return merged.reset_index(drop=True)


def interval_lengths_by_chrom(intervals: pd.DataFrame) -> pd.Series:
    """Sum interval lengths per chromosome.

    Parameters
    ----------
    intervals : pandas.DataFrame
        Must contain ``chrom``, ``start``, ``end`` (0-based half-open).

    Returns
    -------
    pandas.Series
        Total base pairs per chromosome, indexed by ``chrom`` (int64).
    """
    if intervals.empty:
        return pd.Series(dtype="int64")
    length_per_row = intervals["end"].astype("int64") - intervals["start"].astype("int64")
    return length_per_row.groupby(intervals["chrom"]).sum().astype("int64")


def complement_intervals(
    chrom_info: pd.DataFrame,
    excluded_intervals: pd.DataFrame,
) -> pd.DataFrame:
    """Return the chromosome regions not covered by the excluded intervals.

    Parameters
    ----------
    chrom_info : pandas.DataFrame
        Output of :func:`read_chrom_info`, providing chromosome lengths.
    excluded_intervals : pandas.DataFrame
        Regions to subtract, with ``chrom``, ``start``, ``end``. They are
        internally clipped and merged, so they need not be sorted or disjoint.

    Returns
    -------
    pandas.DataFrame
        Columns ``chrom``, ``start``, ``end`` covering everything on each
        chromosome that the excluded intervals do not.
    """
    excluded = merge_intervals(clip_intervals_to_chroms(excluded_intervals, chrom_info))
    excluded_by_chrom = {
        chrom: block[["start", "end"]].to_numpy(dtype=np.int64)
        for chrom, block in excluded.groupby("chrom", observed=True)
    }
    complement_rows = []
    for chrom_value, chrom_length in chrom_info[["chrom", "raw_bp"]].itertuples(
        index=False,
        name=None,
    ):
        chrom = str(chrom_value)
        length = int(chrom_length)
        cursor = 0
        for start, end in excluded_by_chrom.get(chrom, []):
            start = int(start)
            end = int(end)
            if start > cursor:
                complement_rows.append((chrom, cursor, start))
            cursor = max(cursor, end)
        if cursor < length:
            complement_rows.append((chrom, cursor, length))
    return pd.DataFrame(complement_rows, columns=["chrom", "start", "end"])


def bigwig_name_for_chrom(chrom: str, bigwig_chroms: set[str]) -> str | None:
    """Match an Ensembl-style name to the bigWig's own naming.

    Parameters
    ----------
    chrom : str
        Ensembl-style chromosome name (``1``, ``X`` ...).
    bigwig_chroms : set of str
        Chromosome names present in the bigWig.

    Returns
    -------
    str or None
        The matching bigWig chromosome name (``1`` or ``chr1``), or ``None``
        if neither form is present.
    """
    if chrom in bigwig_chroms:
        return chrom
    ucsc_name = f"chr{chrom}"
    if ucsc_name in bigwig_chroms:
        return ucsc_name
    return None


def mappability_sums_by_chrom(
    intervals: pd.DataFrame,
    bigwig_path: str | Path,
    thresholds: tuple[float, ...] = (0.50, 0.70, 0.90),
    chunk_size: int = 1_000_000,
) -> pd.DataFrame:
    """Summarize mappability scores over the given intervals, per chromosome.

    For each chromosome this computes the effective (score-weighted) base
    pairs, ``sum(score)`` over every base, and the count of bases whose score
    meets each threshold. Missing bigWig values are treated as 0. Intervals are
    additionally clipped to the bigWig's own chromosome lengths so that
    coordinate systems that disagree with ``chrom_info`` do not raise.

    Parameters
    ----------
    intervals : pandas.DataFrame
        Regions to scan, with ``chrom``, ``start``, ``end`` (0-based
        half-open).
    bigwig_path : str or pathlib.Path
        Path to the Umap mappability bigWig.
    thresholds : tuple of float, optional
        Score cutoffs for the ``>=`` base counts. Defaults to
        ``(0.50, 0.70, 0.90)``.
    chunk_size : int, optional
        Number of base pairs fetched per bigWig read. Defaults to 1,000,000.

    Returns
    -------
    pandas.DataFrame
        One row per primary chromosome with ``chrom``, ``effective_bp`` and a
        ``t{pct}_bp`` column per threshold.
    """
    import pyBigWig

    bigwig = pyBigWig.open(str(bigwig_path))
    bigwig_chroms = bigwig.chroms()
    bigwig_chrom_names = set(bigwig_chroms.keys())
    records = {
        chrom: {
            "effective_bp": 0.0,
            **{_threshold_label(threshold): 0 for threshold in thresholds},
        }
        for chrom in PRIMARY_CHROMS
    }
    missing_chroms: set[str] = set()

    for chrom_value, start_value, end_value in intervals[["chrom", "start", "end"]].itertuples(
        index=False, name=None
    ):
        chrom = str(chrom_value)
        bigwig_chrom = bigwig_name_for_chrom(chrom, bigwig_chrom_names)
        if bigwig_chrom is None:
            missing_chroms.add(chrom)
            continue
        start = int(start_value)
        # Never read past the length the bigWig itself reports for this chrom.
        end = min(int(end_value), bigwig_chroms[bigwig_chrom])
        for chunk_start in range(start, end, chunk_size):
            chunk_end = min(chunk_start + chunk_size, end)
            values = bigwig.values(bigwig_chrom, chunk_start, chunk_end, numpy=True)
            if values is None:
                continue
            values = np.asarray(values, dtype=np.float32)
            values = np.nan_to_num(values, nan=0.0)
            records[chrom]["effective_bp"] += float(values.sum(dtype=np.float64))
            for threshold in thresholds:
                records[chrom][_threshold_label(threshold)] += int(
                    np.count_nonzero(values >= threshold)
                )
    bigwig.close()

    if missing_chroms:
        warnings.warn(
            "No bigWig data for chromosome(s): "
            f"{', '.join(sorted(missing_chroms))}; treated as fully unmappable.",
            stacklevel=2,
        )

    summary = pd.DataFrame([{"chrom": chrom, **values} for chrom, values in records.items()])
    summary["chrom"] = pd.Categorical(
        summary["chrom"],
        categories=PRIMARY_CHROMS,
        ordered=True,
    )
    return summary.sort_values("chrom").reset_index(drop=True)


def build_chrom_territory(
    chrom_info_path: str | Path,
    gap_path: str | Path,
    cytoband_path: str | Path,
    umap_k100_bw_path: str | Path,
    blacklist_path: str | Path,
    structural_exclude_stains: tuple[str, ...] = DEFAULT_STRUCTURAL_EXCLUDE_STAINS,
    par_policy: str = "reference_hg38_duplicated",
    validate: bool = True,
) -> pd.DataFrame:
    """Assemble the full per-chromosome territory table.

    Combines assembly gaps, structural cytobands, a blacklist and Umap k100
    mappability into a single table describing usable sequence per chromosome
    at several strictness levels.

    Parameters
    ----------
    chrom_info_path : str or pathlib.Path
        UCSC ``chromInfo``/``chrom.sizes`` file.
    gap_path : str or pathlib.Path
        UCSC ``gap`` table.
    cytoband_path : str or pathlib.Path
        UCSC ``cytoBand`` table.
    umap_k100_bw_path : str or pathlib.Path
        Umap k100 mappability bigWig.
    blacklist_path : str or pathlib.Path
        ENCODE-style blacklist BED.
    structural_exclude_stains : tuple of str, optional
        ``cytoBand.gieStain`` values to subtract together with gaps before
        computing ``structural_callable_bp``. Defaults to ``("acen", "gvar",
        "stalk")``. Passing only ``("acen",)`` keeps variable heterochromatin
        and acrocentric stalks in the structural territory.
    par_policy : str, optional
        Metadata describing how pseudoautosomal regions should be interpreted.
        The calculation itself is driven by the supplied bigWig; this value is
        stored in ``DataFrame.attrs`` and in the optional JSON metadata output.
    validate : bool, optional
        If True, run consistency checks after the table is built.

    Returns
    -------
    pandas.DataFrame
        One row per primary chromosome with length, gap-corrected length,
        structural-callable length, mappable base pairs at each threshold and
        the blacklist-subtracted variants.
    """
    chrom_info = read_chrom_info(chrom_info_path)

    gap = read_gap(gap_path)
    gap_intervals = clip_intervals_to_chroms(gap[["chrom", "start", "end"]], chrom_info)
    gap_intervals = merge_intervals(gap_intervals)

    cytoband = read_cytoband(cytoband_path)
    cytoband = cytoband[cytoband["chrom"].isin(PRIMARY_CHROMS)].copy()

    # Structural cytoband exclusions:
    #   acen  = centromere
    #   gvar  = variable heterochromatin
    #   stalk = acrocentric stalk
    structural_cytoband = cytoband[cytoband["gieStain"].isin(structural_exclude_stains)][
        ["chrom", "start", "end"]
    ]
    structural_exclusions = pd.concat(
        [gap_intervals, clip_intervals_to_chroms(structural_cytoband, chrom_info)],
        ignore_index=True,
    )
    structural_exclusions = merge_intervals(
        clip_intervals_to_chroms(structural_exclusions, chrom_info)
    )
    structural_callable_intervals = complement_intervals(
        chrom_info,
        structural_exclusions,
    )

    blacklist = read_bed3(blacklist_path)
    blacklist = merge_intervals(clip_intervals_to_chroms(blacklist, chrom_info))
    structural_plus_blacklist_exclusions = merge_intervals(
        pd.concat([structural_exclusions, blacklist], ignore_index=True)
    )
    structural_no_blacklist_intervals = complement_intervals(
        chrom_info,
        structural_plus_blacklist_exclusions,
    )

    # Length columns.
    territory = chrom_info.copy()
    gap_bp = interval_lengths_by_chrom(gap_intervals)
    territory["gap_bp"] = territory["chrom"].astype(str).map(gap_bp).fillna(0).astype("int64")
    territory["nongap_bp"] = territory["raw_bp"] - territory["gap_bp"]
    structural_callable_bp = interval_lengths_by_chrom(structural_callable_intervals)
    territory["structural_callable_bp"] = (
        territory["chrom"].astype(str).map(structural_callable_bp).fillna(0).astype("int64")
    )

    # Mappability over the structural-callable territory.
    mappability_structural = mappability_sums_by_chrom(
        structural_callable_intervals,
        umap_k100_bw_path,
        thresholds=(0.50, 0.70, 0.90),
    ).rename(
        columns={
            "t50_bp": "mappable_k100_t50_bp",
            "t70_bp": "mappable_k100_t70_bp",
            "t90_bp": "mappable_k100_t90_bp",
            "effective_bp": "effective_mappable_k100_bp",
        }
    )

    # Mappability over the structural-callable territory minus the blacklist.
    mappability_no_blacklist = mappability_sums_by_chrom(
        structural_no_blacklist_intervals,
        umap_k100_bw_path,
        thresholds=(0.90,),
    ).rename(
        columns={
            "t90_bp": "callable_k100_t90_no_blacklist_bp",
            "effective_bp": "effective_callable_k100_no_blacklist_bp",
        }
    )

    territory = territory.merge(mappability_structural, on="chrom", how="left")
    territory = territory.merge(
        mappability_no_blacklist[
            [
                "chrom",
                "callable_k100_t90_no_blacklist_bp",
                "effective_callable_k100_no_blacklist_bp",
            ]
        ],
        on="chrom",
        how="left",
    )

    wanted_columns = [
        "chrom",
        "raw_bp",
        "nongap_bp",
        "structural_callable_bp",
        "mappable_k100_t50_bp",
        "mappable_k100_t70_bp",
        "mappable_k100_t90_bp",
        "effective_mappable_k100_bp",
        "callable_k100_t90_no_blacklist_bp",
        "effective_callable_k100_no_blacklist_bp",
    ]
    territory = territory[wanted_columns].copy()

    territory = cast_numeric_territory_columns_to_int64(territory)

    territory["chrom"] = pd.Categorical(
        territory["chrom"].astype(str),
        categories=PRIMARY_CHROMS,
        ordered=True,
    )
    territory = territory.sort_values("chrom").reset_index(drop=True)
    territory.attrs["par_policy"] = par_policy
    territory.attrs["structural_exclude_stains"] = tuple(structural_exclude_stains)
    if validate:
        validate_territory_table(territory, label_columns=("chrom",))
    return territory


def _require_existing_file(path: str | Path, flag: str) -> Path:
    """Return ``path`` as a Path, raising if it does not point to a file."""
    resolved = Path(path)
    if not resolved.is_file():
        raise FileNotFoundError(f"{flag}: file not found: {resolved}")
    return resolved


def read_xy_regions(path: str | Path) -> pd.DataFrame:
    """Read the X/Y sub-chromosomal region partition.

    The file is 1-based inclusive (matching genome browsers); it is converted to the
    0-based half-open convention used everywhere else in this script.

    Parameters
    ----------
    path : str or pathlib.Path
        Tab-separated file with columns ``chrom, region, region_group, start, end``.

    Returns
    -------
    pandas.DataFrame
        Columns ``chrom``, ``region``, ``region_group``, ``start``, ``end`` (0-based
        half-open), with normalized chromosome names.
    """
    regions = pd.read_csv(path, sep="\t", comment="#")
    regions["chrom"] = normalize_chrom_name(regions["chrom"])
    regions["start"] = regions["start"].astype("int64") - 1
    regions["end"] = regions["end"].astype("int64")
    return regions[["chrom", "region", "region_group", "start", "end"]]


def interval_lengths_by_label(intervals: pd.DataFrame, label_col: str) -> pd.Series:
    """Sum interval lengths grouped by an arbitrary label column.

    Parameters
    ----------
    intervals : pandas.DataFrame
        Must contain ``start``, ``end`` (0-based half-open) and ``label_col``.
    label_col : str
        Column to group by (e.g. ``"region"``).

    Returns
    -------
    pandas.Series
        Total base pairs per label (int64).
    """
    if intervals.empty:
        return pd.Series(dtype="int64")
    length_per_row = intervals["end"].astype("int64") - intervals["start"].astype("int64")
    return length_per_row.groupby(intervals[label_col]).sum().astype("int64")


def intersect_with_regions(
    intervals: pd.DataFrame,
    region_intervals: pd.DataFrame,
) -> pd.DataFrame:
    """Intersect interval blocks with region intervals, tagging each piece with its region.

    Both inputs are 0-based half-open with a ``chrom`` column; ``region_intervals`` also
    carries ``region`` and ``region_group``. Every overlap between an input block and a
    region interval on the same chromosome is emitted as its own piece.

    Parameters
    ----------
    intervals : pandas.DataFrame
        Blocks to split, with ``chrom``, ``start``, ``end``.
    region_intervals : pandas.DataFrame
        Region partition, with ``chrom``, ``region``, ``region_group``, ``start``, ``end``.

    Returns
    -------
    pandas.DataFrame
        Columns ``chrom``, ``start``, ``end``, ``region``, ``region_group``.
    """
    empty = pd.DataFrame(columns=["chrom", "start", "end", "region", "region_group"])
    if intervals.empty or region_intervals.empty:
        return empty
    regions_by_chrom = {
        chrom: block for chrom, block in region_intervals.groupby("chrom", observed=True)
    }
    pieces = []
    for chrom, block in intervals.groupby("chrom", observed=True):
        region_block = regions_by_chrom.get(chrom)
        if region_block is None:
            continue
        block_bounds = block[["start", "end"]].to_numpy(dtype=np.int64)
        for region_start, region_end, region, region_group in region_block[
            ["start", "end", "region", "region_group"]
        ].itertuples(index=False, name=None):
            lower = np.maximum(block_bounds[:, 0], int(region_start))
            upper = np.minimum(block_bounds[:, 1], int(region_end))
            keep = upper > lower
            if not keep.any():
                continue
            pieces.append(
                pd.DataFrame(
                    {
                        "chrom": chrom,
                        "start": lower[keep],
                        "end": upper[keep],
                        "region": region,
                        "region_group": region_group,
                    }
                )
            )
    if not pieces:
        return empty
    return pd.concat(pieces, ignore_index=True)


def mappability_sums_by_label(
    intervals: pd.DataFrame,
    bigwig_path: str | Path,
    labels: list[str],
    label_col: str = "region",
    thresholds: tuple[float, ...] = (0.50, 0.70, 0.90),
    chunk_size: int = 1_000_000,
) -> pd.DataFrame:
    """Summarize mappability over intervals grouped by an arbitrary label.

    Mirrors :func:`mappability_sums_by_chrom`, but the bigWig is still looked up by the
    ``chrom`` column while sums are accumulated per ``label_col`` value.

    Parameters
    ----------
    intervals : pandas.DataFrame
        Regions to scan, with ``chrom``, ``start``, ``end`` and ``label_col``.
    bigwig_path : str or pathlib.Path
        Path to the Umap mappability bigWig.
    labels : list of str
        The label values to report (rows with no intervals report zeros).
    label_col : str, default "region"
        Column whose values group the sums.
    thresholds : tuple of float, optional
        Score cutoffs for the ``>=`` base counts.
    chunk_size : int, optional
        Base pairs fetched per bigWig read.

    Returns
    -------
    pandas.DataFrame
        One row per label with ``label_col``, ``effective_bp`` and a ``t{pct}_bp`` column
        per threshold.
    """
    import pyBigWig

    bigwig = pyBigWig.open(str(bigwig_path))
    bigwig_chroms = bigwig.chroms()
    bigwig_chrom_names = set(bigwig_chroms.keys())
    records = {
        label: {
            "effective_bp": 0.0,
            **{_threshold_label(threshold): 0 for threshold in thresholds},
        }
        for label in labels
    }
    missing_chroms: set[str] = set()

    for chrom_value, start_value, end_value, label in intervals[
        ["chrom", "start", "end", label_col]
    ].itertuples(index=False, name=None):
        chrom = str(chrom_value)
        bigwig_chrom = bigwig_name_for_chrom(chrom, bigwig_chrom_names)
        if bigwig_chrom is None:
            missing_chroms.add(chrom)
            continue
        start = int(start_value)
        end = min(int(end_value), bigwig_chroms[bigwig_chrom])
        for chunk_start in range(start, end, chunk_size):
            chunk_end = min(chunk_start + chunk_size, end)
            values = bigwig.values(bigwig_chrom, chunk_start, chunk_end, numpy=True)
            if values is None:
                continue
            values = np.nan_to_num(np.asarray(values, dtype=np.float32), nan=0.0)
            records[label]["effective_bp"] += float(values.sum(dtype=np.float64))
            for threshold in thresholds:
                records[label][_threshold_label(threshold)] += int(
                    np.count_nonzero(values >= threshold)
                )
    bigwig.close()

    if missing_chroms:
        warnings.warn(
            "No bigWig data for chromosome(s): "
            f"{', '.join(sorted(missing_chroms))}; treated as fully unmappable.",
            stacklevel=2,
        )

    return pd.DataFrame([{label_col: label, **values} for label, values in records.items()])


def build_region_territory(
    chrom_info_path: str | Path,
    gap_path: str | Path,
    cytoband_path: str | Path,
    umap_k100_bw_path: str | Path,
    blacklist_path: str | Path,
    xy_regions_path: str | Path,
    structural_exclude_stains: tuple[str, ...] = DEFAULT_STRUCTURAL_EXCLUDE_STAINS,
    par_policy: str = "reference_hg38_duplicated",
    validate: bool = True,
) -> pd.DataFrame:
    """Assemble a per-region territory table for the X/Y sub-chromosomal partition.

    Mirrors :func:`build_chrom_territory`, but every base-pair column is keyed by the fine
    ``region`` (carrying its coarse ``region_group``): each callable interval set is
    intersected with the region intervals and re-summed, and mappability is scanned over
    those intersections. Because the regions are a disjoint partition of X and Y, any coarser
    level (e.g. ``NPX``) is a groupby-sum of these rows, and each chromosome's per-region
    ``raw_bp`` sums back to its whole-chromosome length.

    Parameters
    ----------
    chrom_info_path, gap_path, cytoband_path, umap_k100_bw_path, blacklist_path
        Same inputs as :func:`build_chrom_territory`.
    xy_regions_path : str or pathlib.Path
        The ``xy_regions.txt`` partition file.

    Returns
    -------
    pandas.DataFrame
        One row per fine region with ``chrom``, ``region``, ``region_group`` and the same
        base-pair columns as the per-chromosome table.
    """
    chrom_info = read_chrom_info(chrom_info_path)
    regions = read_xy_regions(xy_regions_path)
    regions = regions[regions["chrom"].isin(chrom_info["chrom"].astype(str))].copy()
    region_labels = regions["region"].drop_duplicates().tolist()
    region_group_of = dict(zip(regions["region"], regions["region_group"]))
    region_chrom_of = dict(zip(regions["region"], regions["chrom"]))
    region_intervals = regions[["chrom", "region", "region_group", "start", "end"]]

    gap_intervals = merge_intervals(
        clip_intervals_to_chroms(read_gap(gap_path)[["chrom", "start", "end"]], chrom_info)
    )
    nongap_intervals = complement_intervals(chrom_info, gap_intervals)

    cytoband = read_cytoband(cytoband_path)
    cytoband = cytoband[cytoband["chrom"].isin(PRIMARY_CHROMS)].copy()
    structural_cytoband = cytoband[cytoband["gieStain"].isin(structural_exclude_stains)][
        ["chrom", "start", "end"]
    ]
    structural_exclusions = merge_intervals(
        clip_intervals_to_chroms(
            pd.concat(
                [gap_intervals, clip_intervals_to_chroms(structural_cytoband, chrom_info)],
                ignore_index=True,
            ),
            chrom_info,
        )
    )
    structural_callable_intervals = complement_intervals(chrom_info, structural_exclusions)

    blacklist = merge_intervals(clip_intervals_to_chroms(read_bed3(blacklist_path), chrom_info))
    structural_no_blacklist_intervals = complement_intervals(
        chrom_info,
        merge_intervals(pd.concat([structural_exclusions, blacklist], ignore_index=True)),
    )

    raw_by_region = interval_lengths_by_label(region_intervals, "region")
    nongap_by_region = interval_lengths_by_label(
        intersect_with_regions(nongap_intervals, region_intervals), "region"
    )
    structural_callable_by_region = interval_lengths_by_label(
        intersect_with_regions(structural_callable_intervals, region_intervals), "region"
    )

    mappability_structural = mappability_sums_by_label(
        intersect_with_regions(structural_callable_intervals, region_intervals),
        umap_k100_bw_path,
        region_labels,
        thresholds=(0.50, 0.70, 0.90),
    ).rename(
        columns={
            "t50_bp": "mappable_k100_t50_bp",
            "t70_bp": "mappable_k100_t70_bp",
            "t90_bp": "mappable_k100_t90_bp",
            "effective_bp": "effective_mappable_k100_bp",
        }
    )
    mappability_no_blacklist = mappability_sums_by_label(
        intersect_with_regions(structural_no_blacklist_intervals, region_intervals),
        umap_k100_bw_path,
        region_labels,
        thresholds=(0.90,),
    ).rename(
        columns={
            "t90_bp": "callable_k100_t90_no_blacklist_bp",
            "effective_bp": "effective_callable_k100_no_blacklist_bp",
        }
    )

    territory = pd.DataFrame({"region": region_labels})
    territory["chrom"] = territory["region"].map(region_chrom_of)
    territory["region_group"] = territory["region"].map(region_group_of)
    territory["raw_bp"] = territory["region"].map(raw_by_region).fillna(0).astype("int64")
    territory["nongap_bp"] = territory["region"].map(nongap_by_region).fillna(0).astype("int64")
    territory["structural_callable_bp"] = (
        territory["region"].map(structural_callable_by_region).fillna(0).astype("int64")
    )
    territory = territory.merge(mappability_structural, on="region", how="left")
    territory = territory.merge(
        mappability_no_blacklist[
            [
                "region",
                "callable_k100_t90_no_blacklist_bp",
                "effective_callable_k100_no_blacklist_bp",
            ]
        ],
        on="region",
        how="left",
    )

    column_order = [
        "chrom",
        "region",
        "region_group",
        "raw_bp",
        "nongap_bp",
        "structural_callable_bp",
        "mappable_k100_t50_bp",
        "mappable_k100_t70_bp",
        "mappable_k100_t90_bp",
        "effective_mappable_k100_bp",
        "callable_k100_t90_no_blacklist_bp",
        "effective_callable_k100_no_blacklist_bp",
    ]
    territory = territory[column_order]

    territory = cast_numeric_territory_columns_to_int64(territory)

    territory.attrs["par_policy"] = par_policy
    territory.attrs["structural_exclude_stains"] = tuple(structural_exclude_stains)
    if validate:
        validate_territory_table(territory, label_columns=("chrom", "region"))
    return territory


def validate_territory_table(
    territory: pd.DataFrame,
    label_columns: tuple[str, ...] = ("chrom",),
    float_tolerance: float = 1e-6,
) -> None:
    """Validate monotonic and upper-bound relationships between territory columns.

    These checks catch common failure modes: manual PAR extrapolations that exceed the
    callable territory, threshold columns that are not monotonic, and negative values.
    They intentionally warn rather than raise so long production jobs still write their
    tables, but any warning should be inspected before using the output as a denominator.
    """
    if territory.empty:
        return

    id_frame = territory[list(label_columns)].astype(str)
    labels = id_frame.agg(":".join, axis=1)

    numeric_columns = [
        column for column in territory.columns
        if column.endswith("_bp") or column.startswith("effective")
    ]
    for column in numeric_columns:
        if column not in territory.columns:
            continue
        bad = territory[column].astype(float) < -float_tolerance
        for label in labels[bad]:
            warnings.warn(f"Negative value in {column} for {label}.", stacklevel=2)

    chain = [
        "mappable_k100_t90_bp",
        "mappable_k100_t70_bp",
        "mappable_k100_t50_bp",
        "structural_callable_bp",
    ]
    for lower, upper in zip(chain, chain[1:]):
        if lower in territory.columns and upper in territory.columns:
            bad = territory[lower].astype(float) > territory[upper].astype(float) + float_tolerance
            for label, observed, expected in zip(labels[bad], territory.loc[bad, lower], territory.loc[bad, upper]):
                warnings.warn(
                    f"Unexpected territory ordering for {label}: {lower}={observed} > {upper}={expected}.",
                    stacklevel=2,
                )

    upper_bound_checks = [
        ("effective_mappable_k100_bp", "structural_callable_bp"),
        ("callable_k100_t90_no_blacklist_bp", "structural_callable_bp"),
        ("effective_callable_k100_no_blacklist_bp", "structural_callable_bp"),
        ("structural_callable_bp", "nongap_bp"),
        ("nongap_bp", "raw_bp"),
    ]
    for value_col, bound_col in upper_bound_checks:
        if value_col not in territory.columns or bound_col not in territory.columns:
            continue
        bad = territory[value_col].astype(float) > territory[bound_col].astype(float) + float_tolerance
        for label, observed, bound in zip(labels[bad], territory.loc[bad, value_col], territory.loc[bad, bound_col]):
            warnings.warn(
                f"{value_col} exceeds {bound_col} for {label}: {observed} > {bound}.",
                stacklevel=2,
            )


def _check_region_additivity(
    chrom_territory: pd.DataFrame,
    region_territory: pd.DataFrame,
    float_tolerance: float = 1e-6,
) -> None:
    """Warn if X/Y per-region metrics do not sum back to whole chromosomes.

    The region table only covers X/Y in this workflow. Checking every quantitative
    column catches manual edits or PAR policies that silently break additivity.
    """
    if region_territory.empty:
        return

    value_columns = [
        column for column in region_territory.columns
        if column.endswith("_bp") or column.startswith("effective")
    ]
    chrom_index = chrom_territory.copy()
    chrom_index["chrom"] = chrom_index["chrom"].astype(str)
    chrom_index = chrom_index.set_index("chrom")
    region_sums = region_territory.copy()
    region_sums["chrom"] = region_sums["chrom"].astype(str)
    region_sums = region_sums.groupby("chrom", observed=True)[value_columns].sum()

    for chrom, row in region_sums.iterrows():
        if chrom not in chrom_index.index:
            warnings.warn(f"Region table contains chromosome {chrom}, absent from chromosome table.", stacklevel=2)
            continue
        for column in value_columns:
            expected = float(chrom_index.loc[chrom, column])
            actual = float(row[column])
            if abs(actual - expected) > float_tolerance:
                warnings.warn(
                    f"Region {column} for {chrom} sums to {actual:g}, expected {expected:g}. "
                    "Check xy_regions intervals, manual PAR edits, or par_policy.",
                    stacklevel=2,
                )


def chromosome_rows_as_regions(chrom_territory: pd.DataFrame) -> pd.DataFrame:
    """Represent whole-chromosome rows as region-aware rows.

    Autosomes use the chromosome name itself as both ``region`` and
    ``region_group``. This is the schema used by the unified output.
    """
    out = chrom_territory.copy()
    out["chrom"] = out["chrom"].astype(str)
    out.insert(1, "region", out["chrom"].astype(str))
    out.insert(2, "region_group", out["chrom"].astype(str))
    return out[UNIFIED_ID_COLUMNS + [column for column in TERRITORY_BP_COLUMNS if column in out.columns]]



def copy_x_par_mappability_to_y(region_territory: pd.DataFrame) -> pd.DataFrame:
    """Copy X-PAR mappability-derived columns to the homologous Y-PAR rows.

    This is intended for Umap runs on a FASTA where Y-PAR was masked. In that
    setup, PAR1Y/PAR2Y have zero bigWig signal by construction, even though
    their biological sequence is the same as PAR1X/PAR2X. Only mappability and
    callable-mappability columns are copied; raw/nongap/structural territory is
    left as computed from the Y intervals.
    """
    out = region_territory.copy()
    if "region" not in out.columns:
        return out

    par_pairs = {"PAR1X": "PAR1Y", "PAR2X": "PAR2Y"}
    copy_columns = [column for column in MAPPABILITY_BP_COLUMNS if column in out.columns]
    for source_region, target_region in par_pairs.items():
        source_mask = out["region"].astype(str).eq(source_region)
        target_mask = out["region"].astype(str).eq(target_region)
        if not source_mask.any() or not target_mask.any():
            warnings.warn(
                f"Could not copy {source_region} to {target_region}: one of the rows is missing.",
                stacklevel=2,
            )
            continue
        if source_mask.sum() != 1 or target_mask.sum() != 1:
            warnings.warn(
                f"Expected one row each for {source_region}/{target_region}; found "
                f"{source_mask.sum()} and {target_mask.sum()}. Using the first source row.",
                stacklevel=2,
            )
        source_values = out.loc[source_mask, copy_columns].iloc[0]
        out.loc[target_mask, copy_columns] = source_values.to_numpy()
    return cast_numeric_territory_columns_to_int64(out)


def whole_xy_rows_from_regions(region_territory: pd.DataFrame) -> pd.DataFrame:
    """Build whole-X and whole-Y rows by summing their subregions.

    This keeps complete X/Y rows consistent with any region-level post-processing
    such as copying PAR1X/PAR2X mappability values to PAR1Y/PAR2Y.
    """
    value_columns = [column for column in TERRITORY_BP_COLUMNS if column in region_territory.columns]
    rows = []
    for chrom in ["X", "Y"]:
        chrom_regions = region_territory[region_territory["chrom"].astype(str).eq(chrom)]
        if chrom_regions.empty:
            continue
        row = {"chrom": chrom, "region": chrom, "region_group": chrom}
        for column in value_columns:
            row[column] = int(pd.to_numeric(chrom_regions[column], errors="coerce").fillna(0).sum())
        rows.append(row)
    return pd.DataFrame(rows, columns=UNIFIED_ID_COLUMNS + value_columns)




def _prefix_output_path(out_prefix: str | Path, suffix: str) -> Path:
    """Return ``<out_prefix>_<suffix>.txt`` and ensure its parent exists."""
    prefix_path = Path(out_prefix)
    output_path = prefix_path.with_name(prefix_path.name + f"_{suffix}.txt")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    return output_path


def order_label_series(labels: pd.Series, custom_order: list[str]) -> pd.Categorical:
    """Return labels as an ordered categorical with unknown labels at the end."""
    seen = set()
    categories = []
    for label in custom_order:
        label = str(label)
        if label not in seen:
            seen.add(label)
            categories.append(label)
    for label in labels.astype(str):
        if label not in seen:
            seen.add(label)
            categories.append(label)
    return pd.Categorical(labels.astype(str), categories=categories, ordered=True)


def table_from_chrom_rows(chrom_territory: pd.DataFrame) -> pd.DataFrame:
    """Return the chromosome-level output table: 1-22, X and Y."""
    value_columns = [column for column in TERRITORY_BP_COLUMNS if column in chrom_territory.columns]
    out = chrom_territory[["chrom"] + value_columns].copy()
    out["chrom"] = out["chrom"].astype(str)
    out = cast_numeric_territory_columns_to_int64(out)
    out["chrom"] = order_label_series(out["chrom"], PRIMARY_CHROMS)
    out = out.sort_values("chrom").reset_index(drop=True)
    out["chrom"] = out["chrom"].astype(str)
    return out


def autosome_rows_for_single_label(
    chrom_territory: pd.DataFrame,
    label_column: str,
) -> pd.DataFrame:
    """Return autosome chromosome rows renamed to a single label column."""
    value_columns = [column for column in TERRITORY_BP_COLUMNS if column in chrom_territory.columns]
    autosomes = chrom_territory[~chrom_territory["chrom"].astype(str).isin(["X", "Y"])].copy()
    autosomes["chrom"] = autosomes["chrom"].astype(str)
    out = autosomes.rename(columns={"chrom": label_column})[[label_column] + value_columns].copy()
    return cast_numeric_territory_columns_to_int64(out)


def aggregate_region_group_table(
    chrom_territory: pd.DataFrame,
    region_territory: pd.DataFrame,
) -> pd.DataFrame:
    """Return the region-level output table.

    The user-facing ``region`` labels are autosomes plus the coarse X/Y labels
    from ``xy_regions.region_group``:

    1-22 + PAR1X + NPX + PAR2X + PAR1Y + NPY + PAR2Y
    """
    value_columns = [column for column in TERRITORY_BP_COLUMNS if column in region_territory.columns]
    autosomes = autosome_rows_for_single_label(chrom_territory, "region")

    xy = region_territory.copy()
    xy["region"] = xy["region_group"].astype(str)
    xy_grouped = (
        xy.groupby("region", observed=True, sort=False)[value_columns]
        .sum()
        .reset_index()
    )
    out = pd.concat([autosomes, xy_grouped], ignore_index=True)
    out = cast_numeric_territory_columns_to_int64(out)
    region_order = PRIMARY_CHROMS[:22] + ["PAR1X", "NPX", "PAR2X", "PAR1Y", "NPY", "PAR2Y"]
    out["region"] = order_label_series(out["region"], region_order)
    out = out.sort_values("region").reset_index(drop=True)
    out["region"] = out["region"].astype(str)
    return out[["region"] + value_columns]


def aggregate_region_table(
    chrom_territory: pd.DataFrame,
    region_territory: pd.DataFrame,
) -> pd.DataFrame:
    """Return the fine-region output table.

    The user-facing ``region_group`` labels are autosomes plus the fine X/Y
    labels from ``xy_regions.region``:

    1-22 + PAR1X + XAR + XCR + XTR + PAR2X + PAR1Y + NPY + PAR2Y

    The column is named ``region_group`` to match the requested output file,
    although the X labels are the fine ``region`` labels from ``xy_regions``.
    """
    value_columns = [column for column in TERRITORY_BP_COLUMNS if column in region_territory.columns]
    autosomes = autosome_rows_for_single_label(chrom_territory, "region_group")

    xy = region_territory.copy()
    xy["region_group"] = xy["region"].astype(str)
    xy_fine = xy[["region_group"] + value_columns].copy()

    out = pd.concat([autosomes, xy_fine], ignore_index=True)
    out = cast_numeric_territory_columns_to_int64(out)
    region_group_order = (
        PRIMARY_CHROMS[:22]
        + ["PAR1X", "XAR", "XCR", "XTR", "PAR2X", "PAR1Y", "NPY", "PAR2Y"]
    )
    out["region_group"] = order_label_series(out["region_group"], region_group_order)
    out = out.sort_values("region_group").reset_index(drop=True)
    out["region_group"] = out["region_group"].astype(str)
    return out[["region_group"] + value_columns]


def build_three_level_territories(
    chrom_info_path: str | Path,
    gap_path: str | Path,
    cytoband_path: str | Path,
    umap_k100_bw_path: str | Path,
    blacklist_path: str | Path,
    xy_regions_path: str | Path,
    structural_exclude_stains: tuple[str, ...] = DEFAULT_STRUCTURAL_EXCLUDE_STAINS,
    par_policy: str = "reference_hg38_duplicated",
    copy_x_par_to_y: bool = True,
    validate: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build chromosome, coarse-region and fine-region territory tables."""
    chrom_territory_raw = build_chrom_territory(
        chrom_info_path=chrom_info_path,
        gap_path=gap_path,
        cytoband_path=cytoband_path,
        umap_k100_bw_path=umap_k100_bw_path,
        blacklist_path=blacklist_path,
        structural_exclude_stains=structural_exclude_stains,
        par_policy=par_policy,
        validate=validate,
    )
    chrom_territory_raw["chrom"] = chrom_territory_raw["chrom"].astype(str)

    region_territory = build_region_territory(
        chrom_info_path=chrom_info_path,
        gap_path=gap_path,
        cytoband_path=cytoband_path,
        umap_k100_bw_path=umap_k100_bw_path,
        blacklist_path=blacklist_path,
        xy_regions_path=xy_regions_path,
        structural_exclude_stains=structural_exclude_stains,
        par_policy=par_policy,
        validate=validate,
    )
    region_territory["chrom"] = region_territory["chrom"].astype(str)
    region_territory = cast_numeric_territory_columns_to_int64(region_territory)

    if copy_x_par_to_y:
        region_territory = copy_x_par_mappability_to_y(region_territory)
    elif validate:
        _check_region_additivity(chrom_territory_raw, region_territory)

    # Chromosome table: autosomes are the directly computed chromosome rows;
    # X/Y are rebuilt from the corrected subregions so the chromosome table is
    # internally consistent with PAR1Y/PAR2Y copy policy.
    autosome_chrom_rows = chrom_territory_raw[
        ~chrom_territory_raw["chrom"].astype(str).isin(["X", "Y"])
    ].copy()
    xy_whole = whole_xy_rows_from_regions(region_territory).drop(columns=["region", "region_group"])
    chrom_territory = pd.concat([autosome_chrom_rows, xy_whole], ignore_index=True)
    chrom_territory = cast_numeric_territory_columns_to_int64(chrom_territory)
    chrom_table = table_from_chrom_rows(chrom_territory)

    region_group_table = aggregate_region_group_table(chrom_territory, region_territory)
    region_table = aggregate_region_table(chrom_territory, region_territory)

    for table in (chrom_table, region_table, region_group_table):
        table.attrs["par_policy"] = par_policy
        table.attrs["structural_exclude_stains"] = tuple(structural_exclude_stains)
        table.attrs["copy_x_par_to_y"] = bool(copy_x_par_to_y)

    if validate:
        validate_territory_table(chrom_table, label_columns=("chrom",))
        validate_territory_table(region_table, label_columns=("region",))
        validate_territory_table(region_group_table, label_columns=("region_group",))

    return chrom_table, region_table, region_group_table


def make_long_territory(
    territory: pd.DataFrame,
    id_columns: list[str] | tuple[str, ...],
    include_mb: bool = False,
) -> pd.DataFrame:
    """Convert a wide territory table to long format."""
    value_columns = [column for column in TERRITORY_BP_COLUMNS if column in territory.columns]
    long_form = territory.melt(
        id_vars=list(id_columns),
        value_vars=value_columns,
        var_name="territory_metric",
        value_name="bp",
    )
    long_form["bp"] = pd.to_numeric(long_form["bp"], errors="coerce").fillna(0).astype("int64")
    if include_mb:
        long_form["mb"] = long_form["bp"] / 1_000_000
    return long_form


def write_metadata(
    path: str | Path,
    *,
    par_policy: str,
    structural_exclude_stains: tuple[str, ...],
    chrom_info_path: str | Path,
    gap_path: str | Path,
    cytoband_path: str | Path,
    umap_k100_bw_path: str | Path,
    blacklist_path: str | Path,
    xy_regions_path: str | Path,
    copy_x_par_to_y: bool = True,
    out_prefix: str | Path | None = None,
) -> None:
    """Write a small JSON sidecar documenting key interpretation choices."""
    metadata = {
        "par_policy": par_policy,
        "structural_exclude_stains": list(structural_exclude_stains),
        "copy_x_par_to_y": bool(copy_x_par_to_y),
        "out_prefix": str(out_prefix) if out_prefix is not None else None,
        "outputs": {
            "chrom": str(_prefix_output_path(out_prefix, "chrom")) if out_prefix is not None else None,
            "region": str(_prefix_output_path(out_prefix, "region")) if out_prefix is not None else None,
            "region_group": str(_prefix_output_path(out_prefix, "region_group")) if out_prefix is not None else None,
        },
        "inputs": {
            "chrom_info": str(chrom_info_path),
            "gap": str(gap_path),
            "cytoband": str(cytoband_path),
            "umap_k100_bw": str(umap_k100_bw_path),
            "blacklist": str(blacklist_path),
            "xy_regions": str(xy_regions_path),
        },
        "notes": [
            "The supplied bigWig defines mappability values; par_policy records how PAR should be interpreted.",
            "When copy_x_par_to_y is true, PAR1Y/PAR2Y mappability-derived values are copied from PAR1X/PAR2X before aggregation.",
            "Chromosome-level X/Y rows are rebuilt from the corrected region rows so Y is consistent with copied PAR values.",
            "output_region.txt groups X fine regions by xy_regions.region_group, so XAR/XCR/XTR are summed as NPX.",
            "output_region_group.txt uses the fine xy_regions.region labels, despite the requested output column being named region_group.",
            "All territory metric columns are rounded and written as int64.",
        ],
    }
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")


def main() -> None:
    """Parse command-line arguments and write three territory aggregation tables."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chrom-info", required=True)
    parser.add_argument("--gap", required=True)
    parser.add_argument("--cytoband", required=True)
    parser.add_argument("--umap-k100-bw", required=True)
    parser.add_argument("--blacklist", required=True)
    parser.add_argument(
        "--xy-regions",
        required=True,
        help="xy_regions.txt partition used to build X/Y region and region_group outputs.",
    )
    parser.add_argument(
        "--out-prefix",
        required=True,
        help=(
            "Output prefix. For example --out-prefix dir/to/out/output writes "
            "dir/to/out/output_chrom.txt, output_region.txt and output_region_group.txt."
        ),
    )
    parser.add_argument(
        "--write-long",
        action="store_true",
        help=(
            "Also write <out_prefix>_chrom_long.txt, <out_prefix>_region_long.txt and "
            "<out_prefix>_region_group_long.txt."
        ),
    )
    parser.add_argument(
        "--long-include-mb",
        action="store_true",
        help="Add an mb column to long outputs. Disabled by default so bp remains int64-only.",
    )
    parser.add_argument(
        "--no-copy-x-par-to-y",
        action="store_true",
        help=(
            "Do not copy PAR1X/PAR2X mappability-derived values to PAR1Y/PAR2Y. "
            "By default the copy is applied before chromosome/region aggregation."
        ),
    )
    parser.add_argument(
        "--par-policy",
        default="reference_hg38_duplicated",
        choices=PAR_POLICY_CHOICES,
        help=(
            "Metadata label for pseudoautosomal-region interpretation. "
            "The calculation is still driven by the supplied bigWig."
        ),
    )
    parser.add_argument(
        "--structural-exclude-stains",
        default=",".join(DEFAULT_STRUCTURAL_EXCLUDE_STAINS),
        help=(
            "Comma-separated cytoBand gieStain values to exclude together with gaps "
            "before structural_callable_bp is computed. Default: acen,gvar,stalk. "
            "Use acen to keep gvar/stalk regions such as PAR2Y in the structural territory."
        ),
    )
    parser.add_argument(
        "--out-metadata",
        default=None,
        help=(
            "Optional JSON sidecar. If omitted, metadata is written to "
            "<out_prefix>_metadata.json. Pass an empty string to disable."
        ),
    )
    parser.add_argument(
        "--no-validate",
        action="store_true",
        help="Skip territory consistency warnings.",
    )
    args = parser.parse_args()
    structural_exclude_stains = _parse_comma_list(args.structural_exclude_stains)

    chrom_table, region_table, region_group_table = build_three_level_territories(
        chrom_info_path=_require_existing_file(args.chrom_info, "--chrom-info"),
        gap_path=_require_existing_file(args.gap, "--gap"),
        cytoband_path=_require_existing_file(args.cytoband, "--cytoband"),
        umap_k100_bw_path=_require_existing_file(args.umap_k100_bw, "--umap-k100-bw"),
        blacklist_path=_require_existing_file(args.blacklist, "--blacklist"),
        xy_regions_path=_require_existing_file(args.xy_regions, "--xy-regions"),
        structural_exclude_stains=structural_exclude_stains,
        par_policy=args.par_policy,
        copy_x_par_to_y=not args.no_copy_x_par_to_y,
        validate=not args.no_validate,
    )

    chrom_path = _prefix_output_path(args.out_prefix, "chrom")
    region_path = _prefix_output_path(args.out_prefix, "region")
    region_group_path = _prefix_output_path(args.out_prefix, "region_group")
    chrom_table.to_csv(chrom_path, sep="\t", index=False)
    region_table.to_csv(region_path, sep="\t", index=False)
    region_group_table.to_csv(region_group_path, sep="\t", index=False)

    if args.write_long:
        make_long_territory(chrom_table, ["chrom"], include_mb=args.long_include_mb).to_csv(
            _prefix_output_path(args.out_prefix, "chrom_long"), sep="\t", index=False
        )
        make_long_territory(region_table, ["region"], include_mb=args.long_include_mb).to_csv(
            _prefix_output_path(args.out_prefix, "region_long"), sep="\t", index=False
        )
        make_long_territory(region_group_table, ["region_group"], include_mb=args.long_include_mb).to_csv(
            _prefix_output_path(args.out_prefix, "region_group_long"), sep="\t", index=False
        )

    metadata_path = args.out_metadata
    if metadata_path is None:
        metadata_path = str(_prefix_output_path(args.out_prefix, "metadata")).replace(".txt", ".json")
    if metadata_path != "":
        write_metadata(
            metadata_path,
            par_policy=args.par_policy,
            structural_exclude_stains=structural_exclude_stains,
            chrom_info_path=args.chrom_info,
            gap_path=args.gap,
            cytoband_path=args.cytoband,
            umap_k100_bw_path=args.umap_k100_bw,
            blacklist_path=args.blacklist,
            xy_regions_path=args.xy_regions,
            copy_x_par_to_y=not args.no_copy_x_par_to_y,
            out_prefix=args.out_prefix,
        )

    print(f"Wrote {chrom_path}")
    print(f"Wrote {region_path}")
    print(f"Wrote {region_group_path}")
    if metadata_path != "":
        print(f"Wrote {metadata_path}")


if __name__ == "__main__":
    main()
