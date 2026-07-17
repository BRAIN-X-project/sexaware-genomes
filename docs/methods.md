# Methods

## Rationale
The X and Y share regions of identical (PAR) or near-identical (XTR, gametologs,
ampliconic) sequence. Aligning every sample to a single default reference makes
reads from these regions multi-map, deflating mapping quality and creating both
false-positive and false-negative signal that differs systematically between XX
and XY samples. Aligning each sample to a reference matched to its
sex-chromosome complement (SCC) removes this artefact.

## Masking scheme
| sample | action | effect |
|--------|--------|--------|
| XX (no Y) | hard-mask the **entire chrY** | Y-derived misassignment eliminated; chrX diploid |
| XY (has Y) | hard-mask the **chrY PAR** only | PAR reads map only to intact chrX; chrY unique retained |

This is the Olney et al. 2020 consensus. Stricter, optional masking of XTR /
ampliconic Y / gametologs (AJHG 2025) is available via `--extra-mask` and the
`config/extra_mask/` BEDs, but by default those regions are left intact in the
FASTA and used only as downstream blacklists.

## Mappability
For each masked reference we compute single-read and multi-read mappability at
k ∈ {24, 36, 50, 100, 150} with **GenMap** (default engine; exact matching,
`-E 0`). Single-read = uniqueness of the k-mer starting at a position;
multi-read = fraction of the k overlapping k-mers that are unique (Umap's
definition), reproduced by `genmap_to_umap_tracks.py`. Legacy **Umap** (bowtie1,
py2.7) is provided (`run_umap_legacy.sh`) for exact continuity with UCSC/ENCODE
tracks. Because masked Y / Y-PAR bases read as 0 mappability, every downstream
callable set is automatically SCC-consistent.

## Territory / callable accounting
`build_territory_tables.py` (from the original `2A`) turns gaps, structural
cytobands, the ENCODE blacklist and the k100 multi-read bigWig into per-chromosome
and per-X/Y-region usable-sequence tables, with PAR-aware handling so PAR1Y/PAR2Y
mappability can be filled from PAR1X/PAR2X on Y-PAR-masked runs.

## Coordinates
All PAR intervals and their sources are documented in `PAR_coordinates.md`. Human
builds are verified; **mouse PAR is caveat-heavy** (see the warning there and in
the README) and XY masking on mouse is gated behind an explicit flag.
