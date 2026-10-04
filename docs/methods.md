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
ampliconic Y / gametologs (Taravella Oill 2026, `10.1016/j.ajhg.2026.02.019`) is available via `--extra-mask` and by including 
additional BEDs with these regions in the
`config/extra_mask/` folder, but by default those regions are left intact in the
FASTA and used only as downstream blacklists.

## Mappability
For each masked reference we compute single-read and multi-read mappability at
k ∈ {24, 36, 50, 100, 150} with **GenMap** (default engine; exact matching,
`-E 0`). Single-read = uniqueness of the k-mer starting at a position;
multi-read = fraction of the k overlapping k-mers that are unique (Umap's
definition), reproduced by `genmap_to_umap_tracks.py`. 

## Territory / callable accounting
`build_territory_tables.py` turns gaps, structural
cytobands, the ENCODE blacklist and the k100 multi-read bigWig into per-chromosome
and per-X/Y-region usable-sequence tables, with PAR-aware handling so PAR1Y/PAR2Y
mappability can be filled from PAR1X/PAR2X on Y-PAR-masked runs.

## Annotation filtering
Published GTFs follow a single rule, from which everything else is derived:

> **A sex-aware GTF contains no feature that overlaps N-masked sequence in the
> corresponding sex-aware FASTA.**

| variant | mask applied to the FASTA | annotation consequence |
|---|---|---|
| `XX` | whole chrY | every chrY feature removed |
| `XY` | Y-side PAR | none — upstream annotation passed through unchanged |
| `XY-noYPAR` | Y-side PAR | features on the Y-side PAR removed |

Both `XY` flavours are published because the rule and common practice disagree
and the disagreement is defensible either way: `.XY.` is what every pipeline
expects and keeps gene IDs comparable with the rest of the literature, while
`.XY.noYPAR.` is the one that cannot report a count for a gene whose sequence is
all N. Pick `.XY.` unless you specifically want the stricter guarantee.

The filter's input is therefore the `refs/<ASM>/<ASM>.<C>.masked.bed` that
`mask_genome.py` actually wrote, **not** `config/par/*.bed`. The FASTA and the
GTF then cannot silently diverge, and any `--extra-mask` BED (XTR, ampliconic Y,
gametologs) is honoured by the annotation automatically. `--par-bed` remains as
a documented fallback for anyone without the `masked.bed`.

**Overlap rule.** `--overlap-rule any` (default) drops a feature that overlaps
the mask at all; `contained` drops only features entirely inside it. **Drop
unit** `gene` (default) removes the whole gene when any of its features is hit,
so no transcript model survives with exons silently missing. A gene that
straddles the mask boundary is dropped whole, recorded with
`reason=y_par_overlap_partial`, counted in `genes_partially_overlapping` and
logged as a warning; `--on-partial-overlap fail` turns that into exit 3. This is
not hypothetical: in GRCh38 r116 the lncRNA **XGY2** (`ENSG00000290840`,
Y:2,752,083-2,854,641) spans the PAR1 boundary at 2,781,479.

**The drop set is keyed on `(seqname, gene_id)`, never on `gene_id` alone.**
Whether a gene with an X and a Y PAR copy gets one shared ID or two distinct
ones is a per-release convention: Ensembl GRCh38 r116 gives the Y copies their
own IDs (PPP2R3B on Y is `ENSG00000292327`, different from the X's `ENSG00000167393`), but
nothing in the format requires it and the GRCh37 REST API still reports both
locations for a single ID. Under the sharing convention, keying on a bare ID
would delete the **X** copy as well — the exact inverse of what an XY reference
is for, and it would publish looking correct. Keying on the pair is right under
either convention. The invariant that catches it is simply that *chrX must have
the same number of features in and out*; it is asserted on every run and a
failure is exit 4.

**Restricting to the primary assembly is not sex-aware filtering.** Ensembl
release 75 ships only an all-scaffolds GTF, so 204 of its 265 seqnames (MHC
haplotypes, patches, LRG regions) are absent from the `primary_assembly` FASTA
this repo masks. Those rows are dropped under `restrict_to_genome=1` in
`config/annotations.tsv`, and the counts are reported and stored **separately**
(`genes_dropped_off_genome` vs `genes_dropped_sexaware`) — reporting one total
would tell a user their standard `XY` annotation had been sex-aware filtered
when it had not.

### Expected counts
A new release that deviates from these is the signal that the PAR annotation
changed, not that the tool broke. 

| assembly / release | XX genes dropped | XY-noYPAR genes dropped | off-genome |
|---|---|---|---|
| GRCh38 r116 | 672 | **64** (1 straddling: XGY2) | 0 |
| GRCh38 r110 | 601 | **49** (1 straddling: XGY2) | 0 |
| GRCh38 r98  | 522 | **0** | 0 |
| GRCh37 r87  | 495 | **0** | 0 |
| GRCh37 r75  | 495 | **0** | 5,772 |
| GRCm39 r116 | 2,052 | **6** | 0 |
| GRCm39 r110 | 1,633 | **4** | 0 |
| GRCm38 r102 | 1,569 | **5** | 0 |
| GRCm38 r98  | 1,569 | **5** | 0 |



**Mouse `XY-noYPAR` drops real, expressed sequence on a derived boundary — read
this before using it.** The 6 genes removed in GRCm39 r116 are `Mid1-ps1`,
`G530011O06Riky`, `Erdr1y`, `Gm58025`, `Gm21748` and `Gm21742`; GRCm38 r102
removes `Gm21860`, `Mid1-ps1`, `Gm47283`, `Gm21742` and a protein-coding
`Gm21748`. The set is release-dependent: GRCm39 r110 and GRCm38 r98 drop only
`Gm21860`, `Mid1-ps1`, `Gm47283`, `Gm21742` (+ `Gm21748` on GRCm38), because
`Erdr1y` and `Gm58025` were not yet placed there. **`Erdr1y` is an expressed lncRNA**, not a pseudogene artefact. Only
~98.9 kb (16.5%) of the 598.9 kb derived mouse PAR is real sequence in GRCm39's
`primary_assembly` — the rest is already N — so the masking removes roughly
99 kb of genuine, annotated sequence using an interval whose true boundary lies
in *Mid1* intron 3 and **varies between subspecies**. This is why both
`mask_genome.py` and `make_sexaware_gtf.py` refuse mouse `XY-noYPAR` unless you
pass `--allow-approximate-mouse-par` or supply a curated `--par-bed`, and why
`par_status: approximate` is recorded in the metadata and surfaced in
`MANIFEST.tsv`. For mouse, prefer `.XY.`; use `.XY.noYPAR.` only with a curated
boundary.

**Whether the Y PAR is annotated at all is a genebuild-era convention, not an
assembly property.** Three releases of the *same* assembly disagree:

| release | genes on Y | first gene start | last gene end | Y-PAR annotated? |
|---|---|---|---|---|
| GRCh38 r98 (2019)  | 522 | 2,784,749 | 56,855,488 | **no** |
| GRCh38 r110 (2023) | 601 |   253,743 | 57,214,397 | yes (49 genes) |
| GRCh38 r116 (2026) | 672 |   253,743 | 57,215,634 | yes (64 genes) |

GRCh38's PAR1 ends at 2,781,479 and PAR2 starts at 56,887,902, so r98's entire Y
genebuild sits strictly inside NPY even though that assembly's Y PAR is real
ACGT, not N. Ensembl only began placing PAR gene models on the Y copy somewhere
between r98 and r110. Both GRCh37 releases likewise annotate nothing there
(first gene on Y at 2,652,790, past the 2,649,520 PAR1 end; last ending at
59,001,635, before the 59,034,049 PAR2 start), which in GRCh37's case is
*over*determined — its `primary_assembly` has the PAR as all-N as well.

The practical consequence: `GRCh38.XY.noYPAR.ensembl-98`, `GRCh37.XY.noYPAR.*-87`
and `GRCh37.XY.noYPAR.*-75` come out content-identical to their `XY`
counterparts -- identical in every feature line, though not as files: each keeps
its own `#!sexaware-variant` / `#!sexaware-mask-bed` / `#!sexaware-note` header,
so their md5s differ by design. The tool detects the content equivalence,
records `identical_to_standard_xy` in the metadata and surfaces it in
`MANIFEST.tsv`; the files are still published so a script building a filename
programmatically does not hit a 404. **Do not read a
0 in the `XY-noYPAR` column as evidence that the filter ran and found nothing to
do in an assembly where it should have** — check the release.

## Annotation sources and versions
Releases are **pinned**, never floating, and are chosen for what the ecosystem
actually uses — each non-latest row is the Ensembl release behind a widely
deployed GENCODE/Cell Ranger vintage:

| assembly | releases | why |
|---|---|---|
| GRCh38 | 116, 110, 98 | current; r110 = GENCODE v44 / Cell Ranger 2024-A; r98 = GENCODE v32 / Cell Ranger 2020-A |
| GRCh37 | 87, 75 | the grch37 branch is frozen at 87; r75 is the classic GENCODE v19-era pin |
| GRCm39 | 116, 110 | current; r110 |
| GRCm38 | 102, 98 | r102 is the last Ensembl release on GRCm38; r98 = GENCODE vM23 / Cell Ranger mm10-2020-A |

Only the plain `<Species>.<ASM>.<REL>.gtf.gz` flavour is used, because it is the
one that pairs with `dna.primary_assembly`. `chr_patch_hapl_scaff` is excluded
(it adds patch and haplotype scaffolds the primary assembly does not contain)
and `abinitio` is excluded (predictions, not the genebuild);
`download_annotations.sh` refuses either URL before downloading, whatever the
registry says.

The pairing policy is **pin the release, then verify, never assume**: every run
checks the GTF's seqname set against the paired `chrom.sizes` and every feature
end against the contig length. An unknown seqname or an out-of-range coordinate
is exit 2 by default. That is what turns "this GTF goes with this FASTA" from an
assumption into a checked fact — GRCh38 r116, for instance, uses 70 of the
assembly's 194 contigs and exceeds none of their lengths.

T2T-CHM13v2 has **no annotation** here: Ensembl's main branch does not carry it.
Use RefSeq `GCF_009914755.1` if you need one.

## Reproducibility
Outputs are byte-reproducible, so a published md5 is a real identity check
rather than a record of one particular run:

- Source releases are pinned in `config/sources.tsv` and `config/annotations.tsv`.
- gzip output is written with `mtime=0` and no stored filename, since
  `gzip.open()` otherwise stamps both into the header.
- No build timestamp goes into a FASTA or GTF header; it goes to the
  `metadata.json` instead.
- Every `.md5` sidecar and every `MANIFEST.tsv` row holds the md5 of the
  **distributed bytes** (the `.gz` as shipped), which is what `md5sum -c`
  checks and what Zenodo reports per file in its API. `metadata.json` also
  carries `output_fasta_uncompressed_md5` for comparing content across
  providers, where the compressed bytes legitimately differ.

## Coordinates
All PAR intervals and their sources are documented in `PAR_coordinates.md`. Human
builds are verified; **mouse PAR is caveat-heavy** (see the warning there and in
the README) and XY masking on mouse is gated behind an explicit flag.

GTF is **1-based, inclusive**; BED is **0-based, half-open**. Both predicates are
single helpers in `make_sexaware_gtf.py`, each with a unit test:

```python
overlaps_bed      = (gtf_start - 1) <  bed_end   and gtf_end >  bed_start
contained_in_bed  = (gtf_start - 1) >= bed_start and gtf_end <= bed_end
```

### `config/xy_regions.txt` vs `config/par/*.bed` — audited, not a discrepancy
The two files give different numbers for the same PAR and it looks like a bug; it
is not. `xy_regions.txt` (consumed by `build_territory_tables.py`) is 1-based
inclusive and is a **complete partition**: on GRCh38 its X rows tile 1 →
156,040,895 and its Y rows tile 1 → 57,227,415, exactly the contig lengths, with
no gaps, because every base has to be accounted to exactly one region. The PAR
BED is 0-based half-open and is a **mask definition**: it excludes sequence there
is no point masking.

Converted to one convention, every *internal* boundary agrees to the base:

| boundary | `xy_regions.txt` (1-based) | `GRCh38.PAR.bed` (0-based → 1-based) |
|---|---|---|
| PAR1X end / XAR start | 2,781,479 | 2,781,479 |
| PAR2X start | 155,701,383 | 155,701,382 → 155,701,383 |
| PAR2Y start | 56,887,903 | 56,887,902 → 56,887,903 |

The files differ at exactly four places, all of them a chromosome terminus, all of
them 10,000 bp: X start, X end, Y start, Y end. Cross-checked against UCSC's
`gap.txt.gz`, those four blocks are the annotated telomeres —
`chrX 0-10000`, `chrX 156030895-156040895`, `chrY 0-10000`,
`chrY 57217415-57227415` — i.e. all-N by construction. The BED leaves them out
(masking N to N is a no-op); the partition has to include them to stay a
partition, and `build_territory_tables.py` subtracts them again from usable
sequence via the gap table. **The BED is the file used for masking and annotation
filtering**, and it is the one whose intervals land in `masked.bed` and in every
GTF header.
