# Pseudoautosomal region (PAR) coordinates — provenance

All intervals below are **0-based, half-open** (BED convention), matching the
files in `config/par/`. For **XY** references the pipeline masks the **chrY**
intervals; the chrX intervals are shipped for QC and territory accounting only.

## Human — verified

### GRCh38 / hg38
| region | contig | start | end |
|--------|--------|-------|-----|
| PAR1 | chrX | 10000 | 2781479 |
| PAR2 | chrX | 155701382 | 156030895 |
| PAR1 | chrY | 10000 | 2781479 |
| PAR2 | chrY | 56887902 | 57217415 |

X and Y PAR coordinates are identical in GRCh38. Source: GRCh38 assembly /
Ensembl & GATK reference documentation.

### GRCh37 / hg19
| region | contig | start | end |
|--------|--------|-------|-----|
| PAR1 | chrX | 60000 | 2699520 |
| PAR2 | chrX | 154931043 | 155260560 |
| PAR1 | chrY | 10000 | 2649520 |
| PAR2 | chrY | 59034049 | 59363566 |

Note: GRCh37 is the one build where PAR1 **differs** between X (starts 60,001)
and Y (starts 10,001), a ~50 kb offset. Source: GRCh37 assembly / UCSC & GATK
reference documentation.

### T2T-CHM13v2.0
| region | contig | start | end |
|--------|--------|-------|-----|
| PAR1 | chrX | 0 | 2394410 |
| PAR2 | chrX | 153925834 | 154259566 |
| PAR1 | chrY | 0 | 2458320 |
| PAR2 | chrY | 62122809 | 62460029 |

The Y is from HG002; PARs are larger and X/Y sizes differ slightly. Source:
CHM13v2.0 official PAR coordinate file (marbl/CHM13; T2T-Y paper,
biorxiv 2022.12.01.518724).

## Mouse — ⚠️ WARNING

The mouse PAR is a single ~700–960 kb block at the **distal (q-terminal)** end of
chrX/chrY, is **mostly unassembled repetitive sequence**, and the
pseudoautosomal boundary (PAB) is **biologically variable even between mouse
subspecies** (it sits in intron 3 of *Mid1*/*Trim18*). There is no single clean
"canonical" coordinate, so `mask_genome.py` **refuses XY masking on mouse unless
you pass `--allow-approximate-mouse-par` or supply your own curated `--par-bed`.**
XX masking (whole-Y mask) is unaffected and reliable.

### GRCm38 / mm10 — literature-sourced
| region | contig | start | end |
|--------|--------|-------|-----|
| PAR | chrX | 169969759 | 170931299 |
| PAR | chrY | 90745845 | 91644698 |

chrX interval = GRCm38.p5 PAR annotation (≈961,540 bp); the PAB reference
position is chrX:169,986,134. chrY interval from the mm10 assembled PAR extent.
Sources: Morgan & Pardo-Manuel de Villena, *Genetics* 2019 ("Instability of the
Pseudoautosomal Boundary in House Mice", 10.1534/genetics.119.302232;
PMC6553833); GRCm38.p5 annotation.

### GRCm39 / mm39 — DERIVED ESTIMATE from LiftOvering the starting chromosome
| region | contig | start | end |
|--------|--------|-------|-----|
| PAR | chrX | 168752755 | 169376592 |  
| PAR | chrY | 90757114 | 91355967 |

**Not from a primary source.** GRCm39 retiled the PAR (adding *Sts*, *Nlgn4l*,
*Akap17a*, *2510022D24Rik*) so mm10 numbers do not transfer directly. These
intervals are a conservative estimate obtained by preserving the mm10
distance-from-telomere. **Before any XY mouse run on GRCm39, replace this file by
`liftOver`-ing the GRCm38 PAR (mm10→mm39 chain) or by curating against the
current Ensembl/UCSC *Mid1* boundary.** Sources for the retiling: NCBI Insights
"RefSeq annotation of mouse GRCm39" (2020); GenomeRef GRCm39 notes.

### T2T mouse
No stable registered assembly/URL is shipped. Complete mouse assemblies now
exist (Nature Genetics 2025, s41588-025-02367-z, incl. a complete CAST/EiJ PAR),
but coordinates are strain-specific. Use custom mode: `--fasta <t2t_mouse.fa>
--par-bed <your_curated_par.bed> --x-contig chrX --y-contig chrY`.
