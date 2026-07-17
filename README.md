# Sex-aware reference genomes & mappability

Sex-chromosome-complement (SCC) informed reference genomes and single/multi-read
mappability tracks for **T2T-CHM13v2, GRCh38, GRCh37, GRCm39, GRCm38** (and
custom / T2T-mouse via BYO FASTA). Code lives here on GitHub; large binary
products (masked FASTAs, bigWigs) are released on **Zenodo** with checksums.

The approach follows the published consensus:
- Olney et al. 2020, *Biology of Sex Differences* — [10.1186/s13293-020-00312-9](https://doi.org/10.1186/s13293-020-00312-9)
- *Best practices for improving alignment and variant calling on human sex chromosomes*, AJHG 2025 — [PMC12190741](https://pmc.ncbi.nlm.nih.gov/articles/PMC12190741/)

## What "sex-aware" means here
| sample karyotype | reference to use | what is masked |
|---|---|---|
| **XX** (no Y) | `*.XX.fa` | the **entire chrY** |
| **XY** (has Y) | `*.XY.fa` | only the **chrY PAR** (chrX stays intact) |

Masking the **Y** PAR (never the X) sends pseudoautosomal reads unambiguously to
the single intact X copy. Pick the reference by the sample's actual karyotype —
infer it first with `scripts/infer_sex_complement.py` if unsure.

## Contents
```
config/assemblies.tsv          assembly registry (URLs, contig names, PAR status)
config/par/*.bed               verified PAR coordinates (human) / gated mouse
config/extra_mask/             optional stricter masking BEDs (XTR, ampliconic Y…)
config/ploidy/                 generated ploidy/PAR region files
scripts/mask_genome.py         ★ self-contained SCC masking (all organisms)
scripts/download_references.sh fetch source FASTAs + human territory tables
scripts/run_mappability.sh     ★ GenMap single/multi-read for k=24,36,50,100,150
scripts/genmap_to_umap_tracks.py  GenMap → Umap-style single/multi-read
scripts/run_umap_legacy.sh     optional Umap (py2.7/bowtie1) for UCSC-track parity
scripts/build_callable_beds.py callable/mappable BED from a mappability bigWig
scripts/make_ploidy_configs.py PAR + haploid region files for GATK/DeepVariant
scripts/infer_sex_complement.py QC gatekeeper (chrX/chrY coverage screen)
scripts/build_territory_tables.py  usable-sequence accounting (your 2A)
docs/PAR_coordinates.md        every PAR interval + its source
docs/aligner_indexes.md        "How do I build aligner indexes with these genomes?"
docs/methods.md                full methods
envs/                          conda envs + Dockerfile (GenMap default; Umap legacy)
```
★ = the two scripts most users need.

## Quick start (GRCh38, both flavours)
```bash
conda env create -f envs/environment.genmap.yml && conda activate sexaware-genmap

# 1. get the source genome (+ human territory tables/blacklist)
scripts/download_references.sh -a GRCh38 -o data/GRCh38

# 2. build BOTH SCC references (XX: whole-Y masked, XY: Y-PAR masked)
python scripts/mask_genome.py --assembly GRCh38 \
  --fasta data/GRCh38/hg38.fa.gz --complement both --outdir refs/GRCh38 --gzip

# 3. mappability (single + multi-read) for every k, per flavour
for c in XX XY; do
  scripts/run_mappability.sh -f refs/GRCh38/GRCh38.$c.fa.gz \
    -s refs/GRCh38/GRCh38.$c.chrom.sizes -o map/GRCh38 -l GRCh38.$c -t 8
done

# 4. callable BEDs + ploidy files (extra products)
python scripts/build_callable_beds.py \
  --bigwig map/GRCh38/GRCh38.XY/k100/GRCh38.XY.k100.single_read.bw \
  --threshold 1.0 --out callable/GRCh38.XY.k100.callable.bed
python scripts/make_ploidy_configs.py --par-bed config/par/GRCh38.PAR.bed \
  --chrom-sizes refs/GRCh38/GRCh38.XY.chrom.sizes --label GRCh38 --outdir config/ploidy
```

## Custom genome (e.g. T2T mouse or a strain assembly)
```bash
python scripts/mask_genome.py --fasta my_genome.fa.gz \
  --complement both --par-bed my_curated_par.bed \
  --x-contig chrX --y-contig chrY --outdir refs/custom
```

## ⚠️ Mouse PAR — read before using XY masking
The mouse PAR is a single ~700–960 kb block at the distal end of chrX/chrY,
**mostly unassembled**, and its boundary (in *Mid1* intron 3) **varies between
subspecies** — there is no clean canonical coordinate. Therefore:

- **XX masking (whole-Y) is reliable** for mouse; use it freely.
- **XY masking (Y-PAR) is gated**: `mask_genome.py` refuses mouse XY unless you
  pass `--allow-approximate-mouse-par` or supply a curated `--par-bed`.
- **GRCm38/mm10** intervals are literature-sourced (Morgan & Pardo-Manuel de
  Villena, *Genetics* 2019); **GRCm39/mm39** intervals are a **derived estimate**
  — liftOver mm10→mm39 or curate against the current *Mid1* boundary first.

Full provenance and the exact numbers: [`docs/PAR_coordinates.md`](docs/PAR_coordinates.md).

## Mappability engine
GenMap is the default (fast, containerisable, exact-match equivalent to Umap for
uniqueness). Legacy Umap is included only to reproduce pre-existing UCSC/ENCODE
bigWigs exactly. See [`docs/methods.md`](docs/methods.md).

## Building aligner indexes
Indexes are **not shipped** (large, version-specific). Recipes for BWA-MEM2,
bowtie2, STAR, salmon, kallisto and minimap2:
[`docs/aligner_indexes.md`](docs/aligner_indexes.md).

## Releasing to Zenodo
Track only code + light BEDs in git (`.gitignore` excludes big binaries). For a
release: build all references + tracks, then upload the FASTAs, bigWigs and their
`*.md5` to a Zenodo record and link its DOI here. Keep one Zenodo record per
organism to stay under size limits.

## Citation
See [`CITATION.cff`](CITATION.cff). Please also cite Olney 2020, the AJHG 2025
best-practices paper, and Umap/GenMap as appropriate.
