# Sex-aware reference genomes & mappability

Sex-chromosome-complement (SCC) informed reference genomes and single/multi-read
mappability tracks for **T2T-CHM13v2, GRCh38, GRCh37, GRCm39, GRCm38** (and
custom / T2T-mouse via BYO FASTA). Code lives here on GitHub; large binary
products (masked FASTAs, bigWigs) are released on **Zenodo** with checksums.

The approach follows the published consensus:
- Olney et al. 2020, *Biology of Sex Differences* — [10.1186/s13293-020-00312-9](https://doi.org/10.1186/s13293-020-00312-9)
- *Best practices for improving alignment and variant calling on human sex chromosomes* — Taravella Oill et al., *Am J Hum Genet* **113**:782-793 (2026), [10.1016/j.ajhg.2026.02.019](https://doi.org/10.1016/j.ajhg.2026.02.019) 

## What "sex-aware" means here
| sample karyotype | reference to use | annotation to use | what is masked |
|---|---|---|---|
| **XX** (no Y) | `*.XX.fa` | `*.XX.<source>-<rel>.gtf.gz` | the **entire chrY** |
| **XY** (has Y) | `*.XY.fa` | `*.XY.<source>-<rel>.gtf.gz` (or `*.XY.noYPAR.*`) | only the **chrY PAR** (chrX stays intact) |

**The annotation's complement must match the FASTA's.** Feeding an unfiltered
GTF to an XX reference asks your counter for chrY genes whose sequence is all N:
it will return zeros and duplicated symbols rather than an error.

Masking the **Y** PAR (never the X) sends pseudoautosomal reads unambiguously to
the single intact X copy. Pick the reference by the sample's actual karyotype —
infer it first with `scripts/infer_sex_complement.py` if unsure.

## Contents
```
config/assemblies.tsv          assembly registry (contig names, PAR status)
config/sources.tsv             genome FASTA URLs per provider, with pinned releases
config/annotations.tsv         GTF URLs per (assembly, source, release, flavour)
config/par/*.bed               verified PAR coordinates (human) / gated mouse
config/extra_mask/             optional stricter masking BEDs (XTR, ampliconic Y…)
scripts/mask_genome.py         ★ self-contained SCC masking (all organisms)
scripts/make_sexaware_gtf.py   ★ sex-aware annotation, filtered against the FASTA's mask
scripts/download_references.sh fetch source FASTAs + human territory tables
scripts/download_annotations.sh fetch + checksum-verify a registered GTF
scripts/build_manifest.sh      MANIFEST.tsv (checksums + provenance per file)
scripts/run_mappability.sh     ★ GenMap single/multi-read for k=24,36,50,100,150
scripts/genmap_to_umap_tracks.py  GenMap → Umap-style single/multi-read
scripts/run_umap_legacy.sh     optional Umap (py2.7/bowtie1) for UCSC-track parity
                               (no env file shipped; the script header says how to build one)
scripts/build_callable_beds.py callable/mappable BED from a mappability bigWig
scripts/make_ploidy_configs.py PAR + haploid region files for GATK/DeepVariant
scripts/infer_sex_complement.py QC gatekeeper (chrX/chrY coverage screen)
scripts/build_territory_tables.py  usable-sequence accounting per chromosome / X-Y region
scripts/run_all_pipeline.sh    end-to-end rebuild of every published product
docs/PAR_coordinates.md        every PAR interval + its source
docs/aligner_indexes.md        "How do I build aligner indexes with these genomes?"
docs/methods.md                full methods
docs/release_checklist.md      what each Zenodo record contains + the v2 upload plan
tests/                         fixture-based regression suites (plain bash, 90 asserts, <3 s)
envs/                          GenMap conda env + Dockerfile
```
★ = the scripts most users need.

## Genome sources
Download provider is chosen by preference with graceful fallback, via
`--preferred-source` (default `ensembl,ucsc,ncbi`); URLs live in
`config/sources.tsv`. The first listed source that has the assembly wins.

| source | naming | notes |
|--------|--------|-------|
| ensembl | `1/X/Y` | `primary_assembly`; whether the chrY PAR already carries Ns **varies by assembly** (GRCh38: real sequence; GRCh37: all-N; mouse: ~84% N) — so never skip the XY step assuming it is a no-op |
| ucsc | `chr1/chrX/chrY` | soft-masked bigZips; T2T is `hs1` |
| ncbi | RefSeq accessions → renamed to `chrN` | primary source for T2T-CHM13v2 |

Contig naming is normalised by `mask_genome.py`, so the same PAR BEDs work
whatever source you pick. `download_references.sh` prints the exact FASTA path to
feed to `mask_genome.py`.

## Annotation (GTF)
Matching sex-aware annotation is published alongside each genome. 

| assembly | releases | notes |
|---|---|---|
| GRCh38 | **116**, 110, 98 | r110 = GENCODE v44 / Cell Ranger 2024-A; r98 = GENCODE v32 / Cell Ranger 2020-A |
| GRCh37 | **87**, 75 | the grch37 branch is frozen at 87; r75 is the classic GENCODE v19-era pin |
| GRCm39 | **116**, 110 | |
| GRCm38 | **102**, 98 | r102 is the last Ensembl release on GRCm38; r98 = GENCODE vM23 / Cell Ranger mm10-2020-A |
| T2T-CHM13v2 | — | Ensembl's main branch has no annotation for it; use RefSeq `GCF_009914755.1` |

Three variants per release:

| variant | contents | use it when |
|---|---|---|
| `.XX.` | chrY removed entirely | the sample has no Y |
| `.XY.` | standard Ensembl annotation, unchanged | **the usual choice** for XY samples |
| `.XY.noYPAR.` | Y-side PAR features removed too | you want the strict guarantee that no feature sits on N |

`.XY.` is what pipelines expect and keeps gene IDs comparable with the
literature; `.XY.noYPAR.` is the one that cannot report a count for a gene whose
sequence is entirely masked. Pick `.XY.` unless you specifically want the latter.
For **mouse**, prefer `.XY.`: the mouse PAR boundary is approximate, and
`.XY.noYPAR.` there drops ~99 kb of real annotated sequence
 (see [`docs/methods.md`](docs/methods.md)).


Their **checksums still differ**, and that is not a bug: the two files carry
different `#!sexaware-variant`, `#!sexaware-mask-bed` and `#!sexaware-note`
header lines, because they genuinely were produced under different filter
policies and the header is the record of that. Compare the feature bodies, not
the file digests, if you want to confirm the equivalence yourself. Ensembl started annotating the Y PAR on
GRCh38 only between r98 and r110, so this is a property of the release, not of
the assembly.

How many genes each variant removes, measured:

| release | `.XX.` drops | `.XY.noYPAR.` drops |
|---|---|---|
| GRCh38 r116 / r110 / r98 | 672 / 601 / 522 | 64 / 49 / **0** |
| GRCh37 r87 / r75 | 495 / 495 | **0** / **0** |
| GRCm39 r116 / r110 | 2,052 / 1,633 | 6 / 4 |
| GRCm38 r102 / r98 | 1,569 / 1,569 | 5 / 5 |


Build one yourself — the filter reads the mask `mask_genome.py` actually applied,
so the FASTA and the GTF cannot silently disagree:

```bash
# 1. the GTF (checksum-verified against Ensembl's CHECKSUMS)
scripts/download_annotations.sh -a GRCh38 -o data/GRCh38 --release 116

# 2. all three variants, filtered against refs/GRCh38/GRCh38.{XX,XY}.masked.bed
python scripts/make_sexaware_gtf.py --assembly GRCh38 \
  --gtf data/GRCh38/Homo_sapiens.GRCh38.116.gtf.gz \
  --variant all --outdir refs/GRCh38 --emit-t2g
```

Each run also writes `.metadata.json` (full provenance, including the md5 of the
paired FASTA), `.dropped_genes.tsv` (every gene removed, with the reason) and,
with `--emit-t2g`, a `.t2g.tsv` for salmon/kallisto + tximport.

### ⚠️ XX and XY have different gene universes
This is the thing that most often bites people using sex-aware annotation, and
it is silent. `featureCounts`/`HTSeq` on XX samples emit **no chrY rows at all**,
so your XX and XY count matrices have different row indexes — a naive `cbind` or
`merge` will drop genes or fill them with `NA` without complaining, and any
sex-stratified differential expression done on top of that is wrong.

Use the **XY gene set as the canonical index** and fill chrY with zeros for XX
samples. The `.dropped_genes.tsv` shipped with every XX annotation is exactly the
list you need to do that.

## Quick start (GRCh38, both flavours)
```bash
conda env create -f envs/environment.genmap.yml && conda activate sexaware-genmap

# 1. get the source genome (Ensembl by default; +human territory tables/blacklist)
scripts/download_references.sh -a GRCh38 -o data/GRCh38
#   ...or force a different provider (this one yields hg38.fa.gz, not the
#   Homo_sapiens.* name used in step 2 -- use the path step 1 prints):
#   scripts/download_references.sh -a GRCh38 -o data/GRCh38 --preferred-source ucsc

# 2. build BOTH SCC references (XX: whole-Y masked, XY: Y-PAR masked)
#    use the FASTA path that step 1 printed ("==> feed this FASTA ...")
python scripts/mask_genome.py --assembly GRCh38 \
  --fasta data/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
  --complement both --outdir refs/GRCh38 --gzip

# 2b. sex-aware annotation, filtered against the masks step 2 just wrote
scripts/download_annotations.sh -a GRCh38 -o data/GRCh38 --release 116
python scripts/make_sexaware_gtf.py --assembly GRCh38 \
  --gtf data/GRCh38/Homo_sapiens.GRCh38.116.gtf.gz \
  --variant all --outdir refs/GRCh38 --emit-t2g

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
  --chrom-sizes refs/GRCh38/GRCh38.XY.chrom.sizes --label GRCh38 --outdir ploidy

# 5. (optional) territory tables — feed OUR sex-aware k100 multi_read bigWig,
#    NOT UCSC's precomputed track (that one is built on the UNMASKED genome and
#    would report mappability inside the PAR).
python scripts/build_territory_tables.py \
  --chrom-info data/GRCh38/chromInfo.txt.gz --gap data/GRCh38/gap.txt.gz \
  --cytoband data/GRCh38/cytoBand.txt.gz --blacklist data/GRCh38/ENCFF356LFX.hg38.blacklist.bed.gz \
  --umap-k100-bw map/GRCh38/GRCh38.XY/k100/GRCh38.XY.k100.multi_read.bw \
  --xy-regions config/xy_regions.txt \
  --out-prefix territory/GRCh38.XY
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

## Downloads (Zenodo)

| record | DOI | contents |
|---|---|---|
| genomes + annotation | [10.5281/zenodo.21442231](https://doi.org/10.5281/zenodo.21442231) | masked FASTAs, `chrom.sizes`, masks, sex-aware GTFs |
| mappability | [10.5281/zenodo.21461059](https://doi.org/10.5281/zenodo.21461059) | single/multi-read bigWigs, k = 24, 36, 50, 100, 150 |

Assemblies: GRCh38, GRCh37, T2T-CHM13v2 (all ensembl except T2T = ucsc), GRCm39,
GRCm38. `MANIFEST.tsv` at the root of each record lists every file with its
size, md5, sha256 and provenance; it is generated, not kept in git — rebuild it
with `scripts/build_manifest.sh`. [`docs/release_checklist.md`](docs/release_checklist.md)
records what each record actually contains and how a new version is cut.
Verify a download with:

```bash
md5sum -c GRCh38.XY.fa.md5
```

> **⚠️ Erratum (v1 genome record).** The `*.fa.md5` in the first version
> recorded the md5 of the *decompressed* stream while naming the `.fa.gz`, so
> `md5sum -c` reported a mismatch on a perfectly good download. The FASTA content
> was never affected. Corrected from v2 on.

## Citation
See [`CITATION.cff`](CITATION.cff). Cite the **concept** DOIs —
[10.5281/zenodo.21442231](https://doi.org/10.5281/zenodo.21442231) for the
genomes and annotation, [10.5281/zenodo.21461059](https://doi.org/10.5281/zenodo.21461059)
for the mappability tracks — so the citation keeps resolving as new versions are
released. Please also cite Olney 2020, Taravella Oill 2026 (best practices),
Ensembl for the annotation, and Umap/GenMap as appropriate.
