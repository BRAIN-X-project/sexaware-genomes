# How do I build aligner indexes with these genomes?

We deliberately **do not ship prebuilt aligner indexes**, since they are large and
tied to specific aligner versions. Build them yourself from the masked FASTA that
matches your sample's sex-chromosome complement (`*.XX.fa` or `*.XY.fa`). Always
index BOTH flavours and route each sample to the correct one (see
`infer_sex_complement.py`).

**Important**: pick the reference by the sample's karyotype. This means that for XY samples the `*.XY.fa` is necessary and for XX samples the `*.XX.fa` is necessary. For samples containing both
karyotypes you need to separate them and map each one with its own index. We will use the `*.XX.fa` sample as an example in each code.

**The annotation has to match the FASTA.** An `*.XX.*.gtf.gz` against an `*.XY.fa`
(or the reverse) is the single most common way to get silently wrong numbers: the
annotation names features on sequence the FASTA has replaced with N, so those
genes report 0 in every sample regardless of biology, and chrY genes vanish from
some matrices and not others. The pairs that are consistent by construction are:

| FASTA | annotation to pair with it | note |
|---|---|---|
| `GRCh38.XX.fa` | `GRCh38.XX.ensembl-116.gtf.gz` | chrY absent from both |
| `GRCh38.XY.fa` | `GRCh38.XY.noYPAR.ensembl-116.gtf.gz` | nothing annotated on masked sequence |
| `GRCh38.XY.fa` | `GRCh38.XY.ensembl-116.gtf.gz` | stock Ensembl; gene IDs stay comparable with the literature, but the ~64 Y-PAR genes are on N and will read 0 |

Use the stock `.XY.` unless you specifically want the guarantee that no gene in
the annotation sits on all-N sequence. Either way, **the XX and XY gene universes
differ** — read the warning in the README before merging count matrices across
complements.

**Files ship gzipped.** Several of the tools below cannot read `.gz`; decompress
once and reuse:

```bash
zcat GRCh38.XX.fa.gz > GRCh38.XX.fa
zcat GRCh38.XX.ensembl-116.gtf.gz > GRCh38.XX.ensembl-116.gtf
```

**Warning**: this section has not been troubleshooted for all combinations. If you find something failing that requires a correction, please open an issue and I'll address it.

## BWA-MEM2 (DNA)
```bash
bwa-mem2 index GRCh38.XX.fa
bwa-mem2 mem -t 16 GRCh38.XX.fa reads_R1.fq.gz reads_R2.fq.gz | samtools sort -o sample.bam
```

## Bowtie2 (DNA / short reads)
```bash
bowtie2-build --threads 16 GRCh38.XX.fa GRCh38.XX
bowtie2 -x GRCh38.XX -1 R1.fq.gz -2 R2.fq.gz -p 16 | samtools sort -o sample.bam
```

## STAR (RNA-seq)
Rebuild the genome index per flavour. STAR reads neither a gzipped FASTA nor a
gzipped GTF, so decompress both first.
```bash
zcat GRCh38.XX.fa.gz                   > GRCh38.XX.fa
zcat GRCh38.XX.ensembl-116.gtf.gz      > GRCh38.XX.ensembl-116.gtf

STAR --runMode genomeGenerate --genomeDir star_GRCh38_XX \
     --genomeFastaFiles GRCh38.XX.fa \
     --sjdbGTFfile GRCh38.XX.ensembl-116.gtf \
     --sjdbOverhang 100 --runThreadN 16
```
Repeat with `GRCh38.XY.fa` + `GRCh38.XY.ensembl-116.gtf` into `star_GRCh38_XY`.
The Ensembl GTFs here use bare seqnames (`1`, `X`, `Y`), which is what the
Ensembl FASTAs this repo masks use; a `chr`-prefixed GENCODE GTF will build an
index with zero junctions instead of failing loudly.

## Transcriptome FASTA (for salmon / kallisto)
Neither of the next two tools indexes a genome, so you need a transcriptome
first. Extracting it **from the sex-aware FASTA with the matching sex-aware GTF**
is self-consistent by construction: an all-N transcript cannot come out of a
`.noYPAR.` annotation, because nothing in it overlaps masked sequence.

```bash
# gffread is NOT in envs/environment.genmap.yml -- install it separately
# (conda install -c bioconda gffread).
gffread GRCh38.XX.ensembl-116.gtf -g GRCh38.XX.fa -w GRCh38.XX.transcripts.fa
```

Pair it with the `*.t2g.tsv` that `make_sexaware_gtf.py --emit-t2g` writes next
to each GTF: that is the transcript→gene map `tximport` needs, already restricted
to the same complement.

## salmon (decoy-aware, RNA quant)
Use the masked genome as the decoy so PAR/Y reads are absorbed correctly.
```bash
grep '^>' GRCh38.XX.fa | sed 's/^>//; s/ .*//' > decoys.txt
cat GRCh38.XX.transcripts.fa GRCh38.XX.fa > gentrome.fa
salmon index -t gentrome.fa -d decoys.txt -i salmon_GRCh38_XX -k 31 -p 16
salmon quant -i salmon_GRCh38_XX -l A -1 R1.fq.gz -2 R2.fq.gz -p 16 \
             -g GRCh38.XX.ensembl-116.t2g.tsv -o quant_sample
```

## kallisto (RNA quant)
kallisto indexes the transcriptome, so build the sex-aware one above and index
that — dropping Y transcripts for XX and Y-PAR transcripts for XY is exactly what
the shipped GTFs already did.
```bash
kallisto index -i GRCh38.XX.transcripts.idx GRCh38.XX.transcripts.fa
kallisto quant -i GRCh38.XX.transcripts.idx -o quant_sample -t 16 R1.fq.gz R2.fq.gz
# then in R: tximport(files, type="kallisto", tx2gene=read.delim("GRCh38.XX.ensembl-116.t2g.tsv"))
```

## minimap2 (long reads)
```bash
minimap2 -d GRCh38.XX.mmi GRCh38.XX.fa
minimap2 -ax map-ont GRCh38.XX.mmi reads.fq.gz | samtools sort -o sample.bam
```
