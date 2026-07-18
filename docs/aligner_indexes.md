# How do I build aligner indexes with these genomes?

We deliberately **do not ship prebuilt aligner indexes**, since they are large and
tied to specific aligner versions. Build them yourself from the masked FASTA that
matches your sample's sex-chromosome complement (`*.XX.fa` or `*.XY.fa`). Always
index BOTH flavours and route each sample to the correct one (see
`infer_sex_complement.py`).

**Important**: pick the reference by the sample's karyotype. This means that for XY samples the `*.XY.fa` is necessary and for XX samples the `*.XX.fa` is necessary. For samples containing both
karyotypes you need to separate them and map each one with its own index. We will use the `*.XX.fa` sample as an example in each code.

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
Rebuild the genome index per flavour; pass your GTF.
```bash
STAR --runMode genomeGenerate --genomeDir star_GRCh38_XY \
     --genomeFastaFiles GRCh38.XY.fa --sjdbGTFfile gencode.gtf \
     --sjdbOverhang 100 --runThreadN 16
```

## salmon (decoy-aware, RNA quant)
Use the masked genome as the decoy so PAR/Y reads are absorbed correctly.
```bash
grep '^>' GRCh38.XX.fa | sed 's/^>//; s/ .*//' > decoys.txt
cat transcripts.fa GRCh38.XX.fa > gentrome.fa
salmon index -t gentrome.fa -d decoys.txt -i salmon_GRCh38_XX -k 31 -p 16
```

## kallisto (RNA quant)
kallisto indexes the transcriptome; build a **sex-aware transcriptome** first by
dropping Y transcripts for XX and PAR-duplicate Y transcripts for XY, then:
```bash
kallisto index -i GRCh38.XX.transcripts.idx GRCh38.XX.transcripts.fa
```

## minimap2 (long reads)
```bash
minimap2 -d GRCh38.XX.mmi GRCh38.XY.fa
minimap2 -ax map-ont GRCh38.XX.mmi reads.fq.gz | samtools sort -o sample.bam
```
