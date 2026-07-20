# 1. get the source genome (Ensembl by default; +human territory tables/blacklist)
for genome in GRCh38 GRCh37 T2T-CHM13v2 GRCm38 GRCm39; do
    if genome == T2T-CHM13v2; do
        source="ucsc"
    else; do
        source="ensembl"
    done

    scripts/download_references.sh -a $genome -o data/$genome --preferred-source $source

done


# 2. build BOTH SCC references (XX: whole-Y masked, XY: Y-PAR masked)
#    use the FASTA path that step 1 printed ("==> feed this FASTA ...")
for genome in GRCh38 GRCh37 T2T-CHM13v2 GRCm38 GRCm39; do
    python scripts/mask_genome.py --assembly $genome \
    --fasta data/$genome/Homo_sapiens.$genome.dna.primary_assembly.fa.gz \
    --complement both --outdir refs/$genome --gzip
done

for genome in T2T-CHM13v2; do
    python scripts/mask_genome.py --assembly $genome \
    --fasta data/$genome/hs1.fa.gz \
    --complement both --outdir refs/$genome --gzip
done

for genome in GRCm38 GRCm39; do
    python scripts/mask_genome.py --assembly $genome \
    --fasta data/$genome/Mus_musculus.$genome.dna.primary_assembly.fa.gz \
    --complement both --outdir refs/$genome --gzip --allow-approximate-mouse-par
done



# 3. mappability (single + multi-read) for every k, per flavour
for genome in GRCh38 GRCh37 T2T-CHM13v2 GRCm38 GRCm39; do
    for c in XX XY; do
        bash scripts/run_mappability.sh -f refs/$genome/$genome.$c.fa.gz \
    -s refs/$genome/$genome.$c.chrom.sizes -o map/$genome -l $genome.$c -t 16 
    done
done



# 4. callable BEDs + ploidy files (extra products)
for genome in GRCh38 GRCh37 T2T-CHM13v2 GRCm38 GRCm39; do
    for k in 24 36 50 100 150; do
        python scripts/build_callable_beds.py \
            --bigwig map/$genome/$genome.XY/k$k/$genome.XY.k$k.single_read.bw \
            --threshold 1.0 --out callable/$genome.XY.k$k.callable.bed
        python scripts/make_ploidy_configs.py --par-bed config/par/$genome.PAR.bed \
            --chrom-sizes refs/$genome/$genome.XY.chrom.sizes --label $genome --outdir config/ploidy
    done
done


# 5. (optional) territory tables — feed OUR sex-aware k100 multi_read bigWig,
#    NOT UCSC's precomputed track (that one is built on the UNMASKED genome and
#    would report mappability inside the PAR).
python scripts/build_territory_tables.py \
  --chrom-info data/GRCh38/chromInfo.txt.gz --gap data/GRCh38/gap.txt.gz \
  --cytoband data/GRCh38/cytoBand.txt.gz --blacklist data/GRCh38/ENCFF356LFX.hg38.blacklist.bed.gz \
  --umap-k100-bw map/GRCh38/GRCh38.XY/k100/GRCh38.XY.k100.multi_read.bw \
  --out-prefix territory/GRCh38.XY
```
