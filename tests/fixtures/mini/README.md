# `mini` fixture

A hand-written miniature assembly + annotation that exercises every branch of
`make_sexaware_gtf.py` in under a second, with no download.

Geometry (contig lengths in `MINI.*.chrom.sizes`):

```
1                100000
X               3000000
Y               3000000
MT                16569     in the genome, absent from the GTF (expected, never fatal)
HSCHR1_1_CTG3       ----     in the GTF, ABSENT from the genome (the patch-GTF trap)
```

Y mask (`MINI.XY.masked.bed`, 0-based half-open, as mask_genome.py writes it):

```
PAR1_Y_masked   Y    10000 - 2781479
PAR2_Y_masked   Y  2900000 - 2990000
         NPY    Y  2781479 - 2900000   (not masked; must survive XY-noYPAR)
```

Each gene below is `gene` + `transcript` + 2 `exon` lines = 4 features.

| gene_id | seqname | span | why it is here |
|---|---|---|---|
| `ENSG_AUTO1` | 1 | 1000-2000 | plain autosomal control; must survive everything |
| `ENSG_OOR` | 1 | 99000-200000 | **end past the contig length** -> `--on-out-of-range` |
| `ENSG_SHARED` | X | 100000-120000 | **shares its gene_id with the Y copy below.** This is the GRCh37 regression: a drop set keyed on `gene_id` alone deletes this gene, which is the exact opposite of what an XY reference is for |
| `ENSG_SHARED` | Y | 100000-120000 | wholly inside PAR1 -> dropped in XX and XY-noYPAR |
| `ENSG_NONAME` | Y | 200000-210000 | inside PAR1, and carries **no `gene_name`** |
| `ENSG_PAB` | Y | 2780000-2790000 | **straddles the PAB** (PAR1 ends at 2781479). Neither it nor either exon is wholly contained, so `any` drops it as `y_par_overlap_partial` while `contained` keeps it |
| `ENSG_NPY1` | Y | 2800000-2820000 | male-specific region; **must survive** XY-noYPAR |
| `ENSG_PAR2G` | Y | 2910000-2930000 | wholly inside PAR2 -> dropped |
| `ENSG_AUTO2` | 1 | 50000-51000 | **reuses `ENSG_AUTO1`'s `transcript_id`**, so one transcript_id maps to two gene_ids. A t2g table has to pick one, mis-assigning that transcript's reads, so the filter reports the collision instead of collapsing it quietly |
| `ENSG_PATCHG` | HSCHR1_1_CTG3 | 1000-2000 | on a contig absent from the genome -> `--on-unknown-seqname` |

Note that the X and Y copies of `ENSG_SHARED` share a `transcript_id` too, but
*unambiguously* (same gene_id), so collapsing them is harmless and must not warn.

`ENSG_NPY1` also carries a `ref_gene_id` attribute: looking up `gene_id` must not
match the tail of `ref_gene_id`.
