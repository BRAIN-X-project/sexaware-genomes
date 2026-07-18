# Extra (best-practices) masking BEDs

These OPTIONAL BEDs implement the stricter masking from *"Best practices for
improving alignment and variant calling on human sex chromosomes"* (PMC12190741,
AJHG 2025), beyond the standard Olney-2020 PAR masking.

They target regions of **high but not perfect** X–Y homology that still cause
cross-mapping and false-positive variants:

- **XTR** — X-transposed region, ~3.4 Mb, Xq21 ↔ Yp11.2, ~98–99% identity.
- **Ampliconic / palindromic Y** — P1–P8 palindromes, AZF regions (internally
  repetitive, low mappability).
- **Gametologs** — single-copy Y genes with retained X homology.

Although these regions can be useful to mask for sex-specific aligments, 
**it is not always recommended for masking** based on previous literature 
(see `CITATION.cff` file).

## How to use
Two ways, per the paper:
1. **As downstream blacklists** (recommended default): keep the standard masked
   reference and simply *filter/flag* variants that fall in these BEDs.
2. **As extra hard-masking**: pass `--extra-mask <bed>` to `mask_genome.py` to N
   these regions in the FASTA too. This step is more aggressive, as you may lose real signal there.
