# Release checklist (Zenodo)

Everything in this file was measured against the live records and the built
outputs, not assumed. Re-measure before the next release rather than trusting
the numbers here.

## The records

| record | concept DOI (cite this) | v1 version DOI |
|---|---|---|
| genomes + annotation | `10.5281/zenodo.21442231` | `10.5281/zenodo.21442232` |
| mappability | `10.5281/zenodo.21461059` | `10.5281/zenodo.21461060` |

The concept DOI always resolves to the newest version; a version DOI pins one
upload forever. `README.md` and `CITATION.cff` cite the **concept** DOIs. Adding
files to a published record needs Zenodo *New version*, which mints a new version
DOI and seeds the draft with a copy of the previous files.

## What v1 actually contains

Checked against the Zenodo API, because the upload plan was written from the
local tree and the two differ:

- **Genome record, 40 files** — exactly four per assembly per complement:
  `<ASM>.<C>.fa.gz`, `<ASM>.<C>.fa.md5`, `<ASM>.<C>.chrom.sizes`,
  `<ASM>.<C>.masked.bed`. **There is no `metadata.json` in v1.** Those files are
  *new* in v2, not replacements — do not plan around overwriting them.
- **Mappability record, 100 files** — 20 per assembly, for all five assemblies
  **including GRCh37**.

> **The local `map/` tree is missing GRCh37.** It holds 80 bigWigs for GRCh38,
> T2T-CHM13v2, GRCm39 and GRCm38; the published record has 100. So a
> `MANIFEST.tsv` built by walking the local tree **under-describes the
> mappability record by 20 files**. Either recompute GRCh37 mappability before
> generating the manifest, or state in the record description that the manifest
> covers the genome and annotation files only.

## v2 upload

New files:

- 27 sex-aware GTFs — 9 `(assembly, release)` rows × 3 variants, uploaded
  individually because users take one:
  `<ASM>.{XX,XY,XY.noYPAR}.ensembl-<REL>.gtf.gz`
- one `annotations_provenance.tar.gz` holding every `.gtf.gz.md5`,
  `.metadata.json`, `.dropped_genes.tsv` and `.t2g.tsv` — 100+ sidecars in the
  Zenodo file list would be unusable
- the 10 `<ASM>.<C>.metadata.json` for the FASTAs (new in v2)
- `MANIFEST.tsv`

Replacing:

- the 10 `<ASM>.<C>.fa.md5`, which were wrong in v1 (see the erratum below). The
  `.fa.gz` bytes themselves are unchanged — do not re-upload 8 GB of FASTA.

Metadata to set: data license **CC-BY-4.0** (the code stays MIT, see `LICENSE`),
and attribute [Ensembl](https://ensembl.org) as the source of all annotation
under its unrestricted release policy.

## The erratum, verbatim for the description

> The `*.fa.md5` sidecars in version 1 recorded the md5 of the *decompressed*
> stream while naming the `.fa.gz` file, so `md5sum -c` reported a mismatch on a
> perfectly good download. The FASTA content was never affected. Version 2
> corrects every sidecar to the md5 of the distributed bytes and adds
> `MANIFEST.tsv`, which is authoritative. `metadata.json` additionally carries
> `output_fasta_uncompressed_md5` for comparing content across providers, where
> the compressed bytes legitimately differ.

Do not fix this silently. Anyone who ran `md5sum -c` on v1 concluded their
download was corrupt and needs to be told it was not.

## Numbers to check a rebuild against

Deviating from these in a new Ensembl release is the signal that the PAR
annotation changed, not that the tool broke.

| assembly / release | XX dropped | XY-noYPAR dropped | off-genome |
|---|---|---|---|
| GRCh38 r116 | 672 | 64 (1 straddling: XGY2) | 0 |
| GRCh38 r110 | 601 | 49 (1 straddling: XGY2) | 0 |
| GRCh38 r98  | 522 | 0 | 0 |
| GRCh37 r87  | 495 | 0 | 0 |
| GRCh37 r75  | 495 | 0 | 5,772 |
| GRCm39 r116 | 2,052 | 6 | 0 |
| GRCm39 r110 | 1,633 | 4 | 0 |
| GRCm38 r102 | 1,569 | 5 | 0 |
| GRCm38 r98  | 1,569 | 5 | 0 |

Three of those `XY-noYPAR` files (GRCh38 r98, GRCh37 r87, GRCh37 r75) are
content-identical to their `XY` counterpart because the genebuild places no gene
on the Y PAR. That is expected and flagged as `identical_to_standard_xy` in the
metadata and in `MANIFEST.tsv` — see `docs/methods.md`.

"Content-identical" means every feature line is identical. The **files are not**:
each carries its own `#!sexaware-variant` / `#!sexaware-mask-bed` /
`#!sexaware-note` header, so the two have different md5s and appear as distinct
rows in the manifest. Do not "deduplicate" them before upload on the strength of
a checksum comparison, and do not treat the differing digests as a build error.

## Pre-publish gate

```bash
bash tests/test_make_sexaware_gtf.sh            # 52 assertions, < 2 s
bash tests/test_build_manifest.sh               # 38 assertions, < 2 s
scripts/build_manifest.sh -o MANIFEST.tsv       # ~25 min, 265 files, 47.2 GB

# Verify every row against the files on disk. Use the LAST column (the
# repo-relative path), not the first (the basename): published files live in
# refs/<ASM>/ and map/<ASM>/..., and basenames are not unique across them --
# `md5sum -c` on column 1 from any single directory silently fails to find
# almost everything, which reads as corruption rather than as a bad command.
awk -F'\t' 'NR>1 {print $3"  "$NF}' MANIFEST.tsv | md5sum -c --quiet && echo "all 265 verify"
```

`build_manifest.sh` prunes `*genmap_index/` from the walk by default: those are
GenMap build scratch (2.4 GB x 16) that is never published, and hashing them
costs ~40 GB of needless I/O. They are the reason to keep an eye on the
`kind` column -- **a release manifest should contain zero `kind=other` rows**;
anything landing there is a file nothing in the pipeline claims to produce.
Pass `--include-intermediates` only to audit a working tree.

Measured composition of the v2 manifest:

| kind | rows | |
|---|---|---|
| `mappability` | 80 | 4 assemblies x 2 complements x 5 k x 2 track types (the local tree lacks GRCh37 -- see the warning above) |
| `provenance` | 37 | 27 GTF + 10 genome `metadata.json` |
| `checksum_sidecar` | 37 | 27 `.gtf.gz.md5` + 10 `.fa.md5` |
| `annotation` | 27 | the GTFs themselves |
| `t2g` | 27 | |
| `dropped_genes` | 27 | |
| `genome` | 10 | |
| `mask_bed` | 10 | |
| `chrom_sizes` | 10 | |
| **total** | **265** | **0 rows of `kind=other`, 0 rows undescribed** |

Only 37 of those files have a `metadata.json` of their own. The rest inherit
from the product they belong to: the manifest derives each filename's
`<ASM>.<XX|XY>[.noYPAR][.<source>-<release>]` stem and reads that product's
JSON, so a `.fa.md5`, a `.chrom.sizes`, a `.t2g.tsv` and a mappability bigWig
all land fully described. Mappability rows keep the inherited assembly and
variant -- that is genuinely which reference they were computed on -- but their
tool column is reset to GenMap, because the JSON they inherit from was written
by `mask_genome.py` and says nothing about k or about GenMap.

A row with `-` in the `assembly` column is therefore a signal, not a normal
state: it means the filename did not parse as a product of a registered
assembly.

Per-release post-conditions, run on the real outputs and all passing as of
the v2 build (81 assertions over 27 GTFs):

- `XX`: zero features on the Y contig.
- `XY-noYPAR`: zero surviving features overlapping `<ASM>.XY.masked.bed`.
- every variant: the X contig's feature block is **byte-identical** to the source
  GTF's. This is the regression test for keying the drop set on
  `(seqname, gene_id)` rather than `gene_id`.
- every variant: its seqname set is a subset of the paired `chrom.sizes`.
- every `.gtf.gz.md5` verifies with `md5sum -c`.
- the drop set agrees gene-for-gene with an independent `awk` implementation of
  the same overlap predicate driven by the same `masked.bed`.
