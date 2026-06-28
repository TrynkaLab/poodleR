# poodleR

`poodleR` helps estimate donor proportions in pooled sequencing experiments using:

- donor genotypes stored in a VCF
- pooled read counts exported by `bam-readcount`

The package covers the full workflow:

1. Convert donor genotypes to minor-allele dosages with `create_gt_dosages()`.
2. Convert pooled read counts to minor-allele frequency estimates with `estimate_b_from_bam_readcount()`.
3. Run the combined preprocessing step with `calculate_nonRef_vector_and_genotype_df()`.
4. Estimate donor weights with `estimate_weights()`.
5. Use `is_identical_genotype()` and `mod_lsqlincon()` as lower-level utilities.

## Install from GitHub

This package can be installed in any computer running R (was tested with R v. 4.4.1 and Ubuntu 16.04, and MacOS Sequoia 15.7.7), and requires only enough RAM to support the in-memory operations. Package dependencies are detailed and handled in the DESCRIPTION file.

To install, simply run:
```r
install.packages("remotes")

remotes::install_github(
  "TrynkaLab/poodleR",
  dependencies = TRUE,
  build_vignettes = TRUE
)
```

## Example files shipped with the package

The repository includes compressed example inputs in `inst/extdata`:

```r
bam_path <- system.file("extdata", "bam-readcount_example.tsv.gz", package = "poodleR")
vcf_path <- system.file("extdata", "vcf_example.vcf.gz", package = "poodleR")
```

It also includes smaller `.rda` example objects in `data/` for quick experimentation:

- `A`: genotype dosage table
- `b`: pooled allele-frequency table with `b_estimate`
- `b_estimate_dt`: filtered allele-frequency table aligned to `A`
- `bam_readcount`: imported bam-readcount example
- `vcf`: donor genotype VCF as a `vcfR` object

## Quick start

```r
library(poodleR)

bam_path <- system.file("extdata", "bam-readcount_example.tsv.gz", package = "poodleR")
vcf_path <- system.file("extdata", "vcf_example.vcf.gz", package = "poodleR")

prepared <- calculate_nonRef_vector_and_genotype_df(
  bam_readcount_path = bam_path,
  vcf_path = vcf_path,
  min_read_overlap = 5
)

weights <- estimate_weights(
  b = prepared$b_estimate_table$b_estimate,
  A = prepared$ma_dosages
)

weights
sum(weights)
```

## Tutorial for all exported functions

### `create_gt_dosages()`

Convert a `vcfR` object to a donor-by-SNP dosage table:

```r
data(vcf, package = "poodleR")

dosages <- create_gt_dosages(vcf)
head(dosages[, 1:4])
```

Each donor column is numeric with values `0`, `0.5`, or `1`, and `rn` stores the SNP identifier (`CHROM_POS_REF`).

### `estimate_b_from_bam_readcount()`

Convert imported bam-readcount rows to pooled allele frequencies after aligning them to VCF alleles:

```r
data(bam_readcount, package = "poodleR")
data(vcf, package = "poodleR")

colnames(bam_readcount) <- c("CHROM", "POS", "REF", "TOTAL_READS", "DEL", "A", "C", "G", "T", "N")
bam_readcount$CHROM_POS_REF <- paste(bam_readcount$CHROM, bam_readcount$POS, bam_readcount$REF, sep = "_")

vcf_reduced <- data.table::as.data.table(vcf@fix[, 1:5])
vcf_reduced$CHROM_POS_REF <- paste(vcf_reduced$CHROM, vcf_reduced$POS, vcf_reduced$REF, sep = "_")
vcf_reduced <- vcf_reduced[vcf_reduced$CHROM_POS_REF %in% bam_readcount$CHROM_POS_REF, ]
bam_readcount <- bam_readcount[bam_readcount$CHROM_POS_REF %in% vcf_reduced$CHROM_POS_REF, ]

b_estimate_table <- estimate_b_from_bam_readcount(bam_readcount, vcf_reduced)
head(b_estimate_table)
```

### `calculate_nonRef_vector_and_genotype_df()`

Run the full preprocessing pipeline directly from the shipped example files:

```r
prepared <- calculate_nonRef_vector_and_genotype_df(
  bam_readcount_path = bam_path,
  vcf_path = vcf_path,
  min_read_overlap = 5
)

names(prepared)
nrow(prepared$b_estimate_table)
nrow(prepared$ma_dosages)
```

This returns two aligned tables ready for deconvolution:

- `b_estimate_table`: pooled allele frequencies
- `ma_dosages`: donor genotype dosages

### `estimate_weights()`

Estimate donor contributions from the aligned objects:

```r
weights <- estimate_weights(
  b = prepared$b_estimate_table$b_estimate,
  A = prepared$ma_dosages
)

weights
sum(weights)
```

The output is a named numeric vector constrained to sum to 1.

## Learn more

For a longer worked example, see the vignette:

```r
vignette("how_to_use_poodleR", package = "poodleR")
```

## See also
The code to reproduce the results in "Integrated QTL mapping and CRISPR screening in pooled iPSC-derived microglia reveals genetic drivers of neurodegenerative risk" (https://www.biorxiv.org/content/10.1101/2025.08.18.670767v1) can be found in https://github.com/TrynkaLab/Perez_Alcantara_pooled_iMGL
