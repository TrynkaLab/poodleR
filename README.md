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


```r
install.packages("remotes")

remotes::install_github(
  "TrynkaLab/poodleR",
  dependencies = TRUE,
  build_vignettes = TRUE
)
```

## Publishing checklist

Before announcing the package publicly:

1. Install the package dependencies locally.
2. Run `devtools::check()` or `R CMD check`.
3. Confirm the vignette builds and the example workflow runs.
4. Push the repository to a public GitHub repo named `poodleR`.
5. Update the `USERNAME/poodleR` placeholders in this README.
6. Create a release tag such as `v0.1.0`.

This repository now includes a GitHub Actions workflow at `.github/workflows/R-CMD-check.yaml` so every push and pull request can run package checks automatically.

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
