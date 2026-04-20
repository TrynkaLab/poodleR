#' Example genotype dosage table
#'
#' A precomputed minor-allele dosage table with one row per SNP and one column
#' per donor, plus an `rn` SNP identifier column.
#'
#' @format A `data.table` with 8,598 rows and 20 columns.
"A"

#' Example minor-allele frequency vector as a data table
#'
#' A precomputed allele-frequency table derived from pooled read counts. The
#' `b_estimate` column is the input expected by [estimate_weights()].
#'
#' @format A `data.table` with 31,493 rows and 7 columns.
"b"

#' Example filtered bam-readcount allele-frequency table
#'
#' A filtered allele-frequency table aligned to the genotype dosage table. This
#' object includes the `CHROM_POS_REF` identifier used to match rows across
#' inputs.
#'
#' @format A `data.table` with 8,598 rows and 8 columns.
"b_estimate_dt"

#' Example bam-readcount import
#'
#' The first ten columns of a bam-readcount export after import with
#' `data.table::fread()`.
#'
#' @format A `data.table` with 10,000 rows and 10 columns.
"bam_readcount"

#' Example donor VCF object
#'
#' A small `vcfR` object used in examples and tests.
#'
#' @format A `vcfR` object.
"vcf"
