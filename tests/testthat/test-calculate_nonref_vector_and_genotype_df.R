test_that("calculate_nonRef_vector_and_genotype_df returns aligned outputs", {
  skip_if_not_installed("data.table")
  skip_if_not_installed("R.utils")
  skip_if_not_installed("vcfR")
  skip_if_not_installed("stringr")
  skip_if_not_installed("dplyr")

  bam_path <- test_path("../../inst/extdata", "bam-readcount_example.tsv.gz")
  vcf_path <- test_path("../../inst/extdata", "vcf_example.vcf.gz")

  result <- calculate_nonRef_vector_and_genotype_df(
    bam_readcount_path = bam_path,
    vcf_path = vcf_path,
    min_read_overlap = 5
  )

  expect_named(result, c("b_estimate_table", "ma_dosages"))
  expect_s3_class(result$b_estimate_table, "data.table")
  expect_s3_class(result$ma_dosages, "data.table")
  expect_identical(result$b_estimate_table$CHROM_POS_REF, result$ma_dosages$rn)
  expect_true(all(is.na(result$b_estimate_table$b_estimate) | result$b_estimate_table$b_estimate >= 0))
})
