test_that("estimate_b_from_bam_readcount returns the expected aligned columns", {
  skip_if_not_installed("vcfR")
  skip_if_not_installed("data.table")

  load(test_path("../../data", "bam_readcount.rda"))
  load(test_path("../../data", "vcf.rda"))

  colnames(bam_readcount) <- c("CHROM", "POS", "REF", "TOTAL_READS", "DEL", "A", "C", "G", "T", "N")
  bam_readcount$CHROM_POS_REF <- paste(bam_readcount$CHROM, bam_readcount$POS, bam_readcount$REF, sep = "_")

  vcf_reduced <- data.table::as.data.table(vcf@fix[, 1:5])
  vcf_reduced$CHROM_POS_REF <- paste(vcf_reduced$CHROM, vcf_reduced$POS, vcf_reduced$REF, sep = "_")
  vcf_reduced <- vcf_reduced[ALT %in% c("A", "C", "G", "T")]
  vcf_reduced <- vcf_reduced[!duplicated(CHROM_POS_REF)]
  bam_readcount <- bam_readcount[!duplicated(CHROM_POS_REF)]
  vcf_reduced <- vcf_reduced[vcf_reduced$CHROM_POS_REF %in% bam_readcount$CHROM_POS_REF, ]
  bam_readcount <- bam_readcount[bam_readcount$CHROM_POS_REF %in% vcf_reduced$CHROM_POS_REF, ]

  result <- estimate_b_from_bam_readcount(bam_readcount, vcf_reduced)

  expect_identical(
    colnames(result),
    c("total_reads", "A", "C", "T", "G", "b_estimate", "ALT", "CHROM_POS_REF")
  )
  expect_equal(nrow(result), nrow(bam_readcount))
  expect_true(all(is.na(result$b_estimate) | (result$b_estimate >= 0 & result$b_estimate <= 1)))
})
