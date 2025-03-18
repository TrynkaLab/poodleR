test_that("bam_readcounts test file contains the expected columns", {
  expect_equal(ncol(bam_readcount), 10)
})

test_that("The returned data.table has the right columns ", {
  expect_identical(colnames(estimate_b_from_bam_readcount(bam_readcount,vcf)),
                   c("total_reads","A", "C", "T",
                     "G", "b_estimate", "ALT","C_POS_REF"))
  })
