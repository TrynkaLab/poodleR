test_that("create_gt_dosages returns numeric donor dosages with rn ids", {
  skip_if_not_installed("vcfR")
  skip_if_not_installed("data.table")

  data("vcf", package = "poodleR", envir = environment())

  dosage <- create_gt_dosages(vcf)
  donor_cols <- setdiff(colnames(dosage), "rn")

  expect_s3_class(dosage, "data.table")
  expect_true("rn" %in% colnames(dosage))
  expect_true(length(donor_cols) > 0)
  expect_true(all(vapply(dosage[, ..donor_cols], is.numeric, logical(1))))
  expect_false(anyNA(dosage[, ..donor_cols]))
})
