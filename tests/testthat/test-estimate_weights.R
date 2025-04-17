#
# test_that("b contains the expected columns", {
#   load(test_path("../../data","b_estimate_dt.rda"))
#   expect_identical(colnames(b_estimate_dt),
#                c("total_reads","A", "C", "T",
#                  "G", "b_estimate", "ALT","CHROM_POS_REF_ALT"))
#
# })
#
# test_that("nrow(A) equals length(b)", {
#   load(test_path("../../data","b_estimate_dt.rda"))
#   load(test_path("../../data","A.rda"))
#
#   expect_equal(length(b_estimate_dt$b_estimate),
#                nrow(A))
# })
#
# test_that("Estimated weights have length equal to ncol(A)", {
#   load(test_path("../../data","b_estimate_dt.rda"))
#   load(test_path("../../data","A.rda"))
#   expect_equal(length(estimate_weights(b_estimate_dt$b_estimate,A)),
#                ncol(A)-1) # A is a data.table that contains donor names and one "rn" SNP id column
# })
#
# # fix main code so that slightly negative values are set to zero
# test_that("Estimated weights sum 1", {
#   load(test_path("../data","b_estimate_dt.rda"))
#   load(test_path("../data","A.rda"))
#   expect_equal(sum(estimate_weights(b_estimate_dt$b_estimate,A)), 1)
# })
#
# # w is not a named vector
