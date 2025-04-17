#
# test_that("Test that the vcf data has vcfR format", {
#   load(test_path("../../data","vcf.rda"))
#
#   expect_identical(is(vcf), "vcfR")
# })
#
# test_that("Output of create_gt_dosages is a data.table", {
#   load(test_path("../../data","vcf.rda"))
#
#   expect_equal(sum(is(create_gt_dosages(vcf))== "data.table"),1)
# })
#
# test_that("Number of columns is number of donors + 1", {
#   load(test_path("../../data","vcf.rda"))
#
#   expect_equal(ncol(create_gt_dosages(vcf)),
#                ncol(vcf@gt[,setdiff(colnames(vcf@gt),"FORMAT")])+1) # checking with the variable part of the vcf@gt file
# })
#
# test_that("Output of create_gt_dosages contains a rn character column", {
#   load(test_path("../../data","vcf.rda"))
#
#   expect_equal(sum(is(create_gt_dosages(vcf)[[ setdiff(colnames(create_gt_dosages(vcf)),"rn")[1]]] )== "numeric"),
#                1) # testing with just one column - TODO test with all
# })
#
#
