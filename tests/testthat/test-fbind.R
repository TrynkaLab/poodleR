test_that("factor concatenation gives factor not integer", {
  (a <- factor(c("character", "hits", "your", "eyeballs")))
  (b <- factor(c("but", "integer", "where it", "counts")))
  expect_true(is.factor(fbind(a,b)))
})
