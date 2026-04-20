test_that("estimate_weights returns a named simplex solution", {
  skip_if_not_installed("data.table")
  skip_if_not_installed("quadprog")
  skip_if_not_installed("corpcor")

  A <- data.table::data.table(
    donor1 = c(0, 0.5, 1),
    donor2 = c(1, 0.5, 0),
    rn = c("snp1", "snp2", "snp3")
  )
  b <- c(0.2, 0.5, 0.8)

  weights <- estimate_weights(b = b, A = A)

  expect_named(weights, c("donor1", "donor2"))
  expect_equal(sum(weights), 1, tolerance = 1e-8)
  expect_true(all(weights >= 0))
  expect_true(all(weights <= 1))
})

test_that("is_identical_genotype flags invariant rows", {
  mat <- matrix(c(
    0, 0, 0,
    0, 0.5, 1,
    1, 1, 1
  ), byrow = TRUE, ncol = 3)

  expect_identical(is_identical_genotype(mat), c(TRUE, FALSE, TRUE))
})

test_that("mod_lsqlincon solves a bounded least-squares problem", {
  skip_if_not_installed("quadprog")
  skip_if_not_installed("corpcor")

  C <- matrix(c(
    0, 1,
    0.5, 0,
    1, 0.5
  ), byrow = TRUE, ncol = 2)
  d <- c(0.2, 0.5, 0.8)

  solution <- mod_lsqlincon(
    C = C,
    d = d,
    Aeq = matrix(1, nrow = 1, ncol = 2),
    beq = 1,
    lb = 0,
    ub = 1
  )

  expect_length(solution, 2)
  expect_equal(sum(solution), 1, tolerance = 1e-8)
  expect_true(all(solution >= 0))
  expect_true(all(solution <= 1))
})
