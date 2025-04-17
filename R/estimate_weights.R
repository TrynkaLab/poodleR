
#' Check for Identical Genotypes Across Samples
#'
#' This function determines whether all genotypes in each row of a matrix are identical.
#' It is a faster alternative to `apply(A[,..cols], MARGIN = 1, FUN = function(x) var(x) == 0)`,
#' leveraging row-wise variance calculations.
#'
#' @param x A numeric matrix where rows represent loci and columns represent individuals or samples.
#' @param ... Additional arguments passed to `rowMeans()`.
#'
#' @return A logical vector where each element is `TRUE` if all values in the corresponding row are identical, and `FALSE` otherwise.
#'
#' @examples
#' mat <- matrix(c(1,1,1, 2,2,2, 1,2,1), nrow = 3, byrow = TRUE)
#' is_identical_genotype(mat)  # Returns: TRUE, TRUE, FALSE
#'
#' @export
is_identical_genotype <- function(x, ...) {
  # rowVar() == 0, faster version of apply(A[,..cols], MARGIN = 1, FUN = function(x) var(x) == 0)
  (rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)) == 0
}

#' Modify Least Squares Regression to Handle Non-Positive Definite Matrices
#'
#' This function extends standard least squares regression to handle cases where
#' the quadratic programming matrix `Dmat` is not positive definite.
#'
#' The quadprog package requires `Dmat` to be positive definite for `solve.QP()`, but in cases where
#' the number of variants (m) is smaller than the number of donors (n), `Dmat` may not meet this requirement.
#' This function detects such cases and applies a near-positive definite approximation using `Matrix::nearPD()`
#' to ensure stability.
#'
#' @param C A numeric vector or matrix representing the coefficient matrix of the least squares problem.
#'          If a vector, it is converted to a matrix.
#' @param d A numeric vector representing the right-hand side of the equation.
#' @param A (Optional) A numeric matrix representing inequality constraints (size m x n).
#' @param b (Optional) A numeric vector representing inequality constraint thresholds (length m).
#' @param Aeq (Optional) A numeric matrix representing equality constraints (size 1 x n).
#' @param beq (Optional) A numeric vector representing equality constraint thresholds (length 1).
#' @param lb (Optional) A numeric vector specifying lower bounds for variables (length n).
#' @param ub (Optional) A numeric vector specifying upper bounds for variables (length n).
#'
#' @return A numeric vector of length `n`, representing the solution to the constrained least squares problem.
#'
#' @details
#' - Ensures that all input matrices and vectors have valid dimensions.
#' - If `Dmat` (computed as `t(C) %*% C`) is not positive definite, it is replaced with its nearest positive definite approximation.
#' - Uses `quadprog::solve.QP()` to solve the quadratic programming problem.
#' - Inequality constraints `A * x <= b` and equality constraints `Aeq * x = beq` are enforced if provided.
#' - Implements lower and upper bound constraints as additional linear inequalities.
#'
#' @examples
#' \dontrun{
#' C <- matrix(c(1, 2, 3, 4), nrow = 2)
#' d <- c(5, 6)
#' A <- matrix(c(1, 1, -1, 0), nrow = 2)
#' b <- c(1, -2)
#' result <- mod_lsqlincon(C, d, A, b)
#' }
#'
#' @importFrom quadprog solve.QP
#' @importFrom corpcor is.positive.definite
#' @importFrom Matrix nearPD
#' @export

mod_lsqlincon = function (C, d, A = NULL, b = NULL, Aeq = NULL, beq = NULL,
          lb = NULL, ub = NULL)
{
  if (!requireNamespace("quadprog", quietly = TRUE)) {
    stop("quadprog needed for this function to work. Please install it.",
         call. = FALSE)
  }
  stopifnot(is.numeric(C), is.numeric(d))
  if (is.null(A) && !is.null(b) || !is.null(A) && is.null(b))
    stop("If any, both 'A' and 'b' must be NULL.")
  if (is.null(Aeq) && !is.null(beq) || !is.null(Aeq) && is.null(beq))
    stop("If any, both 'Aeq' and 'beq' must be NULL.")
  if (!is.matrix(C))
    C <- matrix(C, 1)
  mc <- nrow(C)
  nc <- ncol(C)
  n <- nc
  if (length(d) != mc)
    stop("Dimensions of 'C' and 'd' do not fit.")
  if (is.null(A) && is.null(Aeq) && is.null(lb) && is.null(ub))
    return(qr.solve(C, d))
  if (!is.null(A)) {
    if (!is.matrix(A))
      A <- matrix(A, 1)
    ma <- nrow(A)
    na <- ncol(A)
    if (na != n)
      stop("Number of columns of 'A' does not fit with 'C'.")
    A <- -A
    b <- -b
  }
  if (is.null(Aeq)) {
    meq <- 0
  }
  else {
    if (!is.matrix(Aeq))
      Aeq <- matrix(Aeq, 1)
    meq <- nrow(Aeq)
    neq <- ncol(Aeq)
    if (neq != n)
      stop("Number of columns of 'Aeq' does not fit with 'C'.")
  }
  if (is.null(lb)) {
    diag_lb <- NULL
  }
  else {
    if (length(lb) == 1) {
      lb <- rep(lb, n)
    }
    else if (length(lb) != n) {
      stop("Length of 'lb' and dimensions of C do not fit.")
    }
    diag_lb <- diag(n)
  }
  if (is.null(ub)) {
    diag_ub <- NULL
  }
  else {
    if (length(ub) == 1) {
      ub <- rep(ub, n)
    }
    else if (length(ub) != n) {
      stop("Length of 'ub' and dimensions of C do not fit.")
    }
    diag_ub <- -diag(n)
    ub <- -ub
  }
  Dmat <- t(C) %*% C

  if(corpcor::is.positive.definite(as.matrix(Dmat))==FALSE){
    message("matrix D in quadratic function is not positive definite!")
    message("Computing positive definite matrix approximation ...")
    Dmat = Matrix::nearPD(as.matrix(Dmat))$mat
  }
  dvec <- t(C) %*% d
  Amat <- rbind(Aeq, A, diag_lb, diag_ub)
  bvec <- c(beq, b, lb, ub)
  rslt <- quadprog::solve.QP(Dmat, dvec, t(Amat), bvec, meq = meq)
  rslt$solution
}


#' Estimate Donor Proportions (Weights) from Minor Allele Frequencies and Genotype Data
#'
#' This function estimates donor proportions (weights) by solving the least squares regression problem `Aw = b`,
#' where `b` represents minor allele frequencies, and `A` is a genotype matrix containing donor minor allele dosages.
#'
#' The solution `w` represents the estimated donor proportions, subject to the constraints that the sum of weights equals 1 (`sum(w) = 1`)
#' and that all weights lie within the range `[0,1]`.
#'
#' @param b A numeric vector of length `k` representing minor allele frequencies, where `k` is the number of donors.
#' @param A A numeric matrix with `k+1` columns:
#'   - `k` named donor columns containing minor allele dosages (0 for homozygotes of the major allele,
#'     1 for homozygotes of the minor allele, 0.5 for heterozygotes).
#'   - An additional column named `"rn"`, which contains SNP identifiers.
#'   The number of rows of `A` must be equal to the length of `b`.
#' @param force_zero A logical value indicating whether to force small negative values in `w` to zero.
#'   Negative values close to zero may arise due to numerical imprecision. Defaults to `TRUE`.
#'
#' @return A numeric vector of length `k` representing the estimated donor proportions (weights), where:
#'   - The sum of all weights equals 1.
#'   - Each weight is constrained between 0 and 1.
#'
#' @details
#' - The function checks that `length(b) == nrow(A)`, ensuring consistency between the frequency vector and genotype data.
#' - It removes SNPs where all genotypes are identical (zero variance), as these do not contribute information.
#' - SNPs with missing values (`NA`) in `b` are also removed.
#' - The function then solves for `w` using the `mod_lsqlincon` function, enforcing the constraints:
#'   - `sum(w) = 1`
#'   - `0 <= w <= 1`
#' - If `force_zero = TRUE`, small negative values in `w` are set to zero.
#'
#' @examples
#' \dontrun{
#' b <- c(0.1, 0.2, 0.3, 0.4)  # Example minor allele frequencies
#' A <- data.frame(rn = c("SNP1", "SNP2", "SNP3", "SNP4"),
#'                 donor1 = c(0, 0.5, 1, 0.5),
#'                 donor2 = c(1, 0, 0.5, 0),
#'                 donor3 = c(0.5, 1, 0, 1))
#' w <- estimate_weights(b, A)
#' print(w)
#' }
#'
#' @importFrom data.table setDT
#' @export
estimate_weights <- function(b, A=gt, force_zero = TRUE){
  # check dimensions
  if(!identical(length(b),nrow(A))) stop("Error: length(b) != nrow(A).\n")

  # ideally check with snp names if they are identical()

  ######### remove rows from A where all genotypes are identical - those are where variance per row equals 0 -------
  # This should have ideally been done in a previous step to save time and memory at this point
  cols = setdiff(colnames(A),"rn")

  to_remove =  is_identical_genotype(A[,..cols])
  message(paste0("We are removing ", sum(to_remove), " identical rows out of ", nrow(A),
                 ": ", floor((sum(to_remove)/nrow(A))*100), "% removed." ))
  A = A[!to_remove,]

  # remove same SNPs from b
  b = b[!to_remove]

  ###### Removing SNPs with NAs - not sampled (not "sequenced") in the tested coverage for this trial -------
  to_remove = is.na(b)
  message(paste0("We are removing ", sum(to_remove), " unsampled SNPs out of the remaining ", length(b),
                 ": ", signif(sum(to_remove)/length(b),digits = 2)*100, "% removed." ))
  b = b[!to_remove]
  A = A[!to_remove,]
  if ( length(b) != nrow(A))
    stop("Number of rows of genotype matrix is not equal to the length of the minor allele frequency vector.")


  # solving for w
  #  Aeqw=beq  where Aeq  is a 1×k matrix of ones and be is a length m -vector, also of ones, contrains the sum of weights to 1
  #  lb=0 and ub=1  constraints our solution to wi∈[0,1]

  w=mod_lsqlincon(C=as.matrix(A[,..cols]),
              d=b,
              Aeq=matrix(1,ncol=length(cols),nrow=1),
              beq=1,
              lb=0,
              ub=1)
  if(force_zero == TRUE ) {
    w[w < 0] = 0
  }

  names(w) = cols
  return(w)
}
