

#' Mark rows where all columns (donors) in a genotype file have identical value
#'
#' Checks per row if all columns have identical numeric value and marks those rows that are identical
#' @param x Data table.
#' @param ... Need to check what this is.
#'
#' @return Vector of row names where genotypes are identical
#' @export
is_identical_genotype <- function(x, ...) {
  # rowVar() == 0, faster version of apply(A[,..cols], MARGIN = 1, FUN = function(x) var(x) == 0)
  (rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)) == 0
}

#' Modify main least square regression function to allow for approximations
#' for this error case:
#' >quadprog::solve.QP(Dmat, dvec, t(Amat), bvec, meq = meq) : 
#' "matrix D in quadratic function is not positive definite!"
#' This can happen when the matrix D is not positive definite, which is a requirement for the quadprog package,
#' and which may be caused by having fewer variants than donors.
#' That is not ideal, but this approximation will allow you to use this package for those extreme cases.
#' @param C a numeric vector of length n  (n = number of donors)
#' @param d a numeric vector of length n
#' @param A a numeric matrix of size m x n (m = number of variants, n = number of donors)
#' @param b a numeric vector of length m
#' @param Aeq a numeric matrix of size 1 x n
#' @param beq a numeric vector of length 1
#' @param lb a numeric vector of length n
#' @param ub a numeric vector of length n
#'
#' @return a numeric vector of length n
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


#' Estimate donor proportions (weights)
#'
#' Estimate donor proportions (weights) from the estimated MAF and the genotype
#' Solve for w in Aw = b using least square regression
#' With constraints sum(w) =1 and w can take values between 0-1.
#' @param b a vector of minor allele frequencies with length k (k = number of donors)
#' @param A a genotype file with k+1 columns. k named donor columns containing minor allele dosages (0 for homozygotes for major allele,
#' 1 for homozygotes of minor allele, 0.5 for heterozygotes).
#' Plus one extra column with "rn" name, containing SNP names. Number of rows of A must equal the length of b.
#' @param force_zero a logical. The least square regression result can return negative values that are very close to zero. If TRUE,
#' this parameter forces those values to zero. Defaults TRUE.
#'
#' @return a vector of proportions (w) with length k, where the sum equals 1.
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
                 ": the ", floor((sum(to_remove)/nrow(A))*100), "%" ))
  A = A[!to_remove,]

  # remove same SNPs from b
  b = b[!to_remove]

  ###### Removing SNPs with NAs - not sampled (not "sequenced") in the tested coverage for this trial -------
  to_remove = is.na(b)
  message(paste0("We are removing ", sum(to_remove), " unsampled SNPs out of the remaining ", length(b),
                 ": the ", signif(sum(to_remove)/length(b),digits = 2)*100, "%" ))
  b = b[!to_remove]
  A = A[!to_remove,]

  length(b) == nrow(A)

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
