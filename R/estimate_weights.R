

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

  w=pracma::lsqlincon(C=as.matrix(A[,..cols]),
              d=b,
              Aeq=matrix(1,ncol=length(cols),nrow=1),
              beq=1,
              lb=0,ub=1)
  if(force_zero == TRUE ) {
    w[w < 0] = 0
  }

  names(w) = cols
  return(w)
}
