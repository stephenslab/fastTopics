# This function implements a sequential co-ordinate descent (SCD)
# update for the factors matrix (the "basis vectors"), H, in which the
# matrix A is approximated by the matrix product W*H. Inputs W and H
# should be dense matrices ("is.matrix" should return TRUE). Input
# "numiter" specifies the number of inner-loop iterations to perform.
# Input argument "e" a non-negative scalar included in the
# computations to prevent NaNs due to division by zero.
#
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#'
scd_update_factors <- function (A, W, H, numiter = 1, nc = 1, e = 1e-16) {
  if (nc == 1) {
    if (is.matrix(A))
      H <- scd_update_factors_rcpp(A,W,H,numiter,e)
    else if (inherits(A,"dgCMatrix"))
      H <- scd_update_factors_sparse_rcpp(A,W,H,numiter,e)
  } else if (nc > 1) {
    if (is.matrix(A))
      H <- scd_update_factors_parallel_rcpp(A,W,H,numiter,e)
    else if (inherits(A,"dgCMatrix"))
      H <- scd_update_factors_sparse_parallel_rcpp(A,W,H,numiter,e)
  }  
  return(H)
}

# This function implements a sequential co-ordinate descent (SCD)
# update for the loadings matrix (the "activations"), W, in which the
# matrix A is approximated by the matrix product W*H. Inputs W and H
# should be dense matrices ("is.matrix" should return TRUE). Input
# "numiter" specifies the number of inner-loop iterations to perform.
# Input argument "e" a non-negative scalar included in the
# computations to prevent NaNs due to division by zero.
#
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#'
scd_update_loadings <- function (A, W, H, numiter = 1, nc = 1, e = 1e-16) {
  A <- t(A)
  W <- t(W)
  H <- t(H)
  if (nc == 1) {
    if (is.matrix(A))
      W <- scd_update_factors_rcpp(A,H,W,numiter,e)
    else if (inherits(A,"dgCMatrix"))
      W <- scd_update_factors_sparse_rcpp(A,H,W,numiter,e)
  } else if (nc > 1) {
    if (is.matrix(A))
      W <- scd_update_factors_parallel_rcpp(A,H,W,numiter,e)
    else if (inherits(A,"dgCMatrix"))
      W <- scd_update_factors_sparse_parallel_rcpp(A,H,W,numiter,e)
  }
  return(t(W))
}
