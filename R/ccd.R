# This function implements the cyclic co-ordinate descent (CCD) update
# for the factors matrix (the "basis vectors"), B, in which the matrix
# V is approximated by the matrix product W*H. Inputs W and H should
# be dense matrices ("is.matrix" should return TRUE). Input argument
# "e" is a non-negative scalar specifying the minimum value of the
# updated factors.
#
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#'
ccd_update_factors <- function (V, W, H, nc = 1, e = 1e-15) {
  if (nc == 1) {
    if (is.matrix(V))
      H <- ccd_update_factors_rcpp(V,W,H,e)
    else if (inherits(V,"dgCMatrix"))
      H <- ccd_update_factors_sparse_rcpp(V,W,H,e)
  } else if (nc > 1) {
    if (is.matrix(V)) {
      # H <- scd_update_factors_parallel_rcpp(A,W,H,numiter,e)
    } else if (inherits(V,"dgCMatrix")) {
      # H <- scd_update_factors_sparse_parallel_rcpp(A,W,H,numiter,e)
    }
  }  
  return(H)
}

# This function implements the cyclic co-ordinate descent (CCD) update
# for the loadings matrix (the "activations"), V, in which the matrix
# V is approximated by the matrix product W*H. Inputs W and H should
# not be sparse matrices ("is.matrix" should return TRUE). Input
# argument "e" is a non-negative scalar specifying the minimum value
# of the updated factors.
#
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#'
ccd_update_loadings <- function (V, W, H, nc = 1, e = 1e-15) {
  V <- t(V)
  W <- t(W)
  H <- t(H)
  if (nc == 1) {
    if (is.matrix(V))
      W <- ccd_update_factors_rcpp(V,H,W,e)
    else if (inherits(V,"dgCMatrix"))
      W <- ccd_update_factors_sparse_rcpp(V,H,W,e)
  } else if (nc > 1) {
    if (is.matrix(V)) {
      # TO DO.
    } else if (inherits(V,"dgCMatrix")) {
      # TO DO.
    }
  }
  return(t(W))
}
