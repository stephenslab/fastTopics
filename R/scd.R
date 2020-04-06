# This function implements a sequential co-ordinate descent (SCD)
# update for the factors matrix (the "basis vectors"), H, in which the
# matrix A is approximated by the matrix product W*H. Inputs W and H
# should be dense matrices ("is.matrix" should return TRUE). Input "j"
# specifies which columns of H to update; by default, all columns are
# updated. Input "numiter" specifies the number of inner-loop
# iterations to perform.  Input argument "e" a non-negative scalar
# included in the computations to prevent NaNs due to division by
# zero.
#
# Note that a single EM update of each factor is performed before
# running the CCD updates (unless runem = FALSE).
#
# Also note that the RcppParallel multithreading (specified by
# argument "nc") will only work correctly if the number of threads is
# set beforehand using RcppParallel::setThreadOptions.
#
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#'
scd_update_factors <- function (A, W, H, j = seq(1,ncol(A)), numiter = 1,
                                nc = 1, e = 1e-16, runem = TRUE) {
  if (!is.numeric(j))
    stop("Input argument \"j\" should be a numeric vector")
  if (nc == 1) {
    if (is.matrix(A)) {
      if (runem)
        H <- pnmfem_update_factors_rcpp(A,H,W,j-1,1)
      H <- scd_update_factors_rcpp(A,W,H,j-1,numiter,e)
    } else if (is.sparse.matrix(A)) {
      if (runem)
        H <- pnmfem_update_factors_sparse_rcpp(A,H,W,j-1,1)
      H <- scd_update_factors_sparse_rcpp(A,W,H,j-1,numiter,e)
    }
  } else if (nc > 1) {
    if (is.matrix(A)) {
      if (runem)
        H <- pnmfem_update_factors_parallel_rcpp(A,H,W,j-1,1)
      H <- scd_update_factors_parallel_rcpp(A,W,H,j-1,numiter,e)
    } else if (is.sparse.matrix(A)) {
      if (runem)
        H <- pnmfem_update_factors_sparse_parallel_rcpp(A,H,W,j-1,1)
      H <- scd_update_factors_sparse_parallel_rcpp(A,W,H,j-1,numiter,e)
    }
  }  
  return(H)
}

# This function implements a sequential co-ordinate descent (SCD)
# update for the loadings matrix (the "activations"), W, in which the
# matrix A is approximated by the matrix product W*H. Inputs W and H
# should be dense matrices ("is.matrix" should return TRUE). Input "i"
# specifies which rows of W to update; by default, all rows are
# updated. Input "numiter" specifies the number of inner-loop
# iterations to perform. Input argument "e" a non-negative scalar
# included in the computations to prevent NaNs due to division by
# zero.
#
# Note that a single EM update of the loadings is performed before
# running the CCD updates (unless runem = FALSE).
#
# Also note that the RcppParallel multithreading (specified by
# argument "nc") will only work correctly if the number of threads is
# set beforehand using RcppParallel::setThreadOptions.
#
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#'
scd_update_loadings <- function (A, W, H, i = seq(1,nrow(A)), numiter = 1,
                                 nc = 1, e = 1e-16, runem = TRUE) {
  if (!is.numeric(i))
    stop("Input argument \"i\" should be a numeric vector")
  A <- t(A)
  W <- t(W)
  H <- t(H)
  if (nc == 1) {
    if (is.matrix(A)) {
     if (runem)
        W <- pnmfem_update_factors_rcpp(A,W,H,i-1,1)
      W <- scd_update_factors_rcpp(A,H,W,i-1,numiter,e)
    } else if (is.sparse.matrix(A)) {
     if (runem)
        W <- pnmfem_update_factors_sparse_rcpp(A,W,H,i-1,1)
      W <- scd_update_factors_sparse_rcpp(A,H,W,i-1,numiter,e)
    }
  } else if (nc > 1) {
    if (is.matrix(A)) {
     if (runem)
        W <- pnmfem_update_factors_parallel_rcpp(A,W,H,i-1,1)
      W <- scd_update_factors_parallel_rcpp(A,H,W,i-1,numiter,e)
    } else if (is.sparse.matrix(A)) {
      if (runem)
        W <- pnmfem_update_factors_sparse_parallel_rcpp(A,W,H,i-1,1)
      W <- scd_update_factors_sparse_parallel_rcpp(A,H,W,i-1,numiter,e)
    }
  }
  return(t(W))
}
