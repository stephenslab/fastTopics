# This function implements the EM updates for the factors matrix, F,
# in which the matrix X is approximated by tcrossprod(L,F). The EM
# updates are equivalent to multiplicative updates, but computation is
# implemented differently. Inputs F and L should be dense matrices
# ("is.matrix" should return TRUE), but for X both dense matrices and
# sparse matrices are supported ("matrix" and "dgCMatrix" classes).
# Input "j" specifies which rows of F to update; by default, all rows
# are updated. Input "numiter" specifies the number of EM updates to
# perform. Input argument "e" is a non-negative scalar specifying the
# minimum value of the updated loadings. A positive value of "e"
# promotes better convergence of the EM updates.
#
# Note that the RcppParallel multithreading (specified by argument
# "nc") will only work correctly if the number of threads is set
# beforehand using RcppParallel::setThreadOptions.
#
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#'
pnmfem_update_factors <- function (X, F, L, j = seq(1,ncol(X)),
                                   numiter = 1, nc = 1) {
  F <- t(F)
  if (nc == 1) {
    if (is.matrix(X))
      F <- pnmfem_update_factors_rcpp(X,F,L,j-1,numiter)
    else if (is.sparse.matrix(X))
      F <- pnmfem_update_factors_sparse_rcpp(X,F,L,j-1,numiter)
  } else if (nc > 1) {
    if (is.matrix(X))
      F <- pnmfem_update_factors_parallel_rcpp(X,F,L,j-1,numiter)
    else if (is.sparse.matrix(X))
      F <- pnmfem_update_factors_sparse_parallel_rcpp(X,F,L,j-1,numiter)
  }
  return(t(F))
}

# This function implements the EM updates for the loadings matrix, L,
# in which the matrix X is approximated by tcrossprod(L,F). The EM
# updates are equivalent to multiplicative updates, but computation is
# implemented differently. Inputs F and L should be dense matrices
# ("is.matrix" should return TRUE), but for X both dense matrices and
# sparse matrices are supported ("matrix" and "dgCMatrix" classes). 
# Input "i" specifies which rows of L to update; by default, all rows
# are updated. Input "numiter" specifies the number of EM updates to
# perform, and input "nc" specifies the number of threads to use in
# the multithreaded updates. Input argument "e" is a non-negative
# scalar specifying the minimum value of the updated loadings. A
# positive value of "e" promotes better convergence of the EM updates.
#
# Note that the RcppParallel multithreading (specified by argument
# "nc") will only work correctly if the number of threads is set
# beforehand using RcppParallel::setThreadOptions.
#
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#'
pnmfem_update_loadings <- function (X, F, L, i = seq(1,nrow(X)),
                                    numiter = 1, nc = 1) {
  X <- t(X)  
  L <- t(L)
  if (nc == 1) {
    if (is.matrix(X))
      L <- pnmfem_update_factors_rcpp(X,L,F,i-1,numiter)
    else if (is.sparse.matrix(X))
      L <- pnmfem_update_factors_sparse_rcpp(X,L,F,i-1,numiter)
  } else if (nc > 1) {
    if (is.matrix(X))
      L <- pnmfem_update_factors_parallel_rcpp(X,L,F,i-1,numiter)
    else if (is.sparse.matrix(X))
      L <- pnmfem_update_factors_sparse_parallel_rcpp(X,L,F,i-1,numiter)
  }
  return(t(L))
}
