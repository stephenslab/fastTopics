# This function implements the EM updates for the loadings matrix, L,
# in which the matrix X is approximated by tcrossprod(L,F). The EM
# updates are equivalent to multiplicative updates, but computation is
# implemented differently.
#
# Inputs L and F should be dense matrices ("is.matrix" should return
# TRUE), but for X both dense matrices and sparse matrices are
# supported ("matrix" and "dgCMatrix" classes).
# 
# Input "numiter" specifies the number of EM updates to perform. Input
# argument "e" is a non-negative scalar specifying the minimum value
# of the updated loadings. A positive value of "e" promotes better
# convergence of the multiplicative updates.
pnmfem_update_factors <- function (X, F, L, numiter = 1, nc = 1, e = 1e-15) {
  if (nc == 1) {
    if (is.matrix(X))
      F <- t(pnmfem_update_factors_rcpp(X,t(F),L,numiter))
    else if (inherits(X,"dgCMatrix"))
      F <- t(pnmfem_update_factors_sparse_rcpp(X,t(F),L,numiter))
  } else if (nc > 1) {
    if (is.matrix(X))
      F <- t(pnmfem_update_factors_parallel_rcpp(X,t(F),L,numiter))
    else if (inherits(X,"dgCMatrix"))
      F <- t(pnmfem_update_factors_sparse_parallel_rcpp(X,t(F),L,numiter))
  }
  return(pmax(F,e))
}

# This function implements the EM updates for the factors matrix, F, in
# which the matrix X is approximated by tcrossprod(L,F). The EM
# updates are equivalent to multiplicative updates, but computation is
# implemented differently.
#
# Inputs L and F should be dense matrices ("is.matrix" should return
# TRUE), but for X both dense matrices and sparse matrices are
# supported ("matrix" and "dgCMatrix" classes).
# 
# Input "numiter" specifies the number of EM updates to perform, and
# input "nc" specifies the number of threads to use in the
# multithreaded updates. Input argument "e" is a non-negative scalar
# specifying the minimum value of the updated loadings. A positive
# value of "e" promotes better convergence of the multiplicative
# updates.
pnmfem_update_loadings <- function (X, F, L, numiter = 1, nc = 1, e = 1e-15) {
  if (nc == 1) {
    if (is.matrix(X))
      L <- t(pnmfem_update_factors_rcpp(t(X),t(L),F,numiter))
    else if (inherits(X,"dgCMatrix"))
      L <- t(pnmfem_update_factors_sparse_rcpp(t(X),t(L),F,numiter))
  } else if (nc > 1) {
    if (is.matrix(X))
      L <- t(pnmfem_update_factors_parallel_rcpp(t(X),t(L),F,numiter))
    else if (inherits(X,"dgCMatrix"))
      L <- t(pnmfem_update_factors_sparse_parallel_rcpp(t(X),t(L),F,numiter))
  }
  return(pmax(L,e))
}
