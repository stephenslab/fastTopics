# This function implements the alternating SQP ("alt-SQP") updates for
# the factors matrix, F, in which the matrix X is approximated by
# tcrossprod(L,F). The EM updates are equivalent to multiplicative
# updates, but computation is implemented differently. Inputs F and
# L should be dense matrices ("is.matrix" should return TRUE), but for
# X both dense matrices and sparse matrices are supported ("matrix"
# and "dgCMatrix" classes). Input "numiter" specifies the number of
# alt-SQP updates to perform.
#
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#'
altsqp_update_factors <-
  function (X, F, L,
            numiter = 1, control = fit_poisson_nmf_control_default()) {
  if (is.na(control$nc))
    control$nc <- 1
  F <- t(F)
  if (control$nc == 1) {
    if (is.matrix(X)) {
      F <- pnmfem_update_factors_rcpp(X,F,L,1)
      F <- altsqp_update_factors_rcpp(X,F,L,numiter,control)
    } else if (is.sparse.matrix(X)) {
      F <- pnmfem_update_factors_sparse_rcpp(X,F,L,1)
      F <- altsqp_update_factors_sparse_rcpp(X,F,L,numiter,control)
    }
  } else if (control$nc > 1) {
    if (is.matrix(X)) {
      F <- pnmfem_update_factors_parallel_rcpp(X,F,L,1)
      F <- altsqp_update_factors_parallel_rcpp(X,F,L,numiter,control)
    } else if (is.sparse.matrix(X)) {
      F <- pnmfem_update_factors_sparse_parallel_rcpp(X,F,L,1)
      F <- altsqp_update_factors_sparse_parallel_rcpp(X,F,L,numiter,control)
    }
  }
  return(t(F))
}

# This function implements the alternating SQP ("alt-SQP") updates for
# the loadings matrix, L, in which the matrix X is approximated by
# tcrossprod(L,F). The EM updates are equivalent to multiplicative
# updates, but computation is implemented differently. Inputs F and L
# should be dense matrices ("is.matrix" should return TRUE), but for X
# both dense matrices and sparse matrices are supported ("matrix" and
# "dgCMatrix" classes). Input "numiter" specifies the number of
# alt-SQP updates to perform.
#
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#'
altsqp_update_loadings <-
  function (X, F, L,
            numiter = 1, control = fit_poisson_nmf_control_default()) {
  if (is.na(control$nc))
    control$nc <- 1
  X <- t(X)
  L <- t(L)
  if (control$nc == 1) {
    if (is.matrix(X)) {
      L <- pnmfem_update_factors_rcpp(X,L,F,1)
      L <- altsqp_update_factors_rcpp(X,L,F,numiter,control)
    } else if (is.sparse.matrix(X)) {
      L <- pnmfem_update_factors_sparse_rcpp(X,L,F,1)
      L <- altsqp_update_factors_sparse_rcpp(X,L,F,numiter,control)
    }
  } else if (control$nc > 1) {
    if (is.matrix(X)) {
      L <- pnmfem_update_factors_parallel_rcpp(X,L,F,1)
      L <- altsqp_update_factors_parallel_rcpp(X,L,F,numiter,control)
    } else if (is.sparse.matrix(X)) {
      L <- pnmfem_update_factors_sparse_parallel_rcpp(X,L,F,1)
      L <- altsqp_update_factors_sparse_parallel_rcpp(X,L,F,numiter,control)
    }
  }
  return(t(L))
}
