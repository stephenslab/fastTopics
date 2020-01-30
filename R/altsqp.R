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
  if (control$nc == 1) {
    if (is.matrix(X))
      F <- t(altsqp_update_factors_rcpp(X,t(F),L,numiter,control))
    else if (is.sparse.matrix(X))
      F <- t(altsqp_update_factors_sparse_rcpp(X,t(F),L,numiter,control))
  } else if (control$nc > 1) {
    if (is.matrix(X))
      F <- t(altsqp_update_factors_parallel_rcpp(X,t(F),L,numiter,control))
    else if (is.sparse.matrix(X))
      F <- t(altsqp_update_factors_sparse_parallel_rcpp(X,t(F),L,numiter,
                                                        control))
  }
  return(F)
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
  if (control$nc == 1) {
    if (is.matrix(X))
      L <- t(altsqp_update_factors_rcpp(t(X),t(L),F,numiter,control))
    else if (is.sparse.matrix(X))
      L <- t(altsqp_update_factors_sparse_rcpp(t(X),t(L),F,numiter,control))
  } else if (control$nc > 1) {
    if (is.matrix(X))
      L <- t(altsqp_update_factors_parallel_rcpp(t(X),t(L),F,numiter,control))
    else if (is.sparse.matrix(X))
      L <- t(altsqp_update_factors_sparse_parallel_rcpp(t(X),t(L),F,numiter,
                                                        control))
  }
  return(L)
}
