# TO DO: Explain here what this function does, and how to use it.
#
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#'
altsqp_update_factors <- function (X, F, L, numiter = 1,
                                   control = altsqp_control_default())
  t(altsqp_update_factors_rcpp(X,t(F),L,numiter,control))

# TO DO: Explain here what this function does, and how to use it.
#
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#'
altsqp_update_loadings <- function (X, F, L, numiter = 1,
                                    control = altsqp_control_default())
  t(altsqp_update_factors_rcpp(t(X),t(L),F,numiter,control))

# These are the default settings used for running the alt-SQP updates.
altsqp_control_default <- function()
  mixsqp_control_default()
