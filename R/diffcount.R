#' @title Differential Count Analysis with a Multinomial Topic Model
#'
#' @description Describe function here.
#'
#' @details Add details here.
#' 
#' @param fit Describe input argument "fit" here.
#'
#' @param X Describe input argument "X" here.
#'
#' @param s Describe input argument "s" here.
#' 
#' @param verbose Describe input argument "verbose" here.
#'
#' @param numiter Describe input argument "numiter" here.
#'
#' @param e Desribe input argument "e" here.
#'
#' @param verbose Describe input argument "verbose" here.
#' 
#' @return Describe the return value here.
#'
#' @importFrom Matrix rowSums
#' 
#' @export
#' 
diff_count_analysis <- function (fit, X, s = rowSums(X), numiter = 40,
                                 e = 1e-15, verbose = TRUE) {

  # CHECK & PROCESS INPUTS
  # ----------------------
  # Check and process input argument "fit".
  if (!(inherits(fit,"poisson_nmf_fit") |
        inherits(fit,"multinom_topic_model_fit")))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
         "\"multinom_topic_model_fit\"")
  if (inherits(fit,"poisson_nmf_fit"))
    fit <- poisson2multinom(fit)

  # Check input argument "X".
  # TO DO.
  
  # Fit the univariate ("single-count") Poisson models.
  if (verbose)
    cat("Fitting univariate Poisson models.\n")
  out <- fit_univar_poisson_models(X,fit$L,s,"em-rcpp",e,numiter,
                                   verbose = verbose)
  F0  <- out$F0
  F1  <- out$F1

  # Compute the log-fold change statistics, including standard errors
  # and z-scores.
  if (verbose)
    cat("Computing log-fold change z-scrores.\n")
  out <- compute_univar_poisson_zscores_fast(X,fit$L,F0,F1,s,e)

  # Return ...
  out$F0 <- F0
  out$F1 <- F1
  class(out) <- c("diff_count_analysis","list")
}
