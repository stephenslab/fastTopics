#' @title Differential Count Analysis with a Multinomial Topic Model
#'
#' @description Describe function here.
#'
#' @details Add details here.
#' 
#' @param fit An object of class \dQuote{poisson_nmf_fit} or
#'   \dQuote{multinom_topic_model_fit}. If a Poisson NMF fit is provided
#'   as input, the corresponding multinomial topic model fit is
#'   automatically recovered using \code{link{poisson2multinom}}.
#'
#' @param X The n x m counts matrix. It can be a sparse matrix (class
#'   \code{"dgCMatrix"}) or dense matrix (class \code{"matrix"}).
#'
#' @param s A numeric vector of length n determining how the rates are
#'   scaled in the Poisson model. See \dQuote{Details} for guidance on
#'   the choice of \code{s}.
#' 
#' @param numiter The number of EM updates performed to compute the
#'   maximum-likelihood estimates of the Poisson model parameters.
#'
#' @param e A small, positive scalar included in some computations to
#'   avoid logarithms of zero and division by zero.
#'
#' @param verbose When \code{verbose = TRUE}, progress information is
#'   printed to the console.
#'
#' @return Describe the return value here.
#'
#' @importFrom Matrix rowSums
#' 
#' @export
#' 
diff_count_analysis <- function (fit, X, s = rowSums(X), numiter = 100,
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

  # Check and process input argument "X".
  if (!((is.numeric(X) & is.matrix(X)) | is.sparse.matrix(X)))
    stop("Input argument \"X\" should be a numeric matrix (a \"matrix\" or ",
         "a \"dgCMatrix\")")
  if (is.matrix(X) & is.integer(X))
    storage.mode(X) <- "double"

  # Check input argument "s".
  if (!(length(s) == nrow(X) & all(s > 0)))
    stop("Input argument \"s\" should be a vector of positive numbers, with ",
         "length(s) = nrow(X)")

  # Get the number of topics (k) and the number of rows in the counts
  # matrix (m).
  m <- nrow(fit$F)
  k <- ncol(fit$F)
  
  # Fit the univariate ("single-count") Poisson models.
  if (verbose)
    cat(sprintf("Fitting %d x %d = %d univariate Poisson models.\n",m,k,m*k))
  out <- fit_univar_poisson_models(X,fit$L,s,"em-rcpp",e,numiter,
                                   verbose = verbose)
  F0 <- out$F0
  F1 <- out$F1

  # Compute the log-fold change statistics, including standard errors
  # and z-scores.
  if (verbose)
    cat("Computing log-fold change statistics.\n")
  out <- compute_univar_poisson_zscores_fast(X,fit$L,F0,F1,s,e)

  # Return the Poisson model MLEs and the log-fold change statistics.
  out <- c(list(F0 = F0,F1 = F1),out)
  class(out) <- c("topic_model_diff_count","list")
  return(out)
}
