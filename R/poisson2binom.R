#' @title Recover Binomial Topic Model Fit from Poisson NMF fit
#'
#' @description Add a brief description
#'
#' @details Describe the binomial topic model in more detail here.
#'  
#' @param X The n x m \dQuote{binary} matrix; all entries of X should be
#'   between 0 and 1. It can be a sparse matrix (class \code{"dgCMatrix"})
#'   or dense matrix (class \code{"matrix"})
#'
#' @param fit Describe input argument "fit" here.
#'
#' @param Describe input argument "numiter" here.
#' 
#' @param verbose Describe input argument "verbose" here.
#' 
#' @return Describe the output here.
#' 
#' @examples
#' # See the vignette for an example.
#' 
#' @export
#'
poisson2binom <- function (X, fit, numiter = 20, verbose = TRUE) {

  # Check input argument "fit".
  if (inherits(fit,"binom_topic_model_fit"))
    return(fit)
  if (!inherits(fit,"poisson_nmf_fit"))
    stop("Input argument \"fit\" should be an object of class ",
         "\"poisson_nmf_fit\"")
  verify.fit(fit)
  if (ncol(fit$F) < 2 | ncol(fit$L) < 2)
    stop("Input matrices \"fit$F\" and \"fit$L\" should have 2 or more",
         "columns")

  # Check and process input argument "X".
  if (!((is.numeric(X) & is.matrix(X)) | is.sparse.matrix(X)))
    stop("Input argument \"X\" should be a numeric matrix (a \"matrix\" or ",
         "a \"dgCMatrix\")")
  verify.fit.and.count.matrix(X,fit)
  if (is.matrix(X) & is.integer(X))
    storage.mode(X) <- "double"

  # Choose U = diag(u) such that L*U is closer to being a matrix of
  # topic proportions.
  ones <- matrix(1,n,1)
  L <- fit$L
  F <- fit$F
  u <- drop(coef(NNLM::nnlm(L,ones)))
  L <- scale.cols(L,u)
  F <- scale.cols(F,1/u)
}
