#' @title Recover Poisson NMF Fit from Multinomial Topic Model Fit
#'
#' @description This function recovers parameter estimates of the
#'   Poisson non-negative matrix factorization (NMF) given parameter
#'   estimates for a multinomial topic model.
#'
#' @param fit An object of class \dQuote{multinom_topic_model_fit},
#'   such as an output from \code{poisson2multinom}.
#'
#' @param X Optional n x m matrix of counts, or pseudocounts. It can
#'   be a sparse matrix (class \code{"dgCMatrix"}) or dense matrix
#'   (class \code{"matrix"}). This only needs to be provided if the
#'   document sizes \code{fit$s} are not available.
#'
#' @return The return value is the list \code{fit}, in which matrices
#'   \code{fit$F} and \code{fit$L} specify the factors and loadings in
#'   the Poisson non-negative matrix factorization; specifically,
#'   the counts matrix is modeled by the low-rank matrix product
#'   \code{tcrossprod(fit$L,fit$F)}.
#'
#' @importFrom Matrix rowSums
#' 
#' @export
#'
multinom2poisson <- function (fit, X) {

  # Check input argument "fit".
  if (!inherits(fit,"multinom_topic_model_fit"))
    stop("Input argument \"fit\" should be an object of class ",
         "\"multinom_topic_model_fit\"")
  verify.fit(fit)
  F <- fit$F
  L <- fit$L

  # Check input argument "X".
  if (!missing(X))
    verify.fit.and.count.matrix(X,fit)    
  
  # Exactly one of X and fit$s should be provided.
  if (sum(c(!missing(X),is.element("s",names(fit)))) != 1)
    stop("Exactly one of \"X\" and \"fit$s\" should be specified")
  
  if (missing(X))
      
    # Process the "scale factors", s.
    s <- as.double(fit$s)
  else

    # Compute maximum-likelihood estimates of the "document sizes", s,
    # from the counts matrix, X.
    s <- as.double(rowSums(X))

  # Recover F and L for the Poisson non-negative matrix factorization.
  out <- rescale.factors(F,s*L)
  
  # Update the "fit" object, and return it.
  fit$F <- out$F
  fit$L <- out$L
  fit$s <- NULL
  class(fit) <- c("poisson_nmf_fit","list")
  return(fit)
}
