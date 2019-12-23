#' @title Get Poisson Non-Negative Matrix Factorization from Multinomial
#'   Topic Model
#'
#' @description This function recovers parameter estimates of the
#'   Poisson non-negative matrix factorization given parameter estimates
#'   for a multinomial topic model. 
#'
#' @param fit A list containing two dense, non-negative matrices,
#'   \code{fit$F} and \code{fit$L}. The former is an m x k matrix of
#'   factors, and the latter is an n x k matrix of loadings. It does not
#'   make sense for a multinomial topic model to have less than two
#'   topics, so an error will be reported when k < 2. An additional
#'   vector \code{fit$s} of length n specifing the "document sizes" may
#'   be provided. If this is not provided, the document sizes will need
#'   to be estimated from the data (see input argument \code{X}).
#'
#' @param X The n x m matrix of counts or pseudocounts. It can be a
#'   sparse matrix (class \code{"dgCMatrix"}) or dense matrix (class
#'   \code{"matrix"}). This only needs to be provided if the document
#'   sizes \code{fit$s} are not made available.
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

  # Verify input argument "fit".
  verify.fit(fit)
  if (inherits(fit,"poisson_nmf"))
    stop("Input argument \"fit\" should not be a Poisson non-negative ",
         "matrix factorization")
  F <- fit$F
  L <- fit$L
  
  # Verify and process input matrix F.
  if (any(colSums(F) <= 0))
    stop("Each column of \"fit$F\" should have at least one positive entry")
  if (is.integer(F))
    storage.mode(F) <- "double"
  
  # Verify and process input matrix L.
  if (any(rowSums(L) <= 0))
    stop("Each row of \"fit$L\" should have at least one positive entry")
  if (is.integer(L))
    storage.mode(L) <- "double"
  
  # Check that k > 1.
  if (ncol(F) < 2)
    stop("Input matrices \"fit$F\" and \"fit$L\" should have at least ",
         "2 columns")

  # Exactly one of fit$s and X should be provided.
  if (sum(c(!missing(X),is.element("s",names(fit)))) != 1)
    stop("Exactly one of \"X\" and \"fit$s\" should be specified")
  
  if (missing(X)) {
      
    # Process the "scale factors", s.
    s <- fit$s
    s <- as.double(s)
  } else {

    # Check the data matrix, X.
    verify.fit.and.count.matrix(X,fit)
    if (is.matrix(X) & is.integer(X))
      storage.mode(X) <- "double"

    # Compute maximum-likelihood estimates of the "document sizes", s,
    # from the counts matrix, X.
    s <- rowSums(X)
  }

  # Recover F and L for the Poisson non-negative matrix factorization.
  out <- rescale.factors(F,s*L)
  
  # Update the "fit" object, and return it.
  fit$F <- out$F
  fit$L <- out$L
  fit$s <- NULL
  class(fit) <- c("poisson_nmf","list")
  return(fit)
}
