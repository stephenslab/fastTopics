#' @title Get Poisson Non-Negative Matrix Factorization from Multinomial
#'   Topic Model
#'
#' @description This function recovers parameter estimates of the
#'   Poisson non-negative matrix factorization given parameter estimates
#'   for a multinomial topic model.
#'
#' @param fit Describe the "fit" argument here.
#'
#' @param X Describe the "X" argument here.
#'
#' @return Describe the return value here.
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

    # Check the vector of "scale factors", s.
    s <- fit$s
    # if (!(is.numeric(s) & length(s) == n))
    #   stop("Input \"fit$s\" should be a numeric vector of length n")
    # if (any(is.infinite(a)) | any(is.na(x)) | any(x <= 0))
    #   stop("Numeric vector \"fit$s\" should contain only finite, ",
    #        "non-missing and positive elements")        
    s <- as.double(s)
  } else {

    # Check the data matrix, X.
    verify.fit.and.count.matrix(X,fit)
    if (is.matrix(X) & is.integer(X))
      storage.mode(X) <- "double"

    # Compute the "size factors", s, from the counts matrix, X.
    # TO DO.
  }

  # Recover F and L for the Poisson non-negative matrix factorization.
  # L <- scale.cols(L,colSums(F))
  # s <- rowSums(L)
  # L <- L / s
  # F <- normalize.cols(F)
  
  # Update the "fit" object, and return it.
  fit$F <- F
  fit$L <- L
  fit$s <- NULL
  class(fit) <- c("poisson_nmf","list")
  return(fit)
}
