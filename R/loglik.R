#' @rdname altsqp
#' 
#' @export
#' 
loglik.poisson <- function (X, fit, e = 1e-15) {

  # Verify and process input matrix X. 
  verify.matrix(X)
  if (is.matrix(X) & is.integer(X))
    storage.mode(X) <- "double"
  
  # Input argument "fit" should be a list with elements "F" and "L".
  verify.fit(fit)
  F <- fit$F
  L <- fit$L

  # Verify and process input matrix F.
  verify.matrix(fit$F)
  if (!is.matrix(F))
    F <- as.matrix(F)
  if (is.integer(F))
    storage.mode(F) <- "double"

  # Verify and process input matrix L.
  verify.matrix(fit$L)
  if (!is.matrix(L))
    L <- as.matrix(L)
  if (is.integer(L))
    storage.mode(L) <- "double"

  # Check that matrices X, F and L are compatible.
  if (!(nrow(L) == nrow(X) &
        nrow(F) == ncol(X) &
        ncol(L) == ncol(F)))
    stop("Dimensions of input arguments \"X\", \"fit$F\" and/or \"fit$L\" ",
         "do not agree")
  
  # Compute the log-likelihood for the Poisson topic model after
  # removing terms that do not depend on L or F.
  return(-cost(X,L,t(F),e))
}

#' @rdname altsqp
#' 
#' @export
#' 
loglik.multinom <- function (X, fit, e = 1e-15) {

  # Verify and process input matrix X.
  verify.matrix(X)
  if (is.matrix(X) & is.integer(X))
    storage.mode(X) <- "double"
  
  # Input argument "fit" should be a list with elements "F" and "L".
  verify.fit(fit)
  F <- fit$F
  L <- fit$L

  # Verify and process input matrix F.
  verify.matrix(fit$F)
  eps <- 2*.Machine$double.eps
  if (any((colSums(F) - 1) > eps))
    stop("Each column of input matrix \"fit$F\" should sum to 1")
  if (!is.matrix(F))
    F <- as.matrix(F)
  if (is.integer(F))
    storage.mode(F) <- "double"

  # Verify and process input matrix L.
  verify.matrix(fit$L)
  if (any((rowSums(L) - 1) > eps))
    stop("Each row of input matrix \"fit$L\" should sum to 1")
  if (!is.matrix(L))
    L <- as.matrix(L)
  if (is.integer(L))
    storage.mode(L) <- "double"

  # Check that matrices X, F and L are compatible.
  if (!(nrow(L) == nrow(X) &
        nrow(F) == ncol(X) &
        ncol(L) == ncol(F)))
    stop(paste("Dimensions of input arguments \"X\", \"fit$F\" and/or",
               "\"fit$L\ do not agree"))
  
  # Compute the log-likelihood for the multinomial topic model after
  # removing terms that do not depend on L or F.
  sum(X * log(tcrossprod(L,F) + e))
}

# Compute the value of the cost function for non-negative matrix
# factorization, in which matrix X is approximated by matrix AB = A*B.
# This is equivalent to the negative Poisson log-likelihood after
# removing terms that do not depend on A or B.
cost <- function (X, A, B, e, version = c("R", "Rcpp")) {
  version <- match.arg(version)
  if (version == "Rcpp" & !is.matrix(X))
    stop("version = \"Rcpp\" is not yet implemented for sparse matrices")
  if (version == "R") {
    AB <- A %*% B
    f  <- sum(AB - X*log(AB + e))
  } else
    f <- cost_rcpp(X,A,B,e)
  return(f)
}
