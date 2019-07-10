#' @rdname altsqp
#' 
#' @export
#' 
loglik.poisson <- function (X, fit, e = 1e-15) {

  # Verify and process input X.
  verify.matrix(X)
  if (is.matrix(X) & is.integer(X))
    storage.mode(X) <- "double"
  
  # Input argument "fit" should be a list with elements "F" and "L".
  verify.fit(fit)
  F <- fit$F
  L <- fit$L

  # Verify and process input F.
  verify.matrix(fit$F)
  if (!is.matrix(F))
    F <- as.matrix(F)
  if (is.integer(F))
    storage.mode(F) <- "double"

  # Verify and process input L.
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
  return(-cost(X,tcrossprod(L,F),e))
}

# Compute the log-likelihood for the multinomial topic model. Input X
# is an n x p matrix of counts, F is a p x K matrix of "factors", and
# L is an n x K matrix of "loadings". Example:
#
#   X <- matrix(0:19,4,5,byrow = TRUE)
#   L <- matrix(0:7,4,2,byrow = TRUE)
#   F <- matrix(0:9,5,2,byrow = TRUE)
#   loglik.multinom(X,F,L)
#
#' @rdname altsqp
#' 
#' @export
#' 
loglik.multinom <- function (X, fit, e = 1e-15) {

  # Verify and process input X.
  verify.matrix(X)
  if (is.matrix(X) & is.integer(X))
    storage.mode(X) <- "double"
  
  # Input argument "fit" should be a list with elements "F" and "L".
  verify.fit(fit)
  F <- fit$F
  L <- fit$L

  # Verify and process input F.
  verify.matrix(fit$F)
  if (!is.matrix(F))
    F <- as.matrix(F)
  if (is.integer(F))
    storage.mode(F) <- "double"

  # Verify and process input L.
  verify.matrix(fit$L)
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
cost <- function (X, AB, e)
  sum(AB - X*log(AB + e))
