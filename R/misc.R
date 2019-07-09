# Scale each column A[,i] by b[i].
scale.cols <- function (A, b)
  t(t(A) * b)

#' @rdname altsqp
#' 
#' @export
#' 
loglik.poisson <- function (X, fit, e = 1e-15) {
  return(0)
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
  sum(X * log(tcrossprod(L,F) + e))
}

# Convert the parameters (factors & loadings) for the Poisson model to
# the factors and loadings for the multinomial model. The return value
# "s" gives the Poisson rates for generating the "document" sizes.
# Example:
#
#   L   <- matrix(0:7,4,2,byrow = TRUE)
#   F   <- matrix(0:9,5,2,byrow = TRUE)
#   out <- poisson2multinom(F,L)
#
#' @rdname altsqp
#' 
#' @export
#' 
poisson2multinom <- function (F, L) {
  L <- t(t(L) * colSums(F))
  s <- rowSums(L)
  L <- L / s
  F <- scale.cols(F)
  return(list(F = F,L = L,s = s))
}

# Compute the value of the cost function for non-negative matrix
# factorization, in which matrix X is approximated by matrix AB = A*B.
# This is equivalent to the negative Poisson log-likelihood after
# removing terms that do not depend on A or B.
cost <- function (X, AB, e)
  sum(AB - X*log(AB + e))
