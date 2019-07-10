# Verify that x is a minimally valid non-negative matrix.
verify.matrix <- function (x, arg.name = deparse(substitute(x))) {
  arg.name <- sprintf("\"%s\"",arg.name)
  msg <- paste("Input argument",arg.name,"should be a non-negative,",
               "numeric matrix with at least two 2 rows and at least",
               "2 columns, and all entries should be finite and non-missing")
  if (!((is.matrix(x) & is.numeric(x)) |
        inherits(x,"dMatrix")))
    stop(msg)
  else if (nrow(x) < 2 |
           ncol(x) < 2 |
           any(x < 0)  |
           any(is.infinite(x)) |
           any(is.na(x)))
    stop(msg)
  return(TRUE)
}

# Scale each column A[,i] by b[i].
scale.cols <- function (A, b)
  t(t(A) * b)

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
  if (!is.list(fit))
    stop("Input argument \"fit\" should be a list")
  if (!all(is.element(c("F","L"),names(fit))))
    stop("Input argument \"fit\" should be a list containing named ",
         "elements \"F\" and \"L\"")
  F <- fit$F
  L <- fit$L

  # Verify and process input F.
  verify.matrix(F)
  if (!is.matrix(F))
    F <- as.matrix(F)
  if (is.integer(F))
    storage.mode(F) <- "double"

  # Verify and process input L.
  verify.matrix(L)
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
  
  # Compute the Poisson log-likelihood after removing terms that do
  # not depend on L or F.
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
