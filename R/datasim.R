# Simulate an n x m matrix X such that X[i,j] is Poisson with rate
# (mean) Y[i,j], such that Y = tcrossprod(L,F), L is an n x k matrix,
# and F is an m x k matrix. The entries of matrix L are drawn
# uniformly at random between 0 and lmax, and the entries of matrix F
# are drawn uniformly at random between 0 and fmax. The return value
# is a list containing the data matrix X, and the factors used to
# simulate the data, L and F. This function is mainly used internally
# to simulate small data sets for testing the package functions.
#
#' @importFrom methods as
#' @importFrom stats runif
#' @importFrom stats rpois
#' 
simulate_poisson_data <- function (n, m, k, lmax = 1, fmax = 1,
                                   sparse = FALSE) {
  F <- matrix(runif(m*k,0,fmax),m,k)
  L <- matrix(runif(n*k,0,lmax),n,k)
  Y <- tcrossprod(L,F)
  X <- matrix(rpois(n*m,Y),n,m)
  if (sparse)
    X <- as(X,"dgCMatrix")
  return(list(X = X,F = F,L = L))
}
