#' @title Simulate count data
#'
#' @description Simulate an n x m matrix X such that X[i,j] is Poisson
#' with rate (mean) Y[i,j], such that Y = tcrossprod(L,F), L is an n x
#' k matrix, and F is an m x k matrix. The entries of matrix L are
#' drawn uniformly at random between 0 and lmax, and the entries of
#' matrix F are drawn uniformly at random between 0 and fmax. The
#' return value is a list containing the data matrix X, and the
#' factors used to simulate the data, L and F. This function is mainly
#' used internally to simulate small data sets for testing the package
#' functions.
#'
#' Only minimal argument checking is performed.
#' 
#' @param n Number of rows in simulated count matrix. 
#'
#' @param m Number of columns in simulated count matrix.
#'
#' @param k Number of factors used to determine Poisson rates.
#'
#' @param fmax Factors are drawn uniformly at random between 0 and
#'   \code{fmax}.
#'
#' @param lmax Loadings are drawn uniformly at random between 0 and
#'   \code{lmax}.
#'
#' @param sparse If \code{sparse = TRUE}, convert the counts matrix to
#'   a sparse matrix in compressed, column-oriented format; see
#'   \code{\link[Matrix]{sparseMatrix}}.
#' 
#' @return 
#'
#' @importFrom methods as
#' @importFrom stats runif
#' @importFrom stats rpois
#' 
simulate_count_data <- function (n, m, fmax = 1, lmax = 1, sparse = FALSE) {
  F <- matrix(runif(m*k,0,fmax),m,k)
  L <- matrix(runif(n*k,0,lmax),n,k)
  X <- matrix(rpois(n*m,tcrossprod(L,F)),n,m)
  if (sparse)
    X <- as(X,"dgCMatrix")

  #
  # TO DO: Add row and column names to X, F and L.
  #
  
  return(list(X = X,F = F,L = L))
}
