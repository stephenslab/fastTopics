#' @title Simulate Count Data from Poisson NMF Model
#'
#' @description Simulate a counts matrix \code{X} such that
#'   \code{X[i,j]} is Poisson with rate (mean) \code{Y[i,j]}, where
#'   \code{Y = tcrossprod(L,F)}, \code{L} is an n x k loadings
#'   ("activations") matrix, and \code{F} is an m x k factors ("basis
#'   vectors") matrix. The entries of matrix \code{L} are drawn
#'   uniformly at random between 0 and \code{lmax}, and the entries of
#'   matrix \code{F} are drawn uniformly at random between 0 and
#'   \code{fmax}.
#'
#'   Note that only minimal argument checking is performed. This
#'   function is mainly used to simulate small data sets for the examples
#'   and package tests.
#' 
#' @param n Number of rows in simulated count matrix. The number of
#'   rows should be at least 2.
#'
#' @param m Number of columns in simulated count matrix. The number of
#'   columns should be at least 2.
#'
#' @param k Number of factors, or "topics", used to determine Poisson
#'   rates. The number of topics should be 1 or more.
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
#' @return The return value is a list containing the data matrix
#'   \code{X} and the factorization, \code{F} and \code{L}, used to
#'   determine the Poisson rates.
#'
#' @importFrom methods as
#' @importFrom stats rpois
#'
#' @export
#' 
simulate_count_data <- function (n, m, k, fmax = 1, lmax = 1, sparse = FALSE) {

  # Check inputs.
  if (!(is.scalar(n) & all(n >= 2)))
    stop("Input argument \"n\" should be 2 or more")
  if (!(is.scalar(m) & all(m >= 2)))
    stop("Input argument \"m\" should be 2 or more")
  if (!(is.scalar(k) & all(k >= 1)))
    stop("Input argument \"k\" should be 1 or more")
  if (!(is.scalar(fmax) & all(fmax > 0)))
    stop("Input argument \"fmax\" should be a positive number")
  if (!(is.scalar(lmax) & all(lmax > 0)))
    stop("Input argument \"lmax\" should be a positive number")
  
  # Simulate the data.
  F <- rand(m,k,0,fmax)
  L <- rand(n,k,0,lmax)
  X <- matrix(rpois(n*m,tcrossprod(L,F)),n,m)
  if (sparse)
    X <- as(X,"dgCMatrix")

  # Add row and column names to outputs.
  rownames(X) <- paste0("i",1:n)
  rownames(L) <- paste0("i",1:n)
  colnames(X) <- paste0("j",1:m)
  rownames(F) <- paste0("j",1:m)
  colnames(F) <- paste0("k",1:k)
  colnames(L) <- paste0("k",1:k)
  
  # Return a list containing: (1) the data matrix, X; (2) the factors
  # matrix, F; and (3) the loadings matrix, L.
  return(list(X = X,F = F,L = L))
}
