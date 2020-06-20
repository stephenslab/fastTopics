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
#' @import Matrix
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
  X <- matrix(as.double(rpois(n*m,tcrossprod(L,F))),n,m)
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

#' @title Add Title Here.
#'
#' @description Add description here.
#'
#' @param n Describe input argument "n" here.
#'
#' @param m Describe input argument "m" here.
#' 
#' @param k Describe input argument "k" here.
#' 
#' @importFrom stats rpois
#' @importFrom MCMCpack rdirichlet
#' 
#' @export
#'
simulate_poisson_gene_data <- function (n, m, k) {
 
  # Simulate the Poisson rates ("factors") according to the following
  # procedure:
  F <- matrix(0,m,k)
  for (j in 1:m) {
    # TO DO.
  }
    
  # For each sample, generate the topic proportions ("loadings")
  # according to the following procedure: (1) the number of nonzero
  # topic proportions is 1 <= n <= k with probability proportional to
  # 2^(-n); (2) sample the indices of the nonzero topics; (3) sample
  # the nonzero topic proportions from the Dirichlet distribution with
  # "prior sample sizes" alpha = 1, so that all topics are equally
  # likely.
  L  <- matrix(0,n,k)
  k1 <- sample(k,n,replace = TRUE,prob = 2^(-seq(1,k)))
  for (i in 1:n) {
    j      <- sample(k,k1[i])
    L[i,j] <- rdirichlet(1,rep(1,k1[i]))
  }
  
  # Add row and column names to outputs.
  rownames(L) <- paste0("i",1:n)
  rownames(F) <- paste0("j",1:m)
  colnames(F) <- paste0("k",1:k)
  colnames(L) <- paste0("k",1:k)
  
  # Return a list containing ...
  return(list(F = F,L = L))
}

#' @title Add Title Here.
#'
#' @description Add description here.
#' 
#' @export
#' 
simulate_multinom_gene_data <- function () {
  # TO DO.
} 
