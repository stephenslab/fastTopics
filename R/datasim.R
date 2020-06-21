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
#' @details Note that only minimal argument checking is performed. This
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
  X <- generate_poisson_nmf_counts(F,L)
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

#' @title Simulate Gene-Expression Data with Absolute Expression Changes
#'
#' @description Simulate count data from a Poisson NMF model, in which
#' topics represent "gene expression programs", and the gene
#' expression programs are characterized by \emph{absolute} changes in
#' expression rates. The way in which the counts are simulated is
#' modeled after gene expression studies in which expression is
#' measured by single-cell RNA sequencing ("RNA-seq") techniques: each
#' row of the counts matrix corresponds a gene expression profile,
#' each column corresponds to a gene, and each matrix element is a
#' "read count" or "UMI count" measuring expression level. Factors are
#' simulated so as to capture "realistic" changes in gene expression
#' across different cell types. See \dQuote{Details} for the procedure
#' used to simulate factors, loadings and counts.
#'
#' @details Here we describe the process for generating the n x k
#'   loadings matrix \code{L} and the m x k factors matrix
#'   \code{F}.
#' 
#'   TO DO: Explain key distinguishing characteristic of the
#'   simulation process: factors can be interpreted as Poisson rates, so
#'   individual counts can be thought of as being generated from a
#'   weighted mixture of "topics" with (possibly) different expression
#'   levels. Therefore, the topics in this simulation are characterised
#'   by \emph{absolute} changes in gene expression levels.
#' 
#'   Once the loadings and factors have been generated, the
#'   counts are simulated from the Poisson NMF model; that is,
#'   \code{X[i,j]} is Poisson with rate \code{Y[i,j]}, where \code{Y =
#'   tcrossprod(L,F)}.
#'
#'   TO DO: Add material here.
#'
#'   Note that only minimal argument checking is performed;
#'   the function is mainly used to test implementation of the
#'   topic-modeling-based differential count analysis.
#'
#' @param n Number of rows in the simulated count matrix. It should be
#'   at least 2.
#'
#' @param m Number of columns in the simulated count matrix. It should
#'   be at least 2.
#'
#' @param k Number of factors, or "topics", used to generate the
#'   individual Poisson rates. The number of topics should be 2 or
#'   greater.
#'
#' @param sparse If \code{sparse = TRUE}, convert the counts matrix to
#'   a sparse matrix in compressed, column-oriented format; see
#'   \code{\link[Matrix]{sparseMatrix}}.
#' 
#' @importFrom methods as
#' 
#' @export
#'
simulate_poisson_gene_data <- function (n, m, k, sparse = FALSE) {

  # Check inputs.
  if (!(is.scalar(n) & all(n >= 2)))
    stop("Input argument \"n\" should be 2 or greater")
  if (!(is.scalar(m) & all(m >= 2)))
    stop("Input argument \"m\" should be 2 or greater")
  if (!(is.scalar(k) & all(k >= 2)))
    stop("Input argument \"k\" should be 2 or greater")

  # Simulate the data.
  F <- generate_poisson_rates(m,k)
  L <- generate_topic_proportions(n,k)
  X <- generate_poisson_nmf_counts(F,L)
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
#' @export
#' 
simulate_multinom_gene_data <- function () {
  # TO DO.
} 

# Generate an n x m counts matrix X such that X[i,j] is Poisson with
# rate Y[i,j], where Y = tcrossprod(L,F).
#
#' @importFrom stats rpois
generate_poisson_nmf_counts <- function (F, L) {
  n <- nrow(L)
  m <- nrow(F)
  return(matrix(as.double(rpois(n*m,tcrossprod(L,F))),n,m))
}

# For each sample (row), generate the topic proportions ("loadings")
# according to the following procedure: (1) the number of nonzero
# topic proportions is 1 <= n <= k with probability proportional to
# 2^(-n); (2) sample the indices of the nonzero topics uniformly at
# random; (3) sample the nonzero topic proportions from the Dirichlet
# distribution with alpha = 1 (so that all topics are equally likely).
#
#' @importFrom MCMCpack rdirichlet
generate_topic_proportions <- function (n, k) {
  L  <- matrix(0,n,k)
  k1 <- sample(k,n,replace = TRUE,prob = 2^(-seq(1,k)))
  for (i in 1:n) {
    j      <- sample(k,k1[i])
    L[i,j] <- rdirichlet(1,rep(1,k1[i]))
  }
  return(L)
}

# Simulate the Poisson rates (these are the factors in the Poisson NMF
# model) according to the following procedure. For each count: (1)
# generate u = abs(r) - 5, where r ~ N(0,2); (2) for each topic k,
# generate the Poisson rate as exp(max(t,-5)), where t ~ 0.95 *
# N(u,s/10) + 0.05 * N(u,s), and s = exp(-u/8).
#
#' @importFrom stats runif
#' @importFrom stats rnorm
generate_poisson_rates <- function (m, k) {
  F <- matrix(0,m,k)
  for (j in 1:m) {
    u <- abs(rnorm(1,0,2)) - 5
    s <- exp(-u/8)
    a <- runif(k)
    z <- (a >= 0.05) * rnorm(k,u,s/10) +
         (a < 0.05)  * rnorm(k,u,s)
    F[j,] <- exp(pmax(-5,z))
  }
  return(F)
}
