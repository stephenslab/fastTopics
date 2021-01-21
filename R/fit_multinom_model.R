#' @title Fit Simple Multinomial Model
#'
#' @description Give short description here.
#'
#' @details The MLE has a closed-form solution; no iterative updates
#'   are needed.
#' 
#' @param cluster A factor specifying a grouping, or clustering, of
#'   the rows of \code{X}; e.g., the \dQuote{cluster} output from
#'   \code{\link[stats]{kmeans}}.
#'
#' @param X The n x m matrix of counts; all entries of X should be
#'   non-negative. It can be a sparse matrix (class \code{"dgCMatrix"})
#'   or dense matrix (class \code{"matrix"}), with some exceptions (see
#'   \sQuote{Details}).
#'
#' @param verbose This is passed as the \dQuote{verbose} argument in
#'   the call to \code{\link{init_poisson_nmf}}.
#' 
#' @param \dots Additional arguments passed to
#'   \code{\link{init_poisson_nmf}}.
#' 
#' @return A multinomial topic model fit.
#' 
#' @seealso \code{\link{fit_topic_model}}
#' 
#' @importFrom Matrix rowSums
#' @importFrom Matrix colSums
#' 
#' @export
#'
fit_multinom_model <- function (cluster, X, verbose = FALSE, ...) {

  # Check the input data matrix.
  verify.count.matrix(X)

  # If necessary, remove all-zero columns from the counts matrix.
  if (any(colSums(X > 0) == 0)) {
    i <- which(colSums(X > 0) >= 1)
    X <- X[,i]
    message(sprintf(paste("One or more columns of X are all zero; after",
                          "removing all-zero columns, %d columns will be",
                          "used for model fitting"),length(i)))
  }

  # Get the number of rows (n) and columns (m) of the data matrix,
  n <- nrow(X)
  m <- ncol(X)
  
  # Check the "cluster" input.
  if (!is.factor(cluster))
    cluster <- factor(cluster)
  if (length(cluster) != n)
    stop("Input argument \"cluster\" should have one entry for each row of ",
         "\"X\"")
  if (any(table(cluster) == 0))
    stop("Each level must appear at least once in factor \"cluster\"")
  
  # Initialize the loadings and factors matrices from the clustering:
  # L[i,j] = 1 if row i is assigned to cluster j, and L[i,j] = 0
  # otherwise. The maximum-likelihood estimates of the factors have a
  # closed-form solution in this case.
  k <- nlevels(cluster)
  F <- matrix(0,m,k)
  L <- matrix(0,n,k)
  rownames(L) <- rownames(X)
  rownames(F) <- colnames(X)
  colnames(L) <- levels(cluster)
  colnames(F) <- levels(cluster)
  for (j in levels(cluster)) {
    i      <- which(cluster == j)
    L[i,j] <- 1
    F[,j]  <- colSums(X[i,])/sum(L[i,j])
  }
  L <- rowSums(X) * L
  
  # Return a multinomial topic model fit.
  return(poisson2multinom(init_poisson_nmf(X,F = F,L = L,
                                           verbose = verbose,...)))
}
