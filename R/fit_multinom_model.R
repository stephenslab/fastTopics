#' @title Fit Simple Multinomial Model
#'
#' @description Fit a simple multinomial model for count data, in
#'   which each sample (\emph{i.e.}, a row of the data matrix \code{X})
#'   is assigned to a cluster. Under this simple multinomial model,
#'   \eqn{x_{ij}} assigned to cluster \eqn{k} is multinomial with sample
#'   size \eqn{s_i = x_{i1} + ... + x_{im}} and multinomial
#'   probabilities \eqn{p_{1k}, ..., p_{mk}}. This is a special case of
#'   the multinomial topic model in which all the mixture proportions
#'   are either 0 or 1. The maximum-likelihood estimates (MLEs) of the
#'   multinomial probabilities have a closed-form solution; no
#'   iterative algorithm is needed to fit this simple model.
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
#' @importFrom Matrix colSums
#' 
#' @export
#'
fit_multinom_model <- function (cluster, X,
                                verbose = c("none","detailed"), ...) {

  # Check the input data matrix.
  verify.count.matrix(X)

  # Check and process input argument "verbose"
  verbose <- match.arg(verbose)
  
  # If necessary, remove all-zero columns from the counts matrix.
  if (any.allzero.cols(X)) {
    X <- remove.allzero.cols(X)
    warning(sprintf(paste("One or more columns of X are all zero; after",
                          "removing all-zero columns, %d columns will be",
                          "used for model fitting"),ncol(X)))
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

  # Return a multinomial topic model fit.
  return(poisson2multinom(init_poisson_nmf(X,F = F,L = L,
                                           verbose = verbose,...)))
}
