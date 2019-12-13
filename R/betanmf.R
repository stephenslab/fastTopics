#
# NOTES:
#  + No convergence checks are performed.
#  + Can be derived as an EM algorithm.
# 
#' @title Multiplicative update rules for non-negative matrix factorization
#'
#' @description This function decomposes the input counts matrix X =
#'   L*F' by nonnegative matrix factorization (NMF) based on the
#'   "divergence" criterion; equivalently, it optimizes the likelihood
#'   under a Poisson model of the data, X. It runs a specified number of
#'   multiplicative updates to fit the L and F matrices. Note that the
#'   multiplicative updates can also be derived---and hence
#'   interpreted---as an EM algorithm.
#'
#' @details This function is mainly for internal use, and should only
#'   be called directly if you really know what you are doing. In
#'   particular, only minimal argument checking is performed; if you are
#'   not careful, you will get poor results are errors that are
#'   difficult to interpret.
#'
#'   This implementation is adapted from the MATLAB code by Daichi
#'   Kitamura \url{http://d-kitamura.net}.
#'
#' @param X The n x m matrix of counts. All entries of X should be
#'   non-negative.
#'
#' @param F This is the initial estimate of the factors (also called
#'   "basis vectors"). It should be an m x k matrix, where m is the
#'   number of columns of X, and k is the rank of the matrix
#'   factorization. All entries of F should be non-negative.
#'
#' @param L This is the initial estimate of the loadings (also called
#'   "activations"). It should an n x k matrix, where 
#'  
#' @param numiter The number of multiplicative updates to run. It
#'   should an n x k matrix, where n is the number of rows of X, and k
#'   is the rank of the matrix factorization. All entries of L should be
#'   non-negative.
#'
#' @param e Describe e here.
#' 
#' @references Lee, D. D. and Seung, H. S. (2001). Algorithms for
#'   non-negative matrix factorization. In Advances in Neural
#'   Information Processing Systems \bold{13}, 556–562.
#'
#' @seealso fit_topics
#' 
#' @keywords internal
#' 
#' @export
#'
betanmf <- function (X, F, L, numiter = 1000, e = 1e-15) {

  # CHECK INPUTS
  # ------------
  # Perfom some very basic checks of the inputs.
  if (!(is.matrix(X) & is.matrix(F) & is.matrix(L)))
    stop("Input arguments \"X\", \"F\" and \"L\" should be matrices; ",
         "see help(matrix) for more information")

  # RESCALE INITIAL ESTIMATES
  # -------------------------
  # Re-scale the initial estimates of the factors and loadings so that
  # they are on same scale on average. This is intended to improve
  # numerical stability of the multiplicative updates.
  r <- sqrt(mean(B)/mean(A))
  A <- r*A
  B <- B/r

  out <- betanmf_helper(X,L,t(F),numiter,e)
  return(list(F = t(out$B),L = out$A))
}

# This implements the core part of the betanmf function.
betanmf_helper <- function (X, A, B, numiter, e) {

  for (i in 1:numiter) {

    # Update the loadings ("activations").
    A <- scale.cols(A * ((X / (A %*% B)) %*% t(B)),1/rowSums(B))
    A <- pmax(A,e)

    # Update the factors ("basis vectors").
    B <- B * (t(A) %*% (X / (A %*% B))) / colSums(A)
    B <- pmax(B,e)
  }

  # Re-scale the final estimates of the factors and loadings so that
  # they are on same scale on average.
  r <- sqrt(mean(B)/mean(A))
  A <- r*A
  B <- B/r

  return(list(A = A,B = B))
}
