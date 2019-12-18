#
# NOTES:
#  + No convergence checks are performed.
#  + Individual updates are very simple, and fast.
#  + But can be very slow to converge, particularly when X is sparse.
#  + Add "safeguard" to prevent factors or loadings from ever being
#    exactly zero---explain why this is done.
#  + 'MUUP' which uses multiplicative updates (this is the approach
#     popularized by Lee and Seung but initially proposed in 
#     Daube-Witherspoon, M. E., & Muehllehner, G. (1986). An iterative 
#     image space reconstruction algorithm suitable for volume ECT. IEEE 
#     Trans. Med. Imaging, 5.
# 
#' @title Multiplicative Update Rules for Non-negative Matrix Factorization
#'
#' @description This function decomposes the input matrix X = L*F' by
#'   nonnegative matrix factorization (NMF) based on the "divergence"
#'   criterion; equivalently, it optimizes the likelihood under a
#'   Poisson model of the count data, X. It runs a specified number of
#'   multiplicative updates (MU) to fit the L and F matrices. Note that
#'   the multiplicative updates can also be derived, and hence
#'   interpreted, as an expectation maximization (EM) algorithm.
#'
#' @details The multiplicative updates are very simple and
#'   fast. However, they can also be very slow to converge to a
#'   stationary point of the objective, particularly when the data are
#'   sparse.
#'
#'   This function is mainly for internal use, and should only
#'   be called directly if you really know what you are doing. In
#'   particular, only minimal argument checking is performed; if you are
#'   not careful, you will get poor results are errors that are
#'   difficult to interpret.
#'
#'   This implementation is adapted from the MATLAB code by Daichi
#'   Kitamura \url{http://d-kitamura.net}.
#'
#' @param X The n x m matrix of counts; all entries of X should be
#'   non-negative. In particular, sparse matrices are not accommodated
#'   in this implementation; \code{is.matrix(X)} must give \code{TRUE},
#'   otherwise an error will be thrown.
#'
#' @param F0 This is the initial estimate of the factors (also called
#'  "basis vectors"). It should be an m x k matrix, where m is the
#'   number of columns of X, and k > 1 is the rank of the matrix
#'   factorization, or the number of topics. All entries of F should be
#'   non-negative.
#'
#' @param L0 This is the initial estimate of the loadings (also called
#'   "activations"). It should an n x k matrix, where n is the number of
#'   rows of X, and k > 1 is the rank of the matrix factorization. All
#'   entries of L should be non-negative.
#'  
#' @param numiter The number of multiplicative updates to run.
#' 
#' @param e A small positive constant used to safeguard the
#'   multiplicative updates. The multiplicative updates are implemented
#'   as \code{F <- pmax(F1,e)} and \code{L <- pmax(L1,e)}, where
#'   \code{F1} and \code{L1} are the factors and loadings matrices
#'   obtained by applying a single multiplicative update rule.
#'
#' @param verbose Describe "verbose" input argument here.
#' 
#' @return A list containing the updated matrices F and L.
#' 
#' @references
#'
#'   Gillis, N. and Glineur, F. (2012). Accelerated multiplicative
#'   updates and hierarchical ALS algorithms for nonnegative matrix
#'   factorization. \emph{Neural Computation} \code{24}, 1085–1105. 
#' 
#'   Lee, D. D. and Seung, H. S. (2001). Algorithms for
#'   non-negative matrix factorization. In \emph{Advances in Neural
#'   Information Processing Systems} \bold{13}, 556–562.
#'
#' @seealso fit_topics
#'
#' @examples
#' 
#' # Add example here.
#' 
#' @keywords internal
#' 
#' @export
#'
betanmf <- function (X, F0, L0, numiter = 1000, e = 1e-15, verbose = TRUE) {

  # Get the number of rows (n) and columns (m) of data matrix, and get
  # the rank of the matrix factorization (k).
  n <- nrow(X)
  m <- ncol(X)
  k <- ncol(F)
    
  # CHECK INPUTS
  # ------------
  # Perfom some very basic checks of the inputs.
  if (!(is.matrix(X) & is.matrix(F) & is.matrix(L)))
    stop("Input arguments \"X\", \"F\" and \"L\" should be matrices; ",
         "see help(matrix) for more information")
  if (k < 2)
    stop("Matrix factorization should have rank at least 2")

  # INITIALIZE ESTIMATES
  # --------------------
  # Initialize the factors and loadings. To prevent the multiplicative
  # updates from getting "stuck", force the initial estimates to be
  # positive.
  F <- pmax(F0,e)
  L <- pmax(L0,e)
  
  # Re-scale the initial estimates of the factors and loadings so that
  # they are on same scale on average. This is intended to improve
  # numerical stability of the multiplicative updates.
  out <- rescale.factors(F,L)
  F   <- out$F
  L   <- out$L

  # Run the multiplicative updates.
  out <- betanmf_helper(X,L,t(F),numiter,e,verbose)

  # Output the updated estimates of the factors and loadings, and ...
  F           <- t(out$B)
  L           <- out$A
  dimnames(F) <- dimnames(F0)
  dimnames(L) <- dimnames(L0)
  return(list(F = F,L = L,progress = progress))
}

# This implements the core part of the betanmf function.
betanmf_helper <- function (X, A, B, numiter, e, verbose) {
  for (i in 1:numiter) {
    out <- betanmf_update(X,A,B,e)
    A   <- out$A
    B   <- out$B
  }
  return(list(A = A,B = B))
}

# Perform a single multiplicative update with safeguarding to promote
# convergence, followed by rescaling.
betanmf_update <- function (X, A, B, e) {
    
  # Update the loadings ("activations").
  A <- scale.cols(A * tcrossprod(X / (A %*% B),B),1/rowSums(B))
  A <- pmax(A,e)

  # Update the factors ("basis vectors").
  B <- B * crossprod(A,X / (A %*% B)) / colSums(A)
  B <- pmax(B,e)

  # Rescale the factors and loadings.
  out <- rescale.factors(t(B),A)
  return(list(A = out$L,B = t(out$F)))
}
