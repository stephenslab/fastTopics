#
# NOTES:
#  + No convergence checks are performed.
#  + Individual updates are very simple, and fast.
#  + Add "safeguard" to prevent factors or loadings from ever being
#    exactly zero---explain why this is done.
#  + Explain re-scaling step.
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
#' @param lowerbound A small positive constant used to safeguard the
#'   multiplicative updates. The multiplicative updates are implemented
#'   as \code{F <- pmax(F1,lowerbound)} and \code{L <- pmax(L1,lowerbound)},
#'   where \code{F1} and \code{L1} are the factors and loadings
#'   matrices obtained by applying a single multiplicative update
#'   rule. Setting \code{lowerbound = 0} is allowed, but the
#'   multiplicative updates are not guaranteed to converge to a
#'   stationary point without this safeguard, and a warning will be
#'   given in this case.
#'
#' @param e A small, non-negative number added to the terms inside the
#'   logarithms to avoid computing logarithms of zero. This prevents
#'   numerical problems at the cost of introducing a very small
#'   inaccuracy in the computation.
#' 
#' @param verbose When \code{verbose = TRUE}, information about the
#'   algorithm's progress is printed to the console at every iteration.
#'
#' @return A list containing updated estimates of the factors, F, and
#'   loadings, L.
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
#' # Simulate a 100 x 200 data set.
#' set.seed(1)
#' X <- simulate_count_data(100,200,3)$X
#' 
#' # Fit a Poisson topic model with k = 3 topics by running 100
#' # multiplicative updates.
#' F0  <- matrix(runif(600),200,3)
#' L0  <- matrix(runif(300),100,3)
#' fit <- betanmf(X,F0,L0,100)
#' 
#' # Plot the improvement in the solution over time.
#' dev.min <- 19661.4155
#' with(fit$progress,
#'      plot(iter,dev - dev.min,type = "l",log = "y",
#'           xlab = "iteration",ylab = "distance to solution"))
#' 
#' @export
#'
betanmf <- function (X, F0, L0, numiter = 1000, lowerbound = 1e-15,
                     e = 1e-15, verbose = TRUE) {

  # CHECK INPUTS
  # ------------
  # Perfom some very basic checks of the inputs.
  if (!(is.matrix(X) & is.matrix(F0) & is.matrix(L0)))
    stop("Input arguments \"X\", \"F0\" and \"L0\" should be matrices; ",
         "see help(matrix) for more information")

  # Get the number of rows (n) and columns (m) of data matrix, and get
  # the rank of the matrix factorization (k).
  n <- nrow(X)
  m <- ncol(X)
  k <- ncol(F0)
  if (k < 2)
    stop("Matrix factorization should have rank at least 2")

  # Check input argument "lowerbound".
  if (lowerbound < 0)
    stop("Input argument \"lowerbound\" should be zero or a positive number")
  if (lowerbound == 0)
    warning("Multiplicative updates may not converge when lowerbound = 0")
  
  # INITIALIZE ESTIMATES
  # --------------------
  # Initialize the estimates of the factors and loadings. To prevent
  # the multiplicative updates from getting "stuck", force the initial
  # estimates to be positive.
  F <- pmax(F0,lowerbound)
  L <- pmax(L0,lowerbound)
  
  # Re-scale the initial estimates of the factors and loadings so that
  # they are on same scale on average. This is intended to improve
  # numerical stability of the multiplicative updates.
  out <- rescale.factors(F,L)
  F   <- out$F
  L   <- out$L

  # Set up the data structure to record the algorithm's progress.
  progress <- as.matrix(data.frame(iter    = 1:numiter,
                                   loglik  = 0,
                                   dev     = 0,
                                   delta.f = 0,
                                   delta.l = 0,
                                   timing  = 0))
  
  # Run the multiplicative updates.
  if (verbose)
    cat("iter      log-likelihood            deviance max|F-F'| max|L-L'|\n")
  out <- betanmf_helper(X,L,t(F),lowerbound,e,progress,verbose)

  # Output the updated estimates of the factors and loadings, and ...
  F           <- t(out$B)
  L           <- out$A
  dimnames(F) <- dimnames(F0)
  dimnames(L) <- dimnames(L0)
  return(list(F = F,L = L,progress = out$progress))
}

# This implements the core part of the betanmf function.
betanmf_helper <- function (X, A, B, lowerbound, e, progress, verbose) {
  loglik.const <- loglik_poisson_const(X)
  dev.const    <- deviance_poisson_const(X)
  numiter      <- nrow(progress)
  for (i in 1:numiter) {
    A0     <- A
    B0     <- B
    timing <- system.time(out <- betanmf_update(X,A,B,lowerbound))
    A      <- out$A
    B      <- out$B
    progress[i,"timing"]  <- timing["elapsed"]
    progress[i,"loglik"]  <- sum(loglik.const - cost(X,A,B,e,"poisson"))
    progress[i,"dev"]     <- sum(dev.const + 2*cost(X,A,B,e,"poisson"))
    progress[i,"delta.f"] <- max(abs(B - B0))
    progress[i,"delta.l"] <- max(abs(A - A0))
    if (verbose)
      cat(sprintf("%4d %+0.12e %+0.12e %0.3e %0.3e\n",
                  i,progress[i,"loglik"],progress[i,"dev"],
                  progress[i,"delta.f"],progress[i,"delta.l"]))
  }
  return(list(A = A,B = B,progress = as.data.frame(progress)))
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

  # Rescale the factors and loadings, then output the updated
  # estimates.
  out <- rescale.factors(t(B),A)
  return(list(A = out$L,B = t(out$F)))
}
