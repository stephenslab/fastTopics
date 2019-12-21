# NOTES:
#
#   + EM updates are equivalent to multiplicative updates (see
#     betanmf), but computation is implemented differently.
#

#' @title Title goes here.
#'
#' @description Provide description here.
#'
#' @details Provide more details go here.
#'
#' @param X The n x m matrix of counts; all entries of X should be
#'   non-negative. Both dense and sparse matrices (the "matrix" and
#'   "dgCMatrix" classes) are supported.
#' 
#' @param F0 This is the initial estimate of the factors (also called
#'   "basis vectors"). It should be an m x k matrix, where m is the
#'   number of columns of X, and k > 1 is the rank of the matrix
#'   factorization, or the number of topics. All entries of F should be
#'   non-negative.
#'
#' @param L0 This is the initial estimate of the loadings (also called
#'   "activations"). It should an n x k matrix, where n is the number of
#'   rows of X, and k > 1 is the rank of the matrix factorization. All
#'   entries of L should be non-negative.
#'  
#' @param numiter The number of EM updates to perform.
#' 
#' @param minval A small positive constant used to safeguard the EM
#'   updates. The EM updates are implemented as \code{F <-
#'   pmax(F1,minval)} and \code{L <- pmax(L1,minval)}, where \code{F1}
#'   and \code{L1} are the factors and loadings matrices obtained by
#'   applying a single EM update. Setting \code{minval = 0} is allowed,
#'   but the EM updates are not guaranteed to converge to a stationary
#'   point without this safeguard, and a warning will be given in this
#'   case.
#'
#' @param e A small, non-negative number added to the terms inside the
#'   logarithms to avoid computing logarithms of zero. This prevents
#'   numerical problems at the cost of introducing a very small
#'   inaccuracy in the computation.
#' 
#' @param verbose When \code{verbose = TRUE}, information about the
#'   algorithm's progress is printed to the console at each iteration.
#'
#' @return \code{pnmfem} returns a list object with the following
#' elements:
#'
#' \item{F}{A dense matrix containing estimates of the factors.}
#'
#' \item{L}{A dense matrix containing estimates of the loadings.}
#'
#' \item{progress}{A data frame containing more detailed information
#'   about the algorithm's progress. The data frame should have
#'   \code{numiter} rows. The columns of the data frame are: "iter", the
#'   iteration number; "loglik", the log-likelihood at the current
#'   factor and loading estimates; "dev", the deviance at the current
#'   factor and loading estimates; "delta.l", the largest change in the
#'   factors matrix; "delta.f", the largest change in the loadings
#'   matrix; and "timing", the elapsed time in seconds (based on
#'   \code{\link{system.time}}).}
#' 
#' @references
#'
#' TO DO: Add reference(s).
#' 
#' @seealso \code{\link{betanmf}}
#' 
#' @examples
#'
#' # TO DO: Add example here.
#' 
#' @export
#' 
pnmfem <- function (X, L0, F0, numiter = 1000, minval = 1e-15,
                    e = 1e-15, verbose = TRUE) {

  # CHECK INPUTS
  # ------------
  # Perfom some very basic checks of the inputs.
  if (!(is.matrix(x) | inherits(x,"dgCMatrix")))
    stop("Input argument \"X\" should be a numeric matrix (a \"matrix\" or ",
         "a \"dgCMatrix\")")
  if (!(is.matrix(F0) & is.matrix(L0)))
    stop("Input arguments \"F0\" and \"L0\" should be numeric matrices; ",
         "see help(matrix) for more information")

  # Get the number of rows (n) and columns (m) of data matrix, and get
  # the rank of the matrix factorization (k).
  n <- nrow(X)
  m <- ncol(X)
  k <- ncol(F0)
  if (k < 2)
    stop("Matrix factorization should have rank at least 2")

  # Check input argument "minval".
  if (minval < 0)
    stop("Input argument \"minval\" should be zero or a positive number")
  if (minval == 0)
    warning("EM updates may not converge when \"minval\" is zero")
  
  # INITIALIZE ESTIMATES
  # --------------------
  # Initialize the estimates of the factors and loadings. To prevent
  # the EM updates from getting "stuck", force the initial estimates
  # to be positive.
  F <- pmax(F0,minval)
  L <- pmax(L0,minval)
  
  # Re-scale the initial estimates of the factors and loadings so that
  # they are on same scale on average. This is intended to improve
  # numerical stability of the EM updates.
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
  
  # Run the EM updates.
  if (verbose)
    cat("iter      log-likelihood            deviance max|F-F'| max|L-L'|\n")
  out <- pnmfem_helper(X,F,L,minval,e,progress,verbose)

  # Return a list containing (1) an estimate of the factors, (2) an
  # estimate of the loadings, and (3) a data frame recording the
  # algorithm's progress at each iteration.
  dimnames(out$F) <- dimnames(F0)
  dimnames(out$L) <- dimnames(L0)
  return(out)
}

# This implements the core part of the pnmfem function.
pnmfem_helper <- function (X, F, L, minval, e, progress, verbose) {
}
