# NOTES:
#
#   + EM updates are equivalent to multiplicative updates, but
#     computation is implemented differently.
#
#   + method = "scd" is based on version 0.4-3 of the NNLM package.
#
#   + method = "ccd" is based on MATLAB code by ...
#

#' @title Fit or Re-fit Poisson Non-negative Matrix Factorization
#'
#' @description This function decomposes the input matrix X = L*F' by
#'   nonnegative matrix factorization (NMF) based on the "divergence"
#'   criterion; equivalently, it optimizes the likelihood under a
#'   Poisson model of the count data, X. It runs a specified number of
#'   multiplicative updates (MU) or expectation maximization (EM)
#'   updates to fit the L and F matrices.
#'
#'   Although the EM updates are mathematically equivalent to the
#'   multiplicative updates, and therefore they share the same
#'   convergence properties, the implementation of EM is quite
#'   different; in particular, the EM updates are more suitable for
#'   sparse counts matrices.
#'
#' @details The multiplicative and EM updates are very simple and
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
#'   The implementation of the multiplicative updates is adapted from
#'   the MATLAB code by Daichi Kitamura \url{http://d-kitamura.net}.
#'
#'   The "safeguard" step preventing the factors and loadings from
#'   exactly reaching zero is motivated by Theorem 1 of Gillis & Glineur
#'   (2012).
#' 
#'   An additional re-scaling step is performed at each iteration to
#'   promote numerical stability.
#'
#'   Since the multiplicative updates are implemented using standard
#'   matrix operations, the speed is heavily dependent on the
#'   BLAS/LAPACK numerical libraries used. In particular, using
#'   optimized implementations such as OpenBLAS or Intel MKL can result
#'   in much improved performance of the multiplcative updates.
#'
#' @param X The n x m matrix of counts; all entries of X should be
#'   non-negative. It can be a sparse matrix (class \code{"dgCMatrix"})
#'   or dense matrix (class \code{"matrix"}), with some exceptions (see
#'   "Details").
#'
#' @param fit0 Describe input argument "fit" here.
#'  
#' @param k An integer 2 or greater giving the matrix rank for a
#'   random initialization of the factors and loadings. (They are
#'   initialized uniformly at random.) This argument should only be
#'   specified if the initial estimates (\code{fit} or \code{F, L})
#'   aren't already provided.
#' 
#' @param numiter The number of multiplicative updates to run.
#' 
#' @param method When \code{method = "em"}, the EM updates will be
#'   performed; when \code{method = "mu"}, the multiplicative updates
#'   will be performed. The multiplicative updates are only implemented
#'   for dense count matrices; if \code{method = "mu"} and \code{X} is a
#'   sparse matrix, and error will be generated.
#'
#' @param control Describe input argument "control" here.
#' 
#' @param minval A small positive constant used to safeguard the
#'   multiplicative updates. The multiplicative updates are implemented
#'   as \code{F <- pmax(F1,minval)} and \code{L <- pmax(L1,minval)},
#'   where \code{F1} and \code{L1} are the factors and loadings matrices
#'   obtained by applying a single multiplicative update rule. Setting
#'   \code{minval = 0} is allowed, but the multiplicative updates are
#'   not guaranteed to converge to a stationary point without this
#'   safeguard, and a warning will be given in this case.
#'
#' @param e A small, non-negative number added to the terms inside the
#'   logarithms to avoid computing logarithms of zero. This prevents
#'   numerical problems at the cost of introducing a very small
#'   inaccuracy in the computation.
#' 
#' @param verbose When \code{verbose = TRUE}, information about the
#'   algorithm's progress is printed to the console at each iteration.
#'
#' @return Both \code{init_poisson_nmf} and \code{fit_poisson_nmf}
#' return an object capturing the optimization algorithm state (for
#' \code{init_poisson_nmf}, this is the initial state). It is a list
#' with the following elements:
#'
#' \item{F}{A matrix containing the factor estimates.}
#'
#' \item{L}{A matrix containing estimates of the loadings.}
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
#'   Gillis, N. and Glineur, F. (2012). Accelerated multiplicative
#'   updates and hierarchical ALS algorithms for nonnegative matrix
#'   factorization. \emph{Neural Computation} \code{24}, 1085–1105. 
#' 
#'   Lee, D. D. and Seung, H. S. (2001). Algorithms for
#'   non-negative matrix factorization. In \emph{Advances in Neural
#'   Information Processing Systems} \bold{13}, 556–562.
#'
#' @seealso \code{\link{pnmfem}}
#'
#' @examples
#' 
#' # Simulate a 100 x 200 data set.
#' suppressWarnings(RNGversion("3.5.0"))
#' set.seed(1)
#' X <- simulate_count_data(100,200,3)$X
#' 
#' # Optimize a Poisson non-negative matrix factorization with k = 3
#' # topics by running 100 multiplicative updates.
#' F0   <- matrix(runif(600),200,3)
#' L0   <- matrix(runif(300),100,3)
#' fit1 <- betanmf(X,F0,L0,100,method = "mu")
#' fit2 <- betanmf(x,F0,L0.100,method = "em")
#' 
#' # Plot the improvement in the solution over time.
#' dev.min <- 19661.4155
#' with(fit$progress,
#'      plot(iter,dev - dev.min,type = "l",log = "y",
#'           xlab = "iteration",ylab = "distance to solution"))
#' 
#' @useDynLib fastTopics
#'
#' @importFrom utils modifyList
#' @importFrom RcppParallel setThreadOptions
#' @importFrom RcppParallel defaultNumThreads
#'
#' @export
#'
fit_poisson_nmf <- function (X, k, fit0, numiter = 100,
                             method = c("altsqp","scd","ccd","em","mu"), 
                             control = list(), verbose = TRUE) {

  # CHECK & PROCESS INPUTS
  # ----------------------
  # Check input argument "X".
  if (!(is.numeric(X) & (is.matrix(X) | is.sparse.matrix(X))))
    stop("Input argument \"X\" should be a numeric matrix (a \"matrix\" or ",
         "a \"dgCMatrix\")")
  if (is.matrix(X) & length(X) > 1e4 & mean(X > 0) < 0.1 & any(method != "mu"))
    message(paste("Input matrix \"X\" has less than 10% nonzero entries;",
                  "consider converting \"X\" to a sparse matrix to reduce",
                  "computational effort"))
  
  # Get the number of rows (n) and columns (m) of the data matrix,
  n <- nrow(X)
  m <- ncol(X)

  # Only one of "k" and "fit0" should be provided. If argument "k" is
  # given, generate a random initialization of the factors and
  # loadings.
  if (!(missing(k) & !missing(fit0) | (!missing(k) & missing(fit0))))
    stop("Provide a rank, \"k\", or an initialization, \"fit0\", but not both")
  if (missing(fit0))
    fit0 <- init_poisson_nmf(X,k = k)
  k <- ncol(fit0$F)

  # Check input argument "fit0".
  if (!(is.list(fit0) & inherits(fit0,"poisson_nmf_fit")))
    stop("Input argument \"fit0\" should be an object of class ",
         "\"poisson_nmf_fit\", such as an output of init_poisson_nmf")
  if (!(is.matrix(F) & is.numeric(F)))
  fit <- fit0
    
  # Check input argument "numiter".
  if (any(numiter < 1))
    stop("Input argument \"numiter\" must be 1 or greater")
      
  # Check and process input argument "method".
  method <- match.arg(method)
  if (method == "mu" & is.sparse.matrix(X)) {
    warning("method = \"mu\" is not implemented for sparse counts matrices; ",
            "attempting to convert \"X\" to a dense matrix")
    X <- as.matrix(X)
  }
  
  # Check and process the optimization settings.
  control <- modifyList(fit_poisson_nmf_control_default(),
                        control,keep.null = TRUE)
  if ((method == "mu" | method == "em") & any(control$minval == 0))
    warning("EM and multiplicative updates may not converge when minval = 0")
  if (is.na(control$nc)) {
    setThreadOptions()
    control$nc <- defaultNumThreads()
  } else
    setThreadOptions(numThreads = control$nc)
  if (control$nc > 1)
    message(sprintf("Setting number of RcppParallel threads to %d",control$nc))

  # INITIALIZE ESTIMATES
  # --------------------
  # Re-scale the initial estimates of the factors and loadings so that
  # they are on same scale on average. This is intended to improve
  # numerical stability of the multiplicative updates.
  out   <- rescale.factors(fit$F,fit$L)
  fit$F <- out$F
  fit$L <- out$L

  # Initialize the estimates of the factors and loadings. To prevent
  # the multiplicative updates from getting "stuck", force the initial
  # estimates to be positive.
  fit$F <- pmax(fit$F,control$minval)
  fit$L <- pmax(fit$L,control$minval)
  
  return(fit)
  
  # Set up the data structure to record the algorithm's progress.
  # progress <- as.matrix(data.frame(iter    = 1:numiter,
  #                                  loglik  = 0,
  #                                  dev     = 0,
  #                                  delta.f = 0,
  #                                  delta.l = 0,
  #                                  timing  = 0))
  
  # Run the multiplicative updates.
  # if (verbose)
  #   cat("iter      log-likelihood            deviance max|F-F'| max|L-L'|\n")
  # out <- betanmf_helper(X,F,L,method,minval,e,progress,verbose)

  # Return a list containing (1) an estimate of the factors, (2) an
  # estimate of the loadings, and (3) a data frame recording the
  # algorithm's progress at each iteration.
  F           <- t(out$B)
  L           <- out$A
  dimnames(F) <- dimnames(F0)
  dimnames(L) <- dimnames(L0)
  out         <- list(F = F,L = L,progress = out$progress)
  class(out)  <- c("poisson_nmf_fit","list")
  return(out)
}

#' @rdname fit_poisson_nmf
#'
#' @param F An optional argument giving is the initial estimate of the
#'   factors (also sometimes called the "basis vectors"). It should be
#'   an m x k matrix, where m is the number of columns in the counts
#'   matrix X, and k > 1 is the rank of the matrix factorization, or,
#'   equivalently, the number of topics. All entries of F should be
#'   non-negative. When not provided, input argument \code{k} should be
#'   specified.
#'
#' @param L An optional argument giving the initial estimate of the
#'   loadings (also sometimes called the "activations"). It should an n
#'   x k matrix, where n is the number of rows in the counts matrix X,
#'   and k > 1 is the rank of the matrix factorization. All entries of L
#'   should be non-negative. When not provided, input argument \code{k}
#'   should be specified.
#'
#' @export
#' 
init_poisson_nmf <- function (X, F, L, k) {

  # Check input X.
  if (!(is.numeric(X) & (is.matrix(X) | is.sparse.matrix(X))))
    stop("Input argument \"X\" should be a numeric matrix (a \"matrix\" or ",
         "a \"dgCMatrix\")")

  # Get the number of rows (n) and columns (m) of the data matrix,
  n <- nrow(X)
  m <- ncol(X)

  # Only one of k and (F,L) should be provided.
  if (!(missing(k) & (!missing(F) & !missing(L)) |
       (!missing(k) & (missing(F) & missing(L)))))
    stop("Provide a rank, k, or an initialization, (F, L), but not both")
  if (missing(k))
    k <- ncol(F)
  if (any(k < 2))
    stop("Matrix factorization rank \"k\" should be 2 or greater")

  # If the factor matrix is not provided, initialize the entries
  # uniformly at random.
  if (missing(F)) {
    F           <- rand(m,k)
    rownames(F) <- colnames(X)
    colnames(F) <- paste0("k",1:k)
  } else if (!(is.matrix(F) & is.numeric(F)))
    stop("Input argument \"F\" should be a numeric matrix (is.matrix(F) ",
         "should return TRUE)")
    
  # If the loading matrix is not provided, initialize the entries
  # uniformly at random.
  if (missing(L)) {
    L           <- rand(n,k)
    rownames(L) <- rownames(X)
    colnames(L) <- paste0("k",1:k)
  } else if (!(is.matrix(L) & is.numeric(L)))
    stop("Input argument \"L\" should be a numeric matrix (is.matrix(L) ",
         "should return TRUE)")

  # Return a list containing (1) the initial estimate of the factors,
  # and (2) the initial estimate of the loadings.
  fit        <- list(F = F,L = L)
  class(fit) <- c("poisson_nmf_fit","list")
  return(fit)
}

# This implements the core part of the betanmf function.
betanmf_main_loop <- function (X, F, L, method, minval, e, progress, verbose) {
  loglik.const <- loglik_poisson_const(X)
  dev.const    <- deviance_poisson_const(X)
  numiter      <- nrow(progress)
  for (i in 1:numiter) {
    F0     <- F
    L0     <- L
    timing <- system.time(out <- betanmf_update(X,A,B,minval),gcFirst = FALSE)
    A      <- out$A
    B      <- out$B
    progress[i,"timing"]  <- timing["elapsed"]
    progress[i,"loglik"]  <- sum(loglik.const - cost(X,A,B,e,"poisson"))
    progress[i,"dev"]     <- sum(dev.const + 2*cost(X,A,B,e,"poisson"))
    progress[i,"delta.f"] <- max(abs(F - F0))
    progress[i,"delta.l"] <- max(abs(L - L0))
    if (verbose)
      cat(sprintf("%4d %+0.12e %+0.12e %0.3e %0.3e\n",
                  i,progress[i,"loglik"],progress[i,"dev"],
                  progress[i,"delta.f"],progress[i,"delta.l"]))
  }
  return(list(F = F,L = L,progress = as.data.frame(progress)))
}

#' @rdname fit_poisson_nmf
#'
#' @export
#' 
fit_poisson_nmf_control_default <- function()
  c(list(num.updates = 4,
         minval      = 1e-15,
         e           = 1e-15,
         nc          = as.integer(NA)),
    mixsqp_control_default())
    
