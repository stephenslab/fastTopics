# NOTES:
#
#   + EM updates are equivalent to multiplicative updates, but
#     computation is implemented differently.
#
#   + method = "scd" is based on version 0.4-3 of the NNLM package.
#
#   + method = "ccd" is based on MATLAB code by ...
#
#   + "minval" is a small positive constant used to safeguard the
#     multiplicative updates. The multiplicative updates are implemented
#     as \code{F <- pmax(F1,minval)} and \code{L <- pmax(L1,minval)},
#     where \code{F1} and \code{L1} are the factors and loadings
#     matrices obtained by applying a single multiplicative update
#     rule. Setting \code{minval = 0} is allowed, but the multiplicative
#     updates are not guaranteed to converge to a stationary point
#     without this safeguard, and a warning will be given in this case.
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
#' \item{Fy}{A matrix containing the extrapolated factor estimates. If
#'   the extrapolation scheme is not used, \code{F} and \code{Fy} will
#'   be the same.}
#'
#' \item{Ly}{A matrix containing the extrapolated estimates of the
#'   loadings. If extrapolation is not used, \code{L} and \code{Ly} will
#'   be the same.}
#'
#' \item{Fbest}{A matrix containing the current best estimates of the
#'   factors. If extrapolation is not used, \code{F} and \code{Fbest}
#'   will be the same.}
#' 
#' \item{Lbest}{A matrix containing the current best estimates of the
#'   loadings. If extrapolation is not used, \code{L} and \code{Lbest}
#'   will be the same.}
#' 
#' \item{loss}{Value of objective ("loss") function computed at the
#'   extrapolated solution for the loadings (\code{Ly}) and the
#'   non-extrapolated solution for the factors (\code{F}). This is used
#'   internally to implement the extrapolated updates.}
#'
#' \item{lossbest}{Value of the objective ("loss" function computed at
#'   the current best estimates of the factors and loadings
#'   (\code{Fbest} and \code{Lbest}. If extrapolation is not used,
#'   \code{loss} and \code{lossbest} will be the same.}
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
#' @examples
#' # Simulate a 80 x 100 data set.
#' library(Matrix)
#' set.seed(1)
#' X <- simulate_count_data(80,100,3)$X
#'
#' # Run 20 EM updates to find a good initialization.
#' fit0 <- fit_poisson_nmf(X,k = 3,numiter = 20,verbose = FALSE)
#' 
#' # The fit_poisson_nmf interface implements 5 different algorithms 
#' # for optimizing a Poisson non-negative matrix factorization. Let's
#' # compare their runtime and abiliity to identify a good "fit".
#' fit.mu     <- fit_poisson_nmf(X,fit0 = fit0,numiter = 400,
#'                               method = "mu",verbose = FALSE)
#' fit.em     <- fit_poisson_nmf(X,fit0 = fit0,numiter = 400,
#'                               method = "em",verbose = FALSE)
#' fit.ccd    <- fit_poisson_nmf(X,fit0 = fit0,numiter = 300,
#'                               method = "ccd",verbose = FALSE)
#' fit.scd    <- fit_poisson_nmf(X,fit0 = fit0,numiter = 300,
#'                               method = "scd",verbose = FALSE)
#' fit.altsqp <- fit_poisson_nmf(X,fit0 = fit0,numiter = 200,
#'                               method = "altsqp",verbose = FALSE)
#'
#' clrs <- c("royalblue","skyblue","firebrick","orange","darkmagenta")
#' plot_progress_poisson_nmf(list(mu     = fit.mu,
#'                                em     = fit.em,
#'                                ccd    = fit.ccd,
#'                                scd    = fit.scd,
#'                                altsqp = fit.altsqp),
#'                           color = clrs)
#'
#' # All optimization algorithms other than the multiplicative updates
#' # can handle sparse matrices as well as dense ones.
#' Y <- as(X,"dgCMatrix")
#' fit.em.sp     <- fit_poisson_nmf(Y,fit0 = fit0,numiter = 400,
#'                                  method = "em",verbose = FALSE)
#' fit.ccd.sp    <- fit_poisson_nmf(Y,fit0 = fit0,numiter = 300,
#'                                  method = "ccd",verbose = FALSE)
#' fit.scd.sp    <- fit_poisson_nmf(Y,fit0 = fit0,numiter = 300,
#'                                  method = "scd",verbose = FALSE)
#' fit.altsqp.sp <- fit_poisson_nmf(Y,fit0 = fit0,numiter = 200,
#'                                  method = "altsqp",verbose = FALSE)
#' 
#' plot_progress_poisson_nmf(list(em        = fit.em,
#'                                ccd       = fit.ccd,
#'                                scd       = fit.scd,
#'                                altsqp    = fit.altsqp,
#'                                em.sp     = fit.em.sp,
#'                                ccd.sp    = fit.ccd.sp,
#'                                scd.sp    = fit.scd.sp,
#'                                altsqp.sp = fit.altsqp.sp),
#'                           color = rep(clrs[-1],times = 2),
#'                           shape = rep(c(20,18),each = 4))
#'
#' # The "extrapolated" updates can sometimes produce much better fits.
#' fit.ccd.ex <-
#'   fit_poisson_nmf(X,fit0 = fit0,numiter = 300,method = "ccd",
#'                   control = list(extrapolate = TRUE),verbose = FALSE)
#' fit.scd.ex <-
#'   fit_poisson_nmf(X,fit0 = fit0,numiter = 300,method = "scd",
#'                   control = list(extrapolate = TRUE),verbose = FALSE)
#' fit.altsqp.ex <-
#'   fit_poisson_nmf(X,fit0 = fit0,numiter = 200,method = "altsqp",
#'                   control = list(extrapolate = TRUE),verbose = FALSE)
#' 
#' plot_progress_poisson_nmf(list(ccd       = fit.ccd,
#'                                scd       = fit.scd,
#'                                altsqp    = fit.altsqp,
#'                                ccd.ex    = fit.ccd.ex,
#'                                scd.ex    = fit.scd.ex,
#'                                altsqp.ex = fit.altsqp.ex),
#'                                color = rep(clrs[3:5],times = 2),
#'                                shape = rep(c(20,18),each = 3))
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
                             method = c("em","altsqp","scd","ccd","mu"), 
                             control = list(), verbose = TRUE) {

  # CHECK & PROCESS INPUTS
  # ----------------------
  # Check input argument "X".
  if (!((is.numeric(X) & is.matrix(X)) | is.sparse.matrix(X)))
    stop("Input argument \"X\" should be a numeric matrix (a \"matrix\" or ",
         "a \"dgCMatrix\")")
  if (is.matrix(X) & length(X) > 1e4 & mean(X > 0) < 0.1 & any(method != "mu"))
    message(paste("Input matrix \"X\" has less than 10% nonzero entries;",
                  "consider converting \"X\" to a sparse matrix to reduce",
                  "computational effort"))
  if (is.matrix(X) & is.integer(X))
    storage.mode(X) <- "double"
  
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
    message(sprintf("Using %d RcppParallel threads.",control$nc))

  # Give an overview of the optimization settings.
  if (verbose) {
    cat(sprintf("Fitting rank-%d Poisson NMF to %d x %d %s matrix.\n",k,n,m,
                ifelse(is.matrix(X),"dense","sparse")))
    if (method == "mu")
      method.text <- "multiplicative"
    else if (method == "em")
      method.text <- "EM"
    else if (method == "scd")
      method.text <- "sequential co-ordinate descent (SCD)"
    else if (method == "ccd")
      method.text <- "cyclic co-ordinate descent (CCD)"
    else if (method == "altsqp")
      method.text <- "alternating SQP (alt-SQP)"
    cat(sprintf("Running %d %s updates, %s extrapolation ",numiter,
        method.text,ifelse(control$extrapolate,"with","without")))
    cat("(fastTopics 0.2-139).\n")
  }
  
  # INITIALIZE ESTIMATES
  # --------------------
  # Re-scale the initial estimates of the factors and loadings so that
  # they are on same scale on average. This is intended to improve
  # numerical stability of the multiplicative updates.
  out   <- rescale.factors(fit$F,fit$L)
  fit$F <- out$F
  fit$L <- out$L

  # Initialize the estimates of the factors and loadings. To prevent
  # the updates from getting "stuck", force the initial estimates to
  # be positive.
  fit$F <- pmax(fit$F,control$minval)
  fit$L <- pmax(fit$L,control$minval)

  # RUN UPDATES
  # -----------
  if (verbose)
    cat("iter log-likelihood      deviance res(KKT)",
        "max|F-F'| max|L-L'| nz(F) nz(L) beta\n")
  fit <- fit_poisson_nmf_main_loop(X,fit,numiter,method,control,verbose)

  # Output the updated "fit".
  fit$progress        <- rbind(fit0$progress,fit$progress)
  dimnames(fit$F)     <- dimnames(fit0$F)
  dimnames(fit$L)     <- dimnames(fit0$L)
  dimnames(fit$Fy)    <- dimnames(fit0$Fy)
  dimnames(fit$Ly)    <- dimnames(fit0$Ly)
  dimnames(fit$Fbest) <- dimnames(fit0$Fbest)
  dimnames(fit$Lbest) <- dimnames(fit0$Lbest)
  class(fit) <- c("poisson_nmf_fit","list")
  return(fit)
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
#' @param beta Initial setting of the extrapolation parameter. This is
#'   \eqn{beta} in Algorithm 3 of Ang & Gillis (2019).
#'
#' @param betamax Initial setting for the upper bound on the
#'   extrapolation parameter. This is \eqn{\bar{\gamma}} in Algorithm 3
#'   of Ang & Gillis (2019).
#'
#' @param e A small, non-negative number added to the terms inside the
#'   logarithms to avoid computing logarithms of zero. This prevents
#'   numerical problems at the cost of introducing a very small
#'   inaccuracy in the computation.
#' 
#' @export
#' 
init_poisson_nmf <- function (X, F, L, k, beta = 0.5, betamax = 0.99,
                              e = fit_poisson_nmf_control_default()$eps) {

  # Check input X.
  if (!((is.numeric(X) & is.matrix(X)) | is.sparse.matrix(X)))
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
  if (is.integer(F))
    storage.mode(F) <- "double"
    
  # If the loading matrix is not provided, initialize the entries
  # uniformly at random.
  if (missing(L)) {
    L           <- rand(n,k)
    rownames(L) <- rownames(X)
    colnames(L) <- paste0("k",1:k)
  } else if (!(is.matrix(L) & is.numeric(L)))
    stop("Input argument \"L\" should be a numeric matrix (is.matrix(L) ",
         "should return TRUE)")
  if (is.integer(L))
    storage.mode(L) <- "double"

  # Compute the value of the objective ("loss") function at the
  # initial estimates of the factors and loading.
  loss <- sum(cost(X,L,t(F),e))

  # Initialize the data frame for keeping track of the algorithm's
  # progress over time.
  progress        <- as.data.frame(matrix(0,0,11))
  names(progress) <- c("loglik","dev","res","delta.f","delta.l","nonzeros.f",
                       "nonzeros.l","extrapolate","beta","betamax","timing")
  
  # Return a list containing: F, an initial estimate of the factors;
  # L, an initial estimate of the loadings; Fy and Ly, the
  # extrapolated estimates of the factors and loadings, which
  # initially are always the same as F and L; Fbest and Lbest, the
  # current best estimates of the factors and loadings, which
  # initially are always the same as F and L; "loss" and "lossbest",
  # the value of the objective or loss function at the initial
  # estimates of the factors and loadings; beta, beta0 and betamax,
  # the initial settings of the extrapolation parameters; and
  # "progress", the initial data frame for keeping track of the
  # algorithm's progress over time.
  fit <- list(F        = F,
              L        = L,
              Fy       = F,
              Ly       = L,
              Fbest    = F,
              Lbest    = L,
              loss     = loss,
              lossbest = loss,
              beta     = beta,
              beta0    = beta,
              betamax  = betamax,
              progress = progress)
  class(fit) <- c("poisson_nmf_fit","list")
  return(fit)
}

# This implements the core part of fit_poisson_nmf.
fit_poisson_nmf_main_loop <- function (X, fit, numiter, method, control,
                                       verbose) {

  # Pre-compute quaantities and set up data structures used in the
  # loop below.
  loglik.const    <- loglik_poisson_const(X)
  dev.const       <- deviance_poisson_const(X)
  progress        <- as.data.frame(matrix(0,numiter,11))
  names(progress) <- c("loglik","dev","res","delta.f","delta.l","nonzeros.f",
                       "nonzeros.l","extrapolate","beta","betamax","timing")

  # Iterate the updates of the factors and loadings.
  for (i in 1:numiter) {
    fit0 <- fit
    t1   <- proc.time()

    # Update the factors and loadings.
    extrapolate <- with(control,extrapolate & (i %% extrapolate.reset != 0))
    if (extrapolate)
      fit <- update_poisson_nmf_extrapolated(X,fit,method,control)
    else
      fit <- update_poisson_nmf(X,fit,method,control)
    t2 <- proc.time()
    
    # Update the "progress" data frame with the log-likelihood,
    # deviance, and other quantities, and report the algorithm's
    # progress to the console if requested. In all cases, the "current
    # best" estimates of the factors and loadings are used to report
    # progress.
    u <- cost(X,fit$Lbest,t(fit$Fbest),control$eps,"poisson")
    progress[i,"loglik"]  <- sum(loglik.const - u)
    progress[i,"dev"]     <- sum(dev.const + 2*u)
    progress[i,"res"]     <- with(poisson_nmf_kkt(X,fit$Fbest,fit$Lbest),
                                  max(abs(rbind(F,L))))
    progress[i,"delta.f"] <- max(abs(fit$Fbest - fit0$Fbest))
    progress[i,"delta.l"] <- max(abs(fit$Lbest - fit0$Lbest))
    progress[i,"nonzeros.f"]  <- mean(fit$Fbest > control$zero.threshold.solution)
    progress[i,"nonzeros.l"]  <- mean(fit$Lbest > control$zero.threshold.solution)
    progress[i,"extrapolate"] <- extrapolate
    progress[i,"beta"]    <- fit$beta
    progress[i,"betamax"] <- fit$betamax
    progress[i,"timing"]  <- t2["elapsed"] - t1["elapsed"]
    if (verbose)
      cat(sprintf("%4d %+0.7e %+0.6e %0.2e %0.3e %0.3e %0.3f %0.3f %0.2f\n",
                  i,progress[i,"loglik"],progress[i,"dev"],progress[i,"res"],
                  progress[i,"delta.f"],progress[i,"delta.l"],
                  progress[i,"nonzeros.f"],progress[i,"nonzeros.l"],
                  extrapolate * progress[i,"beta"]))
  }

  # Output the updated "fit".
  fit$progress <- progress
  return(fit)
}

# This implements a single (non-extrapolated) update of the factors
# and loadings. The output is the updated "fit". Note that, because no
# extrapolation scheme is used, the extrapolated and non-extrapolated
# estimates are the same in the return value; that is, fit$Fy and
# fit$Fbest are the same as fit$F, and fit$Ly and fit$Lbest are the
# same as fit$L. The value of the loss function ("loss", "lossbest")
# is also updated.
update_poisson_nmf <- function (X, fit, method, control) {

  # Update the loadings ("activations"). The factors are forced to be
  # non-negative, or positive; the latter can promote better
  # convergence for some algorithms.
  fit$L <- update_loadings_poisson_nmf(X,fit$F,fit$L,method,control)
  fit$L <- pmax(fit$L,control$minval)

  # Update the factors ("basis vectors"). The loadings are forced to
  # be non-negative, or positive; the latter can promote better
  # convergence for some algorithms.
  fit$F <- update_factors_poisson_nmf(X,fit$F,fit$L,method,control)

  # Re-scale the factors and loadings.
  out   <- rescale.factors(fit$F,fit$L)
  fit$F <- out$F
  fit$L <- out$L

  # The extrapolated and "current best" estimates are the same as the
  # non-extrapolated estimates.
  fit$Fy    <- fit$F
  fit$Ly    <- fit$L
  fit$Fbest <- fit$F
  fit$Lbest <- fit$L

  # Compute the value of the objective ("loss") function at the updated
  # estimates.
  fit$loss     <- sum(cost(X,fit$L,t(fit$F),control$eps))
  fit$lossbest <- fit$loss
  
  # Output the updated "fit".
  return(fit)
}

# TO DO: Give a more detailed explanation of the inputs and outputs;
# in particular, explain what changes in "fit", and what does not.
#
# This implements an extrapolated update of the factors and loadings.
update_poisson_nmf_extrapolated <- function (X, fit, method, control) {

  # Store the value of the objective (loss) function at the current
  # iterate.
  loss0 <- fit$loss

  # Compute the extrapolated update for the loadings ("activations").
  # Note that when beta = 0, Ly = Ln.
  Ln     <- update_loadings_poisson_nmf(X,fit$Fy,fit$Ly,method,control)
  Ln     <- pmax(Ln,control$minval)
  fit$Ly <- pmax(Ln + fit$beta*(Ln - fit$L),control$minval)

  # Compute the extrapolated update for the factors ("basis vectors").
  # Note that when beta = 0, Fy = Fn.
  Fn     <- update_factors_poisson_nmf(X,fit$Fy,fit$Ly,method,control)
  Fn     <- pmax(Fn,control$minval)
  fit$Fy <- pmax(Fn + fit$beta*(Fn - fit$F),control$minval)
  
  # Compute the value of the objective (loss) function at the
  # extrapolated solution for the loadings (Ly) and the
  # non-extrapolated solution for the factors (Fn).
  fit$loss <- sum(cost(X,fit$Ly,t(Fn),control$eps))
  if (fit$beta == 0) {

    # When beta = 0, extrapolation is not used, and the extrapolation
    # parameters are not updated; use the basic coordinate-wise
    # updates for the factors and loadings.
    fit$F <- Fn
    fit$L <- Ln
  } else {

    # Update the extrapolation parameters following Algorithm 3 of
    # Ang & Gillis (2019).
    if (fit$loss > loss0) {

      # The solution did not improve, so restart the extrapolation.
      fit$Fy      <- fit$F
      fit$Ly      <- fit$L
      fit$betamax <- fit$beta0
      fit$beta    <- control$beta.reduce * fit$beta
    } else {
        
      # The solution improved; retain the basic co-ordinate ascent
      # update as well.
      fit$F       <- Fn
      fit$L       <- Ln
      fit$beta    <- min(fit$betamax,control$beta.increase * fit$beta)
      fit$beta0   <- fit$beta
      fit$betamax <- min(0.99,control$betamax.increase * fit$betamax)
    }
  }
  
  # If the solution is improved, update the current best estimates
  # using the non-extrapolated estimates of the factors (Fn) and the
  # extrapolated estimates of the loadings (Ly).
  if (fit$loss < fit$lossbest) {
    fit$Fbest    <- Fn
    fit$Lbest    <- fit$Ly 
    fit$lossbest <- fit$loss

    # Re-scale the factors and loadings.
    out       <- rescale.factors(fit$Fbest,fit$Lbest)
    fit$Fbest <- out$F
    fit$Lbest <- out$L
  }
           
  # Output the updated "fit".
  return(fit)
}

# Implements a single update of the factors matrix.
update_factors_poisson_nmf <- function (X, F, L, method, control) {
  if (method == "mu")
    F <- t(betanmf_update_factors(X,L,t(F)))
  else if (method == "em")
    F <- pnmfem_update_factors(X,F,L,control$numiter,control$nc)
  else if (method == "ccd")
    F <- t(ccd_update_factors(X,L,t(F),control$nc,control$eps))
  else if (method == "scd")
    F <- t(scd_update_factors(X,L,t(F),control$numiter,control$nc,control$eps))
  else if (method == "altsqp")
    F <- altsqp_update_factors(X,F,L,control$numiter,control)
  return(F)
}
  
# Implements a single update of the loadings matrix.
update_loadings_poisson_nmf <- function (X, F, L, method, control) {
  if (method == "mu")
    L <- betanmf_update_loadings(X,L,t(F))
  else if (method == "em")
    L <- pnmfem_update_loadings(X,F,L,control$numiter,control$nc)
  else if (method == "ccd")
    L <- ccd_update_loadings(X,L,t(F),control$nc,control$eps)
  else if (method == "scd")
    L <- scd_update_loadings(X,L,t(F),control$numiter,control$nc,control$eps)
  else if (method == "altsqp")
    L <- altsqp_update_loadings(X,F,L,control$numiter,control)
  return(L)
}

#' @rdname fit_poisson_nmf
#'
#' @export
#' 
fit_poisson_nmf_control_default <- function()
  c(list(numiter           = 1,
         minval            = 1e-15,
         nc                = as.integer(NA),
         extrapolate       = FALSE,
         extrapolate.reset = 20,
         beta.increase     = 1.1,
         beta.reduce       = 0.75,
         betamax.increase  = 1.05),
    mixsqp_control_default())
    
