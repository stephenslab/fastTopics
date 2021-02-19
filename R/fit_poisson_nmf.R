#' @title Fit Non-negative Matrix Factorization to Count Data
#'
#' @description Approximate the input matrix \code{X} by the
#'   non-negative matrix factorization \code{tcrossprod(L,F)}, in which
#'   the quality of the approximation is measured by a
#'   \dQuote{divergence} criterion; equivalently, optimize the
#'   likelihood under a Poisson model of the count data, \code{X}, in
#'   which the Poisson rates are given by \code{tcrossprod(L,F)}.
#'   Function \code{fit_poisson_nmf} runs a specified number of
#'   coordinate-wise updates to fit the L and F matrices.
#'
#' @details In Poisson non-negative matrix factorization (Lee & Seung,
#' 2001), counts \eqn{x_{ij}} in the \eqn{n \times m} matrix, \eqn{X},
#' are modeled by the Poisson distribution: \deqn{x_{ij} \sim
#' \mathrm{Poisson}(\lambda_{ij}).} Each Poisson rate,
#' \eqn{\lambda_{ij}}, is a linear combination of parameters
#' \eqn{f_{jk} \geq 0, l_{ik} \geq 0} to be fitted to the data:
#' \deqn{\lambda_{ij} = \sum_{k=1}^K l_{ik} f_{jk},} in which \eqn{K}
#' is a user-specified tuning parameter specifying the rank of the
#' matrix factorization. Function \code{fit_poisson_nmf} computes
#' maximum-likelihood estimates (MLEs) of the parameters.  For
#' additional mathematical background, and an explanation of how
#' Poisson NMF is connected to topic modeling, see the vignette:
#' \code{vignette(topic = "relationship",package = "fastTopics")}.
#'
#' Using this function requires some care; only minimal argument
#' checking is performed, and error messages may not be helpful.
#'
#' The EM and multiplicative updates are simple and fast, but can be
#' slow to converge to a stationary point. When \code{control$numiter
#' = 1}, the EM and multiplicative updates are mathematically
#' equivalent to the multiplicative updates, and therefore share the
#' same convergence properties. However, the implementation of the EM
#' updates is quite different; in particular, the EM updates are more
#' suitable for sparse counts matrices. The implementation of the
#' multiplicative updates is adapted from the MATLAB code by Daichi
#' Kitamura \url{http://d-kitamura.net}.
#'
#' Since the multiplicative updates are implemented using standard
#' matrix operations, the speed is heavily dependent on the
#' BLAS/LAPACK numerical libraries used. In particular, using
#' optimized implementations such as OpenBLAS or Intel MKL can result
#' in much improved performance of the multiplcative updates.
#'
#' The cyclic co-ordinate descent (CCD) and sequential co-ordinate
#' descent (SCD) updates adopt the same optimization strategy, but
#' differ in the implementation details. In practice, we have found
#' that the CCD and SCD updates arrive at the same solution when
#' initialized "sufficiently close" to a stationary point. The CCD
#' implementation is adapted from the C++ code developed by Cho-Jui
#' Hsieh and Inderjit Dhillon, which is available for download at
#' \url{https://www.cs.utexas.edu/~cjhsieh/nmf}. The SCD
#' implementation is based on version 0.4-3 of the NNLM package.
#'
#' An additional re-scaling step is performed after each update to
#' promote numerical stability.
#'
#' We use three measures of progress for the model fitting: (1)
#' improvement in the log-likelihood (or deviance), (2) change in the
#' model parameters, and (3) the residuals of the Karush-Kuhn-Tucker
#' (KKT) first-order conditions. As the iterates approach a stationary
#' point of the loss function, the change in the model parameters
#' should be small, and the residuals of the KKT system should vanish.
#' Use \code{\link{plot_progress_poisson_nmf}} to plot the improvement
#' in the solution over time.
#'
#' See \code{\link{fit_topic_model}} for additional guidance on model
#' fitting, particularly for large or complex data sets.
#' 
#' The \code{control} argument is a list in which any of the
#' following named components will override the default optimization
#' algorithm settings (as they are defined by
#' \code{fit_poisson_nmf_control_default}):
#' 
#' \describe{
#'
#' \item{\code{numiter}}{Number of "inner loop" iterations to run when
#'   performing and update of the factors or loadings. This must be set
#'   to 1 for \code{method = "mu"} and \code{method = "ccd"}.}
#'
#' \item{\code{nc}}{Number of RcppParallel threads to use for the
#'   updates. When \code{nc} is \code{NA}, the default number of threads
#'   is used; see \code{\link[RcppParallel]{defaultNumThreads}}. This
#'   setting is ignored for the multiplicative upates (\code{method =
#'   "mu"}).}
#'
#' \item{\code{minval}}{A small, positive constant used to safeguard
#'   the multiplicative updates. The safeguarded updates are implemented
#'   as \code{F <- pmax(F1,minval)} and \code{L <- pmax(L1,minval)},
#'   where \code{F1} and \code{L1} are the factors and loadings matrices
#'   obtained by applying an update. This is motivated by Theorem 1 of
#'   Gillis & Glineur (2012). Setting \code{minval = 0} is allowed, but
#'   some methods are not guaranteed to converge to a stationary point
#'   without this safeguard, and a warning will be given in this case.}
#'
#' \item{\code{extrapolate}}{When \code{extrapolate = TRUE}, the
#'   extrapolation scheme of Ang & Gillis (2019) is used.}
#'
#' \item{\code{extrapolate.reset}}{To promote better numerical
#'   stability of the extrapolated updates, they are \dQuote{reset}
#'   every so often. This parameter determines the number of iterations
#'   to wait before resetting.}
#'
#' \item{\code{beta.increase}}{When the extrapolated update improves
#'   the solution, scale the extrapolation parameter by this amount.}
#'
#' \item{\code{beta.reduce}}{When the extrapolaaed update does not
#'   improve the solution, scale the extrapolation parameter by this
#'   amount.}
#'
#' \item{\code{betamax.increase}}{When the extrapolated update
#'   improves the solution, scale the extrapolation parameter by this
#'   amount.}
#'
#' \item{\code{eps}}{A small, non-negative number that is added to the
#'   terms inside the logarithms to sidestep computing logarithms of
#'   zero. This prevents numerical problems at the cost of introducing a
#'   small inaccuracy in the solution. Increasing this number may lead
#'   to faster convergence but possibly a less accurate solution.}
#'
#' \item{\code{zero.threshold}}{A small, non-negative number used to
#'   determine which entries of the solution are exactly zero. Any
#'   entries that are less than or equal to \code{zero.threshold} are
#'   considered to be exactly zero.}}
#'
#' An additional setting, \code{control$init.numiter}, controls the
#' number of sequential co-ordinate descent (SCD) updates that are
#' performed to initialize the loadings matrix when \code{init.method
#' = "topicscore"}.
#'
#' @param X The n x m matrix of counts; all entries of X should be
#'   non-negative. It can be a sparse matrix (class \code{"dgCMatrix"})
#'   or dense matrix (class \code{"matrix"}), with some exceptions (see
#'   \sQuote{Details}).
#'
#' @param k An integer 2 or greater giving the matrix rank. This
#'   argument should only be specified if the initial fit (\code{fit0}
#'   or \code{F, L}) is not provided.
#' 
#' @param fit0 The initial model fit. It should be an object of class
#'   \dQuote{poisson_nmf_fit}, such as an output from
#'   \code{init_poisson_nmf}, or from a previous call to
#'   \code{fit_poisson_nmf}.
#'
#' @param numiter The number of updates of the factors and loadings to
#'   perform.
#'
#' @param update.factors A numeric vector specifying which factors
#'   (rows of \code{F}) to update. By default, all factors are
#'   updated. Note that the rows that are not updated may still change
#'   by rescaling. When \code{NULL}, all factors are fixed. This option
#'   is only implemented for \code{method = "em"} and \code{method =
#'   "scd"}. If another method is selected, the default setting of
#'   \code{update.factors} must be used.
#'
#' @param update.loadings A numeric vector specifying which loadings
#'   (rows of \code{L}) to update. By default, all loadings are
#'   updated. Note that the rows that are not updated may still change
#'   by rescaling. When \code{NULL}, all loadings are fixed. This option
#'   is only implemented for \code{method = "em"} and \code{method =
#'   "scd"}. If another method is selected, the default setting of
#'   \code{update.loadings} must be used.
#' 
#' @param method The method to use for updating the factors and
#'   loadings. Four methods are implemented: multiplicative updates,
#'   \code{method = "mu"}; expectation maximization (EM), \code{method =
#'   "em"}; sequential co-ordinate descent (SCD), \code{method = "scd"};
#'   and cyclic co-ordinate descent (CCD), \code{method = "ccd"}. See
#'   \sQuote{Details} for a detailed description of these methods.
#'
#' @param init.method The method used to initialize the factors and
#'   loadings. When \code{init.method = "random"}, the factors and
#'   loadings are initialized uniformly at random; when
#'   \code{init.method = "topicscore"}, the factors are initialized
#'   using the (very fast) Topic SCORE algorithm (Ke & Wang, 2017), and
#'   the loadings are initialized by running a small number of SCD
#'   updates. This input argument is ignored if initial estimates of the
#'   factors and loadings are already provided via input \code{fit0}, or
#'   inputs \code{F} and \code{L}.
#' 
#' @param control A list of parameters controlling the behaviour of
#'   the optimization algorithm (and the Topic SCORE algorithm if it
#'   is used to initialize the model parameters). See \sQuote{Details}.
#' 
#' @param verbose When \code{verbose = TRUE}, information about the
#'   algorithm's progress is printed to the console at each
#'   iteration. For interpretation of the columns, see the description
#'   of the \code{progress} return value.
#'
#' @return \code{init_poisson_nmf} and \code{fit_poisson_nmf} both
#' return an object capturing the optimization algorithm state (for
#' \code{init_poisson_nmf}, this is the initial state). It is a list
#' with the following elements:
#'
#' \item{F}{A matrix containing the current best estimates of the
#'   factors.}
#' 
#' \item{L}{A matrix containing the current best estimates of the
#'   loadings.}
#' 
#' \item{Fn}{A matrix containing the non-extrapolated factor estimates.
#'   If extrapolation is not used, \code{Fn} and \code{F} will be the
#'   same.}
#'
#' \item{Ln}{A matrix containing the non-extrapolated estimates of the
#'   loadings. If extrapolation is not used, \code{Ln} and \code{L} will
#'   be the same.}
#'
#' \item{Fy}{A matrix containing the extrapolated factor estimates. If
#'   the extrapolation scheme is not used, \code{Fy} and \code{F} will
#'   be the same.}
#'
#' \item{Ly}{A matrix containing the extrapolated estimates of the
#'   loadings. If extrapolation is not used, \code{Ly} and \code{L} will
#'   be the same.}
#'
#' \item{loss}{Value of the objective (\dQuote{loss}) function
#'   computed at the current best estimates of the factors and
#'   loadings.}
#' 
#' \item{loss.fnly}{Value of the objective (\dQuote{loss}) function
#'   computed at the extrapolated solution for the loadings (\code{Ly})
#'   and the non-extrapolated solution for the factors (\code{Fn}). This
#'   is used internally to implement the extrapolated updates.}
#'
#' \item{iter}{The number of the most recently completed iteration.}
#' 
#' \item{beta}{The extrapolation parameter, \eqn{beta} in Algorithm 3
#'   of Ang & Gillis (2019).}
#'
#' \item{betamax}{Upper bound on the extrapolation parameter. This is
#'   \eqn{\bar{\gamma}} in Algorithm 3 of Ang & Gillis (2019).}
#'
#' \item{beta0}{The setting of the extrapolation parameter at the
#'   last iteration that improved the solution.}
#' 
#' \item{progress}{A data frame containing detailed information about
#'   the algorithm's progress. The data frame should have \code{numiter}
#'   rows. The columns of the data frame are: "iter", the iteration
#'   number; "loglik", the log-likelihood at the current best factor and
#'   loading estimates; "dev", the deviance at the current best factor
#'   and loading estimates; "res", the maximum residual of the
#'   Karush-Kuhn-Tucker (KKT) first-order optimality conditions at the
#'   current best factor and loading estimates; "delta.f", the largest
#'   change in the factors matrix; "delta.l", the largest change in the
#'   loadings matrix; "nonzeros.f", the proportion of entries in the
#'   factors matrix that are nonzero; "nonzeros.l", the proportion of
#'   entries in the loadings matrix that are nonzero; "extrapolate",
#'   which is 1 if extrapolation is used, otherwise it is 0; "beta", the
#'   setting of the extrapolation parameter; "betamax", the setting of
#'   the extrapolation parameter upper bound; and "timing", the elapsed
#'   time in seconds (recorded using \code{\link{proc.time}}).}
#' 
#' @references
#'
#'   Ang, A. and Gillis, N. (2019). Accelerating nonnegative matrix
#'   factorization algorithms using extrapolation. \emph{Neural
#'   Computation} \bold{31}, 417–439.
#' 
#'   Cichocki, A., Cruces, S. and Amari, S. (2011). Generalized
#'   alpha-beta divergences and their application to robust nonnegative
#'   matrix factorization. \emph{Entropy} \bold{13}, 134–170.
#' 
#'   Gillis, N. and Glineur, F. (2012). Accelerated multiplicative
#'   updates and hierarchical ALS algorithms for nonnegative matrix
#'   factorization. \emph{Neural Computation} \code{24}, 1085–1105.
#'
#'   Hsieh, C.-J. and Dhillon, I. (2011). Fast coordinate descent
#'   methods with variable selection for non-negative matrix
#'   factorization. In \emph{Proceedings of the 17th ACM SIGKDD
#'   international conference on Knowledge discovery and data mining},
#'   p. 1064-1072
#'
#'   Lee, D. D. and Seung, H. S. (2001). Algorithms for non-negative
#'   matrix factorization. In \emph{Advances in Neural Information
#'   Processing Systems} \bold{13}, 556–562.
#'
#'   Lin, X. and Boutros, P. C. (2018). Optimization and expansion of
#'   non-negative matrix factorization. \emph{BMC Bioinformatics}
#'   \bold{21}, 7.
#'
#'   Ke, Z. & Wang, M. (2017). A new SVD approach to optimal topic
#'   estimation. \emph{arXiv} \url{http://arxiv.org/abs/1704.07016}
#'
#' @seealso \code{\link{fit_topic_model}},
#'   \code{\link{plot_progress_poisson_nmf}}
#' 
#' @examples
#' # Simulate a (sparse) 80 x 100 counts matrix.
#' library(Matrix)
#' set.seed(1)
#' X <- simulate_count_data(80,100,k = 3,sparse = TRUE)$X
#' 
#' # Remove columns (words) that do not appear in any row (document).
#' X <- X[,colSums(X > 0) > 0]
#' 
#' # Run 10 EM updates to find a good initialization.
#' fit0 <- fit_poisson_nmf(X,k = 3,numiter = 10)
#' 
#' # Fit the Poisson NMF model by running 50 EM updates.
#' fit_em <- fit_poisson_nmf(X,fit0 = fit0,numiter = 50,method = "em")
#' 
#' # Fit the Poisson NMF model by running 50 extrapolated SCD updates.
#' fit_scd <- fit_poisson_nmf(X,fit0 = fit0,numiter = 50,method = "scd",
#'                            control = list(extrapolate = TRUE))
#' 
#' # Compare the two fits.
#' fits <- list(em = fit_em,scd = fit_scd)
#' compare_poisson_nmf_fits(fits)
#' plot_progress_poisson_nmf(fits,y = "loglik")
#' plot_progress_poisson_nmf(fits,y = "res")
#' 
#' # Recover the topic model. After this step, the L matrix contains the
#' # mixture proportions ("loadings"), and the F matrix contains the
#' # word frequencies ("factors").
#' fit_multinom <- poisson2multinom(fit_scd)
#'
#' @useDynLib fastTopics
#'
#' @importFrom utils modifyList
#'
#' @export
#'
fit_poisson_nmf <- function (X, k, fit0, numiter = 100,
                             update.factors = seq(1,ncol(X)),
                             update.loadings = seq(1,nrow(X)),
                             method = c("scd","em","mu","ccd"),
                             init.method = c("topicscore","random"),
                             control = list(), verbose = TRUE) {

  # CHECK & PROCESS INPUTS
  # ----------------------
  # Check and process input argument "X".
  verify.count.matrix(X)
  if (is.matrix(X) & length(X) > 1e4 & mean(X > 0) < 0.1 & any(method != "mu"))
    message(paste("Input matrix \"X\" has less than 10% nonzero entries;",
                  "consider converting \"X\" to a sparse matrix to reduce",
                  "computational effort"))
  if (is.matrix(X) & is.integer(X))
    storage.mode(X) <- "double"
  
  # Get the number of rows (n) and columns (m) of the data matrix,
  n <- nrow(X)
  m <- ncol(X)

  # Check input argument "numiter".
  if (any(numiter < 1))
    stop("Input argument \"numiter\" must be 1 or greater")

  # Process input arguments "update.factors" and "update.loadings".
  update.factors  <- sort(update.factors)
  update.loadings <- sort(update.loadings)
  if (length(update.factors) == 0 & length(update.loadings) == 0)
    stop("None of the factors or loadings have been selected for updating")
  
  # Check and process input arguments "method" and "init.method".
  method      <- match.arg(method)
  init.method <- match.arg(init.method)
  if (method == "mu" & is.sparse.matrix(X)) {
    warning("method = \"mu\" is not implemented for sparse counts matrices; ",
            "attempting to convert \"X\" to a dense matrix")
    X <- as.matrix(X)
  }
  if (!(length(update.factors) == m & length(update.loadings) == n))
    if (method == "mu" | method == "ccd")
      stop("All factors and loadings must be updated for method = \"mu\" ",
           "and method = \"ccd\"")
  
  # Check and process the optimization settings.
  control <- modifyList(fit_poisson_nmf_control_default(),
                        control,keep.null = TRUE)
  if (any(control$minval == 0))
    warning("Updates may not converge when minval = 0")
  if (control$numiter > 1 & (method == "mu" | method == "ccd")) {
    warning("multiplicative and CCD updates do not allow ",
            "control$numiter > 1; setting control$numiter = 1")
    control$numiter <- 1
  }
  if ((length(update.factors) == 0 | length(update.loadings) == 0) &
      control$extrapolate)
    stop("control$extrapolate cannot be TRUE when all factors or loadings ",
         "are fixed")
  control$nc <- initialize.multithreading(control$nc)
  
  # Only one of "k" and "fit0" should be provided. If argument "k" is
  # given, generate a random initialization of the factors and
  # loadings.
  if (!(missing(k) & !missing(fit0) | (!missing(k) & missing(fit0))))
    stop("Provide a rank, \"k\", or an initialization, \"fit0\", but not both")
  if (missing(fit0))
    fit0 <- init_poisson_nmf(X,k = k,init.method = init.method,
                             control = control,verbose = verbose)

  # Check input argument "fit0".
  if (!inherits(fit0,"poisson_nmf_fit"))
    stop("Input argument \"fit0\" should be an object of class ",
         "\"poisson_nmf_fit\", such as an output of init_poisson_nmf")
  verify.fit.and.count.matrix(X,fit0)
  fit <- fit0
  k   <- ncol(fit$F)
  
  # Give an overview of the optimization settings.
  if (verbose) {
    cat(sprintf("Fitting rank-%d Poisson NMF to %d x %d %s matrix.\n",k,n,m,
                ifelse(is.matrix(X),"dense","sparse")))
    if (method == "mu")
      method.text <- "multiplicative"
    else if (method == "em")
      method.text <- "EM"
    else if (method == "scd")
      method.text <- "SCD"
    else if (method == "ccd")
      method.text <- "CCD"
    cat(sprintf("Running %d %s updates, %s extrapolation ",numiter,
        method.text,ifelse(control$extrapolate,"with","without")))
    cat("(fastTopics 0.5-2).\n")
  }
  
  # INITIALIZE ESTIMATES
  # --------------------
  # Re-scale the initial estimates of the factors and loadings so that
  # they are on same scale on average. This is intended to improve
  # numerical stability of the updates.
  fit <- rescale.fit(fit)

  # Initialize the estimates of the factors and loadings. To prevent
  # the updates from getting "stuck", force the initial estimates to
  # be positive.
  fit <- safeguard.fit(fit,control$minval)

  # RUN UPDATES
  # -----------
  if (verbose)
    cat("iter   log-likelihood        deviance res(KKT)  |F-F'|  |L-L'|",
        "nz(F) nz(L) beta\n")
  fit <- fit_poisson_nmf_main_loop(X,fit,numiter,update.factors,
                                   update.loadings,method,control,
                                   verbose)

  # Output the updated "fit".
  fit$progress     <- rbind(fit0$progress,fit$progress)
  dimnames(fit$F)  <- dimnames(fit0$F)
  dimnames(fit$L)  <- dimnames(fit0$L)
  dimnames(fit$Fn) <- dimnames(fit0$Fn)
  dimnames(fit$Ln) <- dimnames(fit0$Ln)
  dimnames(fit$Fy) <- dimnames(fit0$Fy)
  dimnames(fit$Ly) <- dimnames(fit0$Ly)
  class(fit) <- c("poisson_nmf_fit","list")
  return(fit)
}

#' @rdname fit_poisson_nmf
#'
#' @param F An optional argument giving is the initial estimate of the
#'   factors (also known as \dQuote{basis vectors}). It should be an m x
#'   k matrix, where m is the number of columns in the counts matrix
#'   \code{X}, and k > 1 is the rank of the matrix factorization
#'   (equivalently, the number of \dQuote{topics}). All entries of
#'   \code{F} should be non-negative. When \code{F} and \code{L} are not
#'   provided, input argument \code{k} should be specified instead.
#'
#' @param L An optional argument giving the initial estimate of the
#'   loadings (also known as \dQuote{activations}). It should an n x k
#'   matrix, where n is the number of rows in the counts matrix
#'   \code{X}, and k > 1 is the rank of the matrix factorization
#'   (equivalently, the number of \dQuote{topics}). All entries of
#'   \code{L} should be non-negative. When \code{F} and \code{L} are not
#'   provided, input argument \code{k} should be specified instead.
#'
#' @param beta Initial setting of the extrapolation parameter. This is
#'   \eqn{beta} in Algorithm 3 of Ang & Gillis (2019).
#'
#' @param betamax Initial setting for the upper bound on the
#'   extrapolation parameter. This is \eqn{\bar{\gamma}} in Algorithm 3
#'   of Ang & Gillis (2019).
#'
#' @importFrom Matrix rowSums
#' 
#' @export
#' 
init_poisson_nmf <-
  function (X, F, L, k, init.method = c("topicscore","random"),
            beta = 0.5, betamax = 0.99, control = list(),
            verbose = TRUE) {

  # Check input X.
  verify.count.matrix(X)

  # Get the number of rows (n) and columns (m) of the data matrix,
  n <- nrow(X)
  m <- ncol(X)

  # Only one of k and (F or L) should be provided.
  if (!(missing(k) & !(missing(F) & missing(L)) |
       (!missing(k) & (missing(F) & missing(L)))))
    stop("Provide a rank (k) or an initialization F and/or L, but not both")
  if (missing(k)) {
    if (missing(F))
      k <- ncol(L)
    else
      k <- ncol(F)
  }
  if (any(k < 2))
    stop("Matrix factorization rank \"k\" should be 2 or greater")

  # Check and process input argument "init.method".
  init.method <- match.arg(init.method)
  
  # Check and process the optimization settings.
  control <- modifyList(fit_poisson_nmf_control_default(),
                        control,keep.null = TRUE)

  # If the factor and loadings matrices are not provided, initialize
  # them using the chosen method. If the factor and loadings matrices
  # are provided, run a few simple checks to verify that they are
  # valid initial estimates.
  if (missing(F) | missing(L)) {

    if (init.method == "random") {

      # Initialize the factors and loadings uniformly at random.
      if (missing(F))
        F <- rand(m,k)
      if (missing(L))
        L <- rand(n,k)
    } else if (init.method == "topicscore") {
      if (!(missing(F) & missing(L)))
        stop("init.method = \"topicscore\" can only be used when both L ",
             "and F are not provided")
      
      # Initialize the factors using the "Topic SCORE" algorithm.
      if (verbose)
        cat("Initializing factors using Topic SCORE algorithm.\n")
      F <- tryCatch(topic_score(X,k),
                    error = function (e) {
                      warning("Topic SCORE failure occurred; falling back ",
                              "to init.method == \"random\"")
                      return(NULL)
                    })
      if (is.null(F)) {
        
        # The Topic SCORE algorithm failed; generate a random
        # initialization instead.
        if (verbose)
          cat("Topic SCORE failure occurred; using random initialization",
              "instead.\n")
        F <- rand(m,k)
        L <- rand(n,k)
    } else {
        
        # Fit the loadings by running a small number of sequential
        # co-ordinate descent (SCD) updates.
        if (verbose)
          cat(sprintf("Initializing loadings by running %d SCD updates.\n",
                      control$init.numiter))
        control$nc <- initialize.multithreading(control$nc)
        s <- rowSums(X)
        L <- matrix(s,n,k)
        L <- scd_update_loadings(X,L,t(F),1:n,control$init.numiter,
                                 control$nc,control$eps)
      }
    }
    rownames(F) <- colnames(X)
    rownames(L) <- rownames(X)
    colnames(F) <- paste0("k",1:k)
    colnames(L) <- paste0("k",1:k)
  } else {

    # Check the provided factor matrix.
    if (!(is.matrix(F) & is.numeric(F)))
      stop("Input argument \"F\" should be a numeric matrix (is.matrix(F) ",
           "should return TRUE)")
    if (is.integer(F))
      storage.mode(F) <- "double"

    # Check the provided loadings matrix.
    if (!(is.matrix(L) & is.numeric(L)))
      stop("Input argument \"L\" should be a numeric matrix (is.matrix(L) ",
           "should return TRUE)")
    if (is.integer(L))
      storage.mode(L) <- "double"
  }

  # Force the initial estimates to be positive.
  F <- pmax(F,control$minval)
  L <- pmax(L,control$minval)

  # Compute the value of the objective ("loss") function at the
  # initial estimates of the factors and loading.
  loss <- sum(cost(X,L,t(F),control$eps))

  # Initialize the data frame for keeping track of the algorithm's
  # progress over time.
  progress        <- as.data.frame(matrix(0,0,12))
  names(progress) <- c("iter","loglik","dev","res","delta.f","delta.l",
                       "nonzeros.f","nonzeros.l","extrapolate","beta",
                       "betamax","timing")

  # Return a list containing: F, an initial estimate of the factors;
  # L, an initial estimate of the loadings; Fn and Ln, the
  # non-extrapolated estimates of the factors and loadings, which
  # initially are always the same as F and L; Fy and Ly, the
  # extrapolated estimates of the factors and loadings, which
  # initially are always the same as F and L; "iter", the current
  # iteration number (initially set to zero); "loss" and "loss.fnly",
  # the value of the objective or loss function at the current best
  # and partially extrapolated estimates of the factors and loadings;
  # beta, beta0 and betamax, the initial settings of the extrapolation
  # parameters; and "progress", the initial data frame for keeping
  # track of the algorithm's progress over time.
  fit <- list(F         = F,
              L         = L,
              Fn        = F,
              Ln        = L,
              Fy        = F,
              Ly        = L,
              loss      = loss,
              loss.fnly = loss,
              iter      = 0,
              beta      = beta,
              beta0     = beta,
              betamax   = betamax,
              progress  = progress)
  class(fit) <- c("poisson_nmf_fit","list")
  return(fit)
}

# This implements the core part of fit_poisson_nmf.
fit_poisson_nmf_main_loop <- function (X, fit, numiter, update.factors,
                                       update.loadings, method, control,
                                       verbose) {
    
  # Pre-compute quantities and set up data structures used in the
  # loop below.
  loglik.const    <- sum(loglik_poisson_const(X))
  dev.const       <- sum(deviance_poisson_const(X))
  progress        <- as.data.frame(matrix(0,numiter,12))
  names(progress) <- c("iter","loglik","dev","res","delta.f","delta.l",
                       "nonzeros.f","nonzeros.l","extrapolate","beta",
                       "betamax","timing")

  # Iterate the updates of the factors and loadings.
  for (i in 1:numiter) {
    fit0 <- fit
    t1   <- proc.time()

    # Update the factors and loadings.
    if (control$extrapolate &
        fit$beta > 0 &
        i %% control$extrapolate.reset != 0) {
        
      # Perform an "extrapolated" update of the factors and loadings.
      extrapolate <- TRUE
      fit <- update_poisson_nmf_extrapolated(X,fit,update.factors,
                                             update.loadings,method,
                                             control)
    } else {

      # Perform a basic coordinate-wise update of the factors and
      # loadings.
      extrapolate <- FALSE
      fit <- update_poisson_nmf(X,fit,update.factors,update.loadings,
                                method,control)
    }
    t2 <- proc.time()

    # Update the iteration number.
    fit$iter <- fit$iter + 1
    
    # Update the "progress" data frame with the log-likelihood,
    # deviance, and other quantities, and report the algorithm's
    # progress to the console if requested. In all cases, the "current
    # best" estimates of the factors and loadings are used to report
    # progress.
    progress[i,"iter"]        <- fit$iter
    progress[i,"loglik"]      <- loglik.const - fit$loss
    progress[i,"dev"]         <- dev.const + 2*fit$loss
    progress[i,"res"]         <- with(poisson_nmf_kkt(X,fit$F,fit$L),
                                      max(abs(rbind(F[update.factors,],
                                                    L[update.loadings,]))))
    progress[i,"delta.f"]     <- max(abs(fit$F - fit0$F))
    progress[i,"delta.l"]     <- max(abs(fit$L - fit0$L))
    progress[i,"beta"]        <- fit$beta
    progress[i,"betamax"]     <- fit$betamax
    progress[i,"timing"]      <- t2["elapsed"] - t1["elapsed"]
    progress[i,"nonzeros.f"]  <- mean(fit$F > control$zero.threshold)
    progress[i,"nonzeros.l"]  <- mean(fit$L > control$zero.threshold)
    progress[i,"extrapolate"] <- extrapolate
    if (verbose)
      cat(sprintf("%4d %+0.9e %+0.8e %0.2e %0.1e %0.1e %0.3f %0.3f %0.2f\n",
                  fit$iter,progress[i,"loglik"],progress[i,"dev"],
                  progress[i,"res"],progress[i,"delta.f"],
                  progress[i,"delta.l"],progress[i,"nonzeros.f"],
                  progress[i,"nonzeros.l"],extrapolate * progress[i,"beta"]))
  }

  # Output the updated "fit".
  fit$progress <- progress
  return(fit)
}

# This implements a single (non-extrapolated) update of the factors
# and loadings.
#
# The output is the updated "fit". Because no extrapolation scheme is
# used, the extrapolated and non-extrapolated estimates are the same
# in the return value; that is, Fy and Fn are the same as F, and Ly
# and Ln are the same as L. The value of the loss function ("loss",
# "loss.fnly") is also updated.
#
# Note that "update.factors" and "update.loadings" is ignored for
# method = "mu" and method = "ccd".
update_poisson_nmf <- function (X, fit, update.factors, update.loadings,
                                method, control) {

  # Update the loadings ("activations"). The factors are forced to be
  # non-negative, or positive; the latter can promote better
  # convergence for some algorithms.
  if (length(update.loadings) > 0) {
    fit$L <- update_loadings_poisson_nmf(X,fit$F,fit$L,update.loadings,
                                         method,control)
    fit$L <- pmax(fit$L,control$minval)
  }

  # Update the factors ("basis vectors"). The loadings are forced to
  # be non-negative, or positive; the latter can promote better
  # convergence for some algorithms.
  if (length(update.factors) > 0) {
    fit$F <- update_factors_poisson_nmf(X,fit$F,fit$L,update.factors,
                                        method,control)
    fit$F <- pmax(fit$F,control$minval)
  }
  
  # The extrapolated and "current best" estimates are the same as the
  # non-extrapolated estimates.
  fit$Fy <- fit$F
  fit$Ly <- fit$L
  fit$Fn <- fit$F
  fit$Ln <- fit$L

  # Re-scale the factors and loadings.
  fit <- rescale.fit(fit)

  # Compute the value of the objective ("loss") function at the updated
  # estimates.
  fit$loss      <- sum(cost(X,fit$L,t(fit$F),control$eps))
  fit$loss.fnly <- fit$loss
  
  # Output the updated "fit".
  return(fit)
}

# This implements an extrapolated update of the factors and loadings.
#
# The output is the updated "fit". The updated parts of the "fit" are:
#
#   F          "best current" factor estimates
#   Fn         non-extrapolated factor estimates
#   Fy         extrapolated factor estimates
#   L          "best current" loadings estimates
#   Ln         non-extrapolated loadings estimates
#   Ly         extrapolated loadings estimates
#   loss       value of the objective at the "best current" estimates 
#   loss.fnly  value of the objective at (Fn, Ly)
#   beta       extrapolation parameter
#   betamax    upper bound on the extrapolation parameter
#   beta0      extrapolation parameter setting from last improvement
#
# Note that "update.factors" and "update.loadings" is ignored for
# method = "mu" and method = "ccd".
update_poisson_nmf_extrapolated <- function (X, fit, update.factors,
                                             update.loadings, method,
                                             control) {

  # Store the value of the objective (loss) function at the current
  # iterate (Fn, Ly).
  loss0.fnly <- fit$loss.fnly

  # Compute the extrapolated update for the loadings ("activations").
  # Note that when beta is zero, Ly = Ln.
  Ln     <- update_loadings_poisson_nmf(X,fit$Fy,fit$Ly,update.loadings,
                                        method,control)
  Ln     <- pmax(Ln,control$minval)
  fit$Ly <- pmax(Ln + fit$beta*(Ln - fit$Ln),control$minval)

  # Compute the extrapolated update for the factors ("basis vectors").
  # Note that when beta = 0, Fy = Fn.
  Fn     <- update_factors_poisson_nmf(X,fit$Fy,fit$Ly,update.factors,
                                       method,control)
  Fn     <- pmax(Fn,control$minval)
  fit$Fy <- pmax(Fn + fit$beta*(Fn - fit$Fn),control$minval)
  
  # Compute the value of the objective (loss) function at the
  # extrapolated solution for the loadings (Ly) and the
  # non-extrapolated solution for the factors (Fn).
  fit$loss.fnly <- sum(cost(X,fit$Ly,t(Fn),control$eps))

  # Update the extrapolation parameters following Algorithm 3 of
  # Ang & Gillis (2019).
  if (fit$loss.fnly >= loss0.fnly) {

    # The solution did not improve, so restart the extrapolation
    # scheme.
    fit$Fy      <- fit$Fn
    fit$Ly      <- fit$Ln
    fit$betamax <- fit$beta0
    fit$beta    <- control$beta.reduce * fit$beta
  } else {
        
    # The solution improved; retain the basic co-ordinate ascent
    # update as well.
    fit$Fn      <- Fn
    fit$Ln      <- Ln
    fit$beta    <- min(fit$betamax,control$beta.increase * fit$beta)
    fit$beta0   <- fit$beta
    fit$betamax <- min(0.99,control$betamax.increase * fit$betamax)
  }
 
  # If the solution improves the "current best" estimate, update the
  # current best estimate using the non-extrapolated estimates of the
  # factors (Fn) and the extrapolated estimates of the loadings (Ly).
  if (fit$loss.fnly < fit$loss) {
    fit$F    <- Fn
    fit$L    <- fit$Ly 
    fit$loss <- fit$loss.fnly
  }

  # Re-scale the factors and loadings.
  fit <- rescale.fit(fit)
  
  # Output the updated "fit".
  return(fit)
}

# Implements a single update of the factors matrix. Note that input
# argument "j", specifying which columns of F to update, is ignored for
# method = "mu" and method = "ccd".
update_factors_poisson_nmf <- function (X, F, L, j, method, control) {
  numiter <- control$numiter
  nc      <- control$nc
  eps     <- control$eps
  if (method == "mu")
    F <- t(betanmf_update_factors(X,L,t(F)))
  else if (method == "em")
    F <- pnmfem_update_factors(X,F,L,j,numiter,nc)
  else if (method == "ccd")
    F <- t(ccd_update_factors(X,L,t(F),nc,eps))
  else if (method == "scd")
    F <- t(scd_update_factors(X,L,t(F),j,numiter,nc,eps))
  return(F)
}
  
# Implements a single update of the loadings matrix. Note that input
# argument "i", specifying which rows of L to update, is ignored for
# method = "mu" and method = "ccd".
update_loadings_poisson_nmf <- function (X, F, L, i, method, control) {
  numiter <- control$numiter
  nc      <- control$nc
  eps     <- control$eps
  if (method == "mu")
    L <- betanmf_update_loadings(X,L,t(F))
  else if (method == "em")
    L <- pnmfem_update_loadings(X,F,L,i,numiter,nc)
  else if (method == "ccd")
    L <- ccd_update_loadings(X,L,t(F),nc,eps)
  else if (method == "scd")
    L <- scd_update_loadings(X,L,t(F),i,numiter,nc,eps)
  return(L)
}

# Rescale the extrapolated, non-extrapolated, and "current best"
# factors and loadings.
rescale.fit <- function (fit) {

  # Rescale the "current best" factors and loadings.
  out    <- rescale.factors(fit$F,fit$L)
  fit$F  <- out$F
  fit$L  <- out$L

  # Rescale the non-extrapolated factors and loadings.
  out    <- rescale.factors(fit$Fn,fit$Ln)
  fit$Fn <- out$F
  fit$Ln <- out$L

  # Rescale the extrapolated factors and loadings.
  out    <- rescale.factors(fit$Fy,fit$Ly)
  fit$Fy <- out$F
  fit$Ly <- out$L

  return(fit)
}

# Force the factors and loadings to be positive.
safeguard.fit <- function (fit, minval) {
  fit$F  <- pmax(fit$F,minval)
  fit$L  <- pmax(fit$L,minval)
  fit$Fn <- pmax(fit$Fn,minval)
  fit$Ln <- pmax(fit$Ln,minval)
  fit$Fy <- pmax(fit$Fy,minval)
  fit$Ly <- pmax(fit$Ly,minval)
  return(fit)
}

#' @rdname fit_poisson_nmf
#'
#' @export
#' 
fit_poisson_nmf_control_default <- function()
  list(numiter           = 4,
       init.numiter      = 10,
       minval            = 1e-15,
       eps               = 1e-8,
       zero.threshold    = 1e-6,
       nc                = as.integer(NA),
       extrapolate       = FALSE,
       extrapolate.reset = 20,
       beta.increase     = 1.1,
       beta.reduce       = 0.75,
       betamax.increase  = 1.05)
