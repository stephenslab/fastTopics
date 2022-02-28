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
#'   loadings (also known as \dQuote{activations}). It should be an n x k
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
            verbose = c("detailed","none")) {

  # Check input X.
  verify.count.matrix(X)

  # Get the number of rows (n) and columns (m) of the data matrix,
  n <- nrow(X)
  m <- ncol(X)

  # Only one of k and (F or L) should be provided.
  if (!(missing(k) & !(missing(F) & missing(L)) |
       (!missing(k) & (missing(F) & missing(L)))))
    stop("Provide a rank (k) or an initialization of F and/or L, but not both")
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

  # Check and progress input argument "verbose".
  verbose <- match.arg(verbose)
  
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
      if (verbose == "detailed")
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
        if (verbose == "detailed")
          cat("Topic SCORE failure occurred; using random initialization",
              "instead.\n")
        F <- rand(m,k)
        L <- rand(n,k)
      } else {
        
        # Fit the loadings by running a small number of sequential
        # co-ordinate descent (SCD) updates.
        if (verbose == "detailed")
          cat(sprintf("Initializing loadings by running %d SCD updates.\n",
                      control$init.numiter))
        control$nc <- initialize.multithreading(control$nc,verbose != "none")
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
  progress        <- as.data.frame(matrix(0,0,13))
  names(progress) <- c("iter","loglik","loglik.multinom","dev","res",
                       "delta.f","delta.l","nonzeros.f","nonzeros.l",
                       "extrapolate","beta","betamax","timing")

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

#' @rdname fit_poisson_nmf
#'
#' @param clusters A factor specifying a grouping, or clustering, of
#'   the rows of \code{X}.
#'
#' @param \dots Additional arguments passed to \code{init_poisson_nmf}.
#'
#' @importFrom Matrix rowSums
#' @importFrom Matrix colSums
#'
#' @export
#'
init_poisson_nmf_from_clustering <- function (X, clusters, ...) {

  # Check input X.
  if (!((is.numeric(X) & is.matrix(X)) | is.sparse.matrix(X)))
    stop("Input argument \"X\" should be a numeric matrix (a \"matrix\" or ",
         "a \"dgCMatrix\")")

  # Get the number of rows (n) and columns (m) of the data matrix,
  n <- nrow(X)
  m <- ncol(X)

  # Check "clusters" input.
  if (!(is.factor(clusters) & length(clusters) == n))
    stop("Input argument \"clusters\" should be a factor with one entry ",
         "for each row of \"X\"")
  if (any(table(clusters) == 0))
    stop("Each level must appear at least once in factor \"clusters\"")

  # Initialize the loadings matrix from the clustering.
  k <- nlevels(clusters)
  L <- matrix(0,n,k)
  rownames(L) <- rownames(X)
  colnames(L) <- levels(clusters)
  for (j in levels(clusters)) {
    i      <- which(clusters == j)
    L[i,j] <- 1
  }
  L <- rowSums(X) * L

  # Initialize the factors matrix; the MLEs are available in
  # closed-form in this case.
  F <- matrix(0,m,k)
  rownames(F) <- colnames(X)
  colnames(F) <- levels(clusters)
  for (j in levels(clusters)) {
    i     <- which(clusters == j)
    F[,j] <- colSums(X[i,])/sum(L[i,j])
  }

  # Output the Poisson NMF fit object.
  return(init_poisson_nmf(X,F,L,...))
}
