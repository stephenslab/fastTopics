#' @title Simple Interface for Fitting a Multinomial Topic Model
#'
#' @description Fits a multinomial topic model to the count data,
#'   hiding most of the complexities of model fitting. The default
#'   optimization settings used here are intended to work well in a wide
#'   range of data sets, although some fine-tuning may be needed for
#'   more difficult cases. For full control, use \code{fit_poisson_nmf}.
#'
#' @details With the default settings, the model fitting is
#' accomplished in four steps: (1) initialize the Poisson NMF model
#' fit (\code{init_poisson_nmf}); (2) perform the main model fitting
#' step by running 100 EM updates using \code{fit_poisson_nmf}; (3)
#' refine the fit by running 100 extrapolated SCD updates, again using
#' \code{fit_poisson_nmf}; and (4) recover the multinomial topic model
#' by calling \code{poisson2multinom}.
#'
#' This two-stage fitting approach is based on our findings that the
#' EM algorithm initially makes rapid progress toward a solution, but
#' its convergence slows considerably as the iterates approach a
#' solution. Close to a solution, we have found that other algorithms
#' make much more rapid progress than EM; in particularly, we found
#' that the extrapolated SCD updates usually performed best). For
#' larger data sets, more updates in the main model fitting and
#' refinement steps may be needed to obtain a good fit.
#'
#' For larger data sets, more than 200 updates may be needed to obtain
#' a good fit.
#'
#' @param X The n x m matrix of counts; all entries of X should be
#'   non-negative. It can be a sparse matrix (class \code{"dgCMatrix"})
#'   or dense matrix (class \code{"matrix"}).
#'
#' @param k The number of topics. Must be 2 or greater.
#'
#' @param numiter.main Number of updates of the factors and loadings
#'   to perform in the main model fitting step. Should be 1 or more.
#'
#' @param numiter.refine Number of updates of the factors and loadings
#'   to perform in the model refinement step.
#'
#' @param method.main The method to use for updating the factors and
#'   loadings in the main model fitting step. Passed as argument
#'   "method" to \code{\link{fit_poisson_nmf}}.
#'
#' @param method.refine The method to use for updating the factors in
#'   evthe model refinement step. Passed as argument "method"
#'   to \code{\link{fit_poisson_nmf}}.
#'
#' @param init.method The method used to initialize the factors and
#'   loadings. See \code{\link{init_poisson_nmf}} for details.
#'
#' @param control.init A list of parameters controlling the behaviour
#'   of the optimization and Topic SCORE method in the call to
#'   \code{init_poisson_nmf}. This is passed as argument "control" to
#'   \code{init_poisson_nmf}.
#'
#' @param control.main A list of parameters controlling the behaviour
#'   of the optimization in the main model fitting step. This is passed
#'   as argument "control" to \code{fit_poisson_nmf}.
#'
#' @param control.refine A list of parameters controlling the
#'   behaviour of the of the optimization algorithm in the model
#'   refinement step. This is passed as argument "control" to
#'   \code{fit_poisson_nmf}.
#'
#' @param verbose When \code{verbose = "progressbar"} or \code{verbose
#'   = "detailed"}, information about the progress of the model fitting
#'   is printed to the console. See \code{\link{fit_poisson_nmf}} for
#'   more information.
#' 
#' @seealso \code{\link{init_poisson_nmf}},
#'   \code{\link{fit_poisson_nmf}},
#'   \code{\link{poisson2multinom}},
#'   \code{\link{fit_multinom_model}}
#' 
#' @return A multinomial topic model fit; see
#'   \code{\link{poisson2multinom}} and \code{\link{fit_poisson_nmf}}
#'   for details. Note that outputted likelihoods and deviances in
#'   \code{progress} are evaluated with respect to the equivalent
#'   Poisson NMF model.
#'
#' @references
#' Dey, K. K., Hsiao, C. J. and Stephens, M. (2017). Visualizing the
#' structure of RNA-seq expression data using grade of membership
#' models. \emph{PLoS Genetics} \bold{13}, e1006599.
#'
#' Blei, D. M., Ng, A. Y. and Jordan, M. I. (2003). Latent Dirichlet
#' allocation. \emph{Journal of Machine Learning Research} \bold{3},
#' 993-1022.
#'
#' Hofmann, T. (1999). Probabilistic latent semantic indexing. In
#' \emph{Proceedings of the 22nd International ACM SIGIR Conference},
#' 50-57. doi:10.1145/312624.312649
#' 
#' @examples
#' library(Matrix)
#' set.seed(1)
#' X <- simulate_count_data(80,100,k = 3,sparse = TRUE)$X
#' fit <- fit_topic_model(X,k = 3)
#' print(summary(fit))
#' 
#' @export
#'
fit_topic_model <-
  function (X, k, numiter.main = 100, numiter.refine = 100, method.main = "em",
            method.refine = "scd", init.method = c("topicscore","random"),
            control.init = list(), control.main = list(numiter = 4),
            control.refine = list(numiter = 4,extrapolate = TRUE),
            verbose = c("progressbar","detailed","none")) {
      
  # Check the input data matrix.
  verify.count.matrix(X)
  
  # Check and process input argument "verbose".
  verbose <- match.arg(verbose)
  
  # If necessary, remove all-zero columns from the counts matrix.
  if (any.allzero.cols(X)) {
    X <- remove.allzero.cols(X)
    warning(sprintf(paste("One or more columns of X are all zero; after",
                          "removing all-zero columns, %d columns will be",
                          "used for model fitting"),ncol(X)))
  }

  # Initialize the Poisson NMF model fit.
  fit <- init_poisson_nmf(X,k = k,init.method = init.method,
                          control = control.init,
                          verbose = ifelse(verbose == "none",
                                           "none","detailed"))

  # Perform the main model fitting step.
  fit <- fit_poisson_nmf(X,fit0 = fit,numiter = numiter.main,
                         method = method.main,control = control.main,
                         verbose = verbose)
  
  # Perform the model refinement step.
  if (numiter.refine > 0) {
    if (verbose != "none")
      cat("Refining model fit.\n")
    fit <- fit_poisson_nmf(X,fit0 = fit,numiter = numiter.refine,
                           method = method.refine,control = control.refine,
                           verbose = verbose)
  }
  
  # Output the multinomial topic model fit.
  return(poisson2multinom(fit))
}
