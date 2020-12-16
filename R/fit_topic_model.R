#' @title Simple Interface for Fitting a Multinomial Topic Model
#'
#' @description Fits a multinomial topic model to the count data,
#'   hiding most of the complexities of model fitting. The default
#'   optimization settings used here are intended to work well in a wide
#'   range of data sets, although some fine-tuning may be needed for
#'   more difficult cases.
#'
#' @details TO DO: Add details here.
#'
#' @param X The n x m matrix of counts; all entries of X should be
#'   non-negative. It can be a sparse matrix (class \code{"dgCMatrix"})
#'   or dense matrix (class \code{"matrix"}).
#'
#' @param k The number of topics. Must be 2 or greater.
#'
#' @param init.method The method used to initialize the factors and
#'   loadings. See \code{\link{init_poisson_nmf}} for details.
#'
#' @param control.init A list of parameters controlling the behaviour
#'   of the optimization and Topic SCORE method in the call to
#'   \code{init_poisson_nmf}.
#'
#' @param verbose When \code{verbose = TRUE}, information about the
#'   progress of the model fitting is printed to the console. See
#'   \code{\link{fit_poisson_nmf}} for an explanation of the output.
#' 
#' @seealso \code{\link{init_poisson_nmf}},
#'   \code{\link{fit_poisson_nmf}}, \code{\link{poisson2multinom}}.
#' 
#' @return A multinomial topic model fit; see
#'   \code{\link{poisson2multinom}} and \code{\link{fit_poisson_nmf}}
#'   for details.
#' 
#' @examples
#' library(Matrix)
#' set.seed(1)
#' X <- simulate_count_data(80,100,k = 3,sparse = TRUE)$X
#' 
#' @export
#'
fit_topic_model_simple <- function (X, k,
                                    init.method = c("topicscore","random"),
                                    control.init = list(),
                                    verbose = TRUE) {

  # Check inputs.
  if (!((is.numeric(X) & is.matrix(X)) | is.sparse.matrix(X)))
    stop("Input argument \"X\" should be a numeric matrix (a \"matrix\" or ",
         "a \"dgCMatrix\")")

  # If necessary, remove all-zero columns from the counts matrix.
  if (any(colSums(X > 0) == 0)) {
    i <- which(colSums(X > 0) >= 1)
    X <- X[,i]
    message(sprintf(paste("One or more columns of X are all zero; after",
                          "removing all-zero columns, %d columns will be",
                          "used for model fitting"),length(i)))
  }

  # Initialize the Poisson NMF model fit.
  fit0 <- init_poisson_nmf(X,k = k,init.method = init.method,
                           control = control.init,
                           verbose = verbose)

  # Temporary.
  fit <- fit0

  # Output the multinomial topic model fit.
  return(poisson2multinom(fit))
}
