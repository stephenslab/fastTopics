#' @title Add Title Here.
#'
#' @description Add description here.
#'
#' @param X Describe input argument "X" here.
#' 
#' @param fit0 Describe input argument "fit0" here.
#'
#' @param estimate.alpha Describe input argument "estimate.alpha" here.
#'
#' @param iter.max Describe input argument "numiter" here.
#' 
#' @param minval Describe input argument "minval" here.
#'
#' @param control Describe input argument "control" here.
#'
#' @param \dots Describe input arguments "..." here.
#' 
#' @examples
#' library(Matrix)
#' set.seed(1)
#' X <- simulate_count_data(80,100,k = 4,sparse = TRUE)$X
#' fit <- fit_topic_model(X,k = 4)
#' lda <- fit_lda(X,fit)
#'
#' @importFrom tm as.DocumentTermMatrix
#' @importFrom topicmodels LDA
#' 
#' @export
#' 
fit_lda <- function (X, fit0, estimate.alpha = TRUE, iter.max = 20,
                     minval = 1e-8,
                     control = list(estimate.alpha = estimate.alpha,
                                    alpha = 1,verbose = 1,keep = 1,
                                    em = list(iter.max = iter.max,tol = 0),
                                    var = list(iter.max = 10,tol = 0)), 
                     ...) {
    
  # Check and process input argument "fit0", if provided.
  if (!(inherits(fit0,"poisson_nmf_fit") |
      inherits(fit0,"multinom_topic_model_fit")))
  stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
       "\"multinom_topic_model_fit\"")
  if (inherits(fit0,"poisson_nmf_fit"))
    fit0 <- poisson2multinom(fit0)
  fit0$F <- normalize.cols(pmax(fit0$F,minval))
  fit0$L <- normalize.rows(pmax(fit0$L,minval))

  # Get the number of topics.
  k <- ncol(fit0$F)

  # Check and process input argument "X".
  if (!((is.numeric(X) & is.matrix(X)) | is.sparse.matrix(X)))
    stop("Input argument \"X\" should be a numeric matrix (a \"matrix\" or ",
         "a \"dgCMatrix\")")
  verify.fit.and.count.matrix(X,fit0)
  X <- as.DocumentTermMatrix(X,weighting = c("term frequency","tf"))

  # Initialize the LDA model.
  control0 <- control
  control0$em <- list(iter.max = 4,tol = 0)
  control0$estimate.alpha <- FALSE
  control0$verbose <- 0
  if (control$verbose > 0)
    cat("Initializing LDA model.\n")
  lda0       <- LDA(X,k,control = control0,...)
  lda0@beta  <- t(log(fit0$F))
  lda0@gamma <- fit0$L
  
  # Fit the LDA model.
  if (control$verbose > 0)
      cat("Fitting LDA model.\n")
  lda <- LDA(X,k,model = lda0,control = control,...)
  return(lda)
}
