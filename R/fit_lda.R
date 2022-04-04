#' @title Fit LDA Model Initialized with Maximum-LikeLihood Estimates
#'
#' @description Fit a Latent Dirichlet Allocation (LDA) model using
#'   \code{\link{LDA}} from the \code{topicmodels} package, in which the
#'   model parameters are initialized from maximum-likelihood estimates
#'   of a Poisson NMF or multinomial topic model.
#'
#' @param The n x m counts matrix. It can be a sparse matrix (class
#'   \code{"dgCMatrix"}) or dense matrix (class \code{"matrix"}).
#' 
#' @param fit0 An object of class \dQuote{poisson_nmf_fit} or
#'   \dQuote{multinom_topic_model_fit}. If a Poisson NMF fit is provided
#'   as input, the corresponding multinomial topic model fit is
#'   automatically recovered using \code{\link{poisson2multinom}}.
#' 
#' @param estimate.alpha Indicates if the parameter \dQuote{alpha} is
#'   fixed or estimated; see \code{help("TopicModelcontrol-class")} for
#'   details.
#' 
#' @param iter.max The maximum number of variational EM iteration to
#'   perform; see \code{help("TopicModelcontrol-class")} for
#'   details.
#' 
#' @param minval All model parameters in \code{fit0} smaller than
#'   \code{minval} are set to \code{minval} before initializing the LDA
#'   fit.
#'
#' @param control Control parameters passed to
#'   \code{\link[topicmodels]{LDA}}.
#'
#' @param \dots Additional arguments passed to
#'   \code{\link[topicmodels]{LDA}}.
#' 
#' @examples
#' library(Matrix)
#' library(modeltools)
#' set.seed(1)
#' X <- simulate_count_data(80,100,k = 4,sparse = TRUE)$X
#' fit <- fit_topic_model(X,k = 4)
#' lda <- fit_lda(X,fit)
#'
#' # Dirichlet prior for topic proportions.
#' lda@alpha
#'
#' # Compare LDA and maximum-likelihood estimates.
#' plot(fit$F,t(posterior(lda)$terms),pch = 20)
#' plot(fit$L,posterior(lda)$topics,pch = 20)
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
  if (!(requireNamespace("tm",quietly = TRUE) &
        requireNamespace("topicmodels",quietly = TRUE)))
    stop("mixKWDual requires packages tm and topicmodels")
    
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
  X <- tm::as.DocumentTermMatrix(X,weighting = c("term frequency","tf"))

  # Initialize the LDA model.
  control0 <- control
  control0$em <- list(iter.max = 4,tol = 0)
  control0$estimate.alpha <- FALSE
  control0$verbose <- 0
  if (control$verbose > 0)
    cat("Initializing LDA model.\n")
  lda0       <- topicmodels::LDA(X,k,control = control0,...)
  lda0@beta  <- t(log(fit0$F))
  lda0@gamma <- fit0$L
  
  # Fit the LDA model.
  if (control$verbose > 0)
      cat("Fitting LDA model.\n")
  lda <- topicmodels::LDA(X,k,model = lda0,control = control,...)
  return(lda)
}
