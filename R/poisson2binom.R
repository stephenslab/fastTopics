#' @title Recover Binomial Topic Model Fit from Poisson NMF fit
#'
#' @description Approximately recover the maximum-likelihood estimates
#' of the binomial topic model from a Poisson NMF model fit (see
#' \code{\link{fit_poisson_nmf}}). The Poisson NMF model is a good
#' approximation when the binary matrix \code{X} is large and ones are
#' infrequent.
#'
#' @details The binomial topic model for binary matrix \eqn{X} with
#' entries \eqn{x_{ij}} is \deqn{x_{ij} \sim \mathrm{Binom}(1,
#' p_{ij})} where the binomial probabilities are given by the linear
#' combinations \deqn{p_{ij} = \sum_{k=1}^K l_{ik} f_{jk},} in which
#' the the \eqn{l_{ik}}'s are the topic proportions and the
#' \eqn{f_{jk}}'s are the topic-specific frequencies of ones.
#'
#' Note the special case of the binomial with sample size of 1 is also
#' called the Bernoulli distribution.
#' 
#' @param X The n x m \dQuote{binary} matrix; all entries of X should
#'   be between 0 and 1 (including 0 and 1). It can be a sparse matrix
#'   (class \code{"dgCMatrix"}) or dense matrix (class \code{"matrix"}).
#'
#' @param fit An object of class \dQuote{poisson_nmf_fit}, such as an
#'   output from \code{fit_poisson_nmf}. It does not make sense for a
#'   binomial topic model to have less than 2 topics, so an error
#'   will be reported when k < 2, where k is the rank of the matrix
#'   factorization.
#'
#' @param numem Number of EM refinement steps. Increasing this number
#'   may improve the quality of the estimates.
#'
#' @param umin A small positive number near zero used to prevent
#'   division by zero in the rescaling of the Poisson NMF factorization.
#' 
#' @param verbose When \code{verbose = TRUE}, progress information is
#'   printed to the console.
#' 
#' @return The return value is the list \code{fit}, in which
#'   \code{fit$F} and \code{fit$L} are the parameters of the binomial
#'   topic model; specifically, \code{fit$L[i,]} gives the topic
#'   proportions for document/sample i, and \code{fit$F[,k]} gives
#'   the frequencies for topic k.
#' 
#' @examples
#' # See the vignette for an example.
#' 
#' @export
#'
poisson2binom <- function (X, fit, numem = 0, umin = 1e-4, verbose = TRUE) {
  if (!requireNamespace("NNLM",quietly = TRUE))
    stop("poisson2binom requires the NNLM package")

  # Check input argument "fit".
  if (inherits(fit,"binom_topic_model_fit"))
    return(fit)
  if (!inherits(fit,"poisson_nmf_fit"))
    stop("Input argument \"fit\" should be an object of class ",
         "\"poisson_nmf_fit\"")
  verify.fit(fit)
  if (ncol(fit$F) < 2 | ncol(fit$L) < 2)
    stop("Input matrices \"fit$F\" and \"fit$L\" should have 2 or more",
         "columns")

  # Check and process input argument "X".
  if (!((is.numeric(X) & is.matrix(X)) | is.sparse.matrix(X)))
    stop("Input argument \"X\" should be a numeric matrix (a \"matrix\" or ",
         "a \"dgCMatrix\")")
  verify.fit.and.count.matrix(X,fit)
  if (any(X < 0) | any(X > 1))
    warning("Input argument \"X\" should be a \"binary\" matrix ",
            "(that is, all entries should range from 0 and 1)")
  if (is.matrix(X) & is.integer(X))
    storage.mode(X) <- "double"

  # Choose U = diag(u) such that L*U is closer to being a matrix of
  # topic proportions.
  ones <- matrix(1,n,1)
  L <- fit$L
  F <- fit$F
  if (verbose)
    cat("Rescaling L and F using non-negative linear regression (nnlm).\n")
  u   <- drop(coef(NNLM::nnlm(L,ones)))
  u   <- pmax(u,umin)
  L   <- scale.cols(L,u)
  L   <- normalize.rows(L)
  F   <- scale.cols(F,1/u)
  n   <- nrow(L)
  fit <- list(F = F,L = L,s = rep(1,n),progress = NA)
  names(fit$s) <- rownames(L)
  
  # Refine the binomial topic model fit by performing several EM updates.
  if (numem > 0) {
    cat("Performing",numem,"EM updates to refine the fit.\n")
    progress <- as.data.frame(cbind(1:numem,0,0))
    names(progress) <- c("iter","delta.f","delta.l")
    if (verbose)
    cat("iter  |F - F'|  |L - L'|\n")
    for (i in 1:numem) {
      fit0 <- fit
      fit  <- fit_binom_topic_model_em(X,fit,numem)
      progress[i,"delta.f"] <- max(abs(fit0$F - fit$F))
      progress[i,"delta.l"] <- max(abs(fit0$L - fit$L))
      if (verbose)
        cat(sprintf("%4d %0.3e %0.3e\n",i,progress[i,"delta.f"],
                    progress[i,"delta.f"]))
    }
    fit$progress <- progress
  }
  
  # Return the Binomial topic model fit.
  class(fit) <- c("binom_topic_model_fit","list")
  return(fit)
}

# Perform a single EM udpate for fiitting the binomial topic model to
# binary data matrix X. This code is adapted from the meth_tpxEM
# function in the methClust package by Kushal Dey.
fit_binom_topic_model_em <- function (X, fit, numiter) {
  if (!is.matrix(X))
    X <- as.matrix(X)
    
  # Make sure no parameters are exactly zero or exactly one.
  e <- 1e-8
  L <- fit$L
  F <- fit$F
  F <- clamp(F,e,1 - e)
  L <- clamp(L,e,1 - e)
  L <- normalize.rows(L)

  # Perform the E step.
  A  <- X/tcrossprod(L,F)
  M  <- (A %*% F) * L
  Mt <- crossprod(A,L) * F
  A  <- (1 - X)/tcrossprod(L,1 - F)
  U  <- (A %*% (1 - F)) * L
  Ut <- crossprod(A,L) * (1 - F)

  # Perform the M step.
  L <- normalize.rows(M + U)
  F <- Mt/(Mt + Ut)
  return(list(F = F,L = L))
}
