#' @rdname likelihood
#'
#' @title NMF and Topic Model Likelihoods and Deviances
#'
#' @description Compute log-likelihoods and deviances for assessing
#'   fit of a topic model or a non-negative matrix factorization (NMF).
#'
#' @details Function \code{cost} computes loss functions proportional
#'   to the negative log-likelihoods, and is mainly for internal use to
#'   quickly compute log-likelihoods and deviances; it should not be
#'   used directly unless you know what you are doing. In particular,
#'   little argument checking is performed by \code{cost}.
#'
#' @param X The n x m matrix of counts or pseudocounts. It can be a
#'   sparse matrix (class \code{"dgCMatrix"}) or dense matrix (class
#'   \code{"matrix"}).
#'
#' @param fit A Poisson NMF or multinomial topic model fit, such as an
#'   output from \code{\link{fit_poisson_nmf}} or
#'   \code{\link{fit_topic_model}}. 
#'
#' @param A The n x k matrix of loadings. It should be a dense matrix.
#'
#' @param B The k x m matrix of factors. It should be a dense matrix.
#'
#' @param e A small, non-negative number added to the terms inside the
#'   logarithms to avoid computing logarithms of zero. This prevents
#'   numerical problems at the cost of introducing a very small
#'   inaccuracy in the computation.
#'
#' @param family If \code{model = "poisson"}, the loss function values
#'   corresponding to the Poisson non-negative matrix factorization are
#'   computed; if \code{model = "multinom"}, the multinomial topic model
#'   loss function values are returned.
#'
#' @param version When \code{version == "R"}, the computations are
#'   performed entirely in R; when \code{version == "Rcpp"}, an Rcpp
#'   implementation is used. The R version is typically faster when
#'   \code{X} is a dense matrix, whereas the Rcpp version is faster and
#'   more memory-efficient when \code{X} is a large, sparse matrix. When
#'   not specified, the most suitable version is called depending on
#'   whether \code{X} is dense or sparse.
#'
#' @return A numeric vector with one entry per row of \code{X}.
#'
#' @examples
#'
#' # Generate a small counts matrix.
#' set.seed(1)
#' out <- simulate_count_data(10,20,3)
#' X   <- out$X
#' fit <- out[c("F","L")]
#' class(fit) <- c("poisson_nmf_fit","list")
#' 
#' # Compute the Poisson log-likelihoods and deviances.
#' data.frame(loglik   = loglik_poisson_nmf(X,fit),
#'            deviance = deviance_poisson_nmf(X,fit))
#' 
#' # Compute multinomial log-likelihoods.
#' loglik_multinom_topic_model(X,fit)
#' 
#' @export
#' 
loglik_poisson_nmf <- function (X, fit, e = 1e-8)
  loglik_helper(X,fit,"loglik.poisson",e)

#' @rdname likelihood
#' 
#' @export
#' 
loglik_multinom_topic_model <- function (X, fit, e = 1e-8)
  loglik_helper(X,fit,"loglik.multinom",e)

#' @rdname likelihood
#' 
#' @export
#' 
deviance_poisson_nmf <- function (X, fit, e = 1e-8)
  loglik_helper(X,fit,"deviance.poisson",e)

#' @rdname likelihood
#'
#' @importFrom Rcpp evalCpp
#' 
#' @export
#'
cost <- function (X, A, B, e = 1e-8, family = c("poisson","multinom"),
                  version) {

  # For the special case of a rank-1 factorization, A and B may need
  # to be vectors, in which case we need to force them to be matrices.
  if (!is.matrix(A))
    A <- matrix(A)
  if (!is.matrix(B))
    B <- matrix(B,1,ncol(X))
    
  # Check and process "model" and "version" input arguments.
  family  <- match.arg(family)
  poisson <- family == "poisson"
  if (missing(version)) {
    if (is.matrix(X))
      version <- "R"
    else
      version <- "Rcpp"
  }

  # Compute the terms in the log-likelihoods that depend on A or B.
  if (version == "R") {
    AB <- A %*% B
    f  <- rowSums(poisson*AB - X*log(AB + e))
  } else if (version == "Rcpp") {
    if (is.matrix(X))
      f <- drop(cost_rcpp(X,A,B,e,poisson))
    else
      f <- drop(cost_sparse_rcpp(X,A,B,e,poisson))
  }
  names(f) <- rownames(X)
  return(f)
}

# This function provides the core implementation for calculation of
# log-likelihoods and deviances.
loglik_helper <- function (X, fit,
  output.type = c("loglik.poisson","loglik.multinom","deviance.poisson"), e) {

  # Verify and process input "output.type".
  output.type <- match.arg(output.type)
      
  # Verify and process inputs "X" and "fit". 
  if (!(inherits(fit,"poisson_nmf_fit") |
        inherits(fit,"multinom_topic_model_fit")))
    stop("Input argument \"fit\" should be an object of class ",
         "\"multinom_topic_model_fit\" or \"multinom_topic_model_fit\"")
  if (inherits(fit,"multinom_topic_model_fit"))
    fit <- multinom2poisson(fit)
  verify.fit.and.count.matrix(X,fit)
  if (is.matrix(X) & is.integer(X))
    storage.mode(X) <- "double"
  F <- fit$F
  L <- fit$L

  # Compute the log-likelihood or deviance.
  if (output.type == "loglik.poisson" | output.type == "loglik.multinom") {
    f <- loglik_poisson_const(X) - cost(X,L,t(F),e,"poisson")
    if (output.type == "loglik.multinom")
      f <- f - loglik_size_factors(X,fit$F,fit$L)
  } else if (output.type == "deviance.poisson")
    f <- deviance_poisson_const(X) + 2*cost(X,L,t(F),e,"poisson")
  return(f)
}

# Compute the constant terms in the Poisson log-likelihoods.
#
#' @importFrom Matrix rowSums
#
loglik_poisson_const <- function (X) {
  if (is.matrix(X))
    return(-rowSums(lgamma(X + 1)))
  else
    return(-rowSums(apply.nonzeros(X,function (x) lgamma(x + 1))))
}

# Compute the constant terms in the multinomial log-likelihoods.
loglik_multinom_const <- function (X)
  lgamma(rowSums(X) + 1) + loglik_poisson_const(X)

# Compute the constant terms in the Poisson deviances.
#
#' @importFrom Matrix rowSums
#
deviance_poisson_const <- function (X) {
  e <- .Machine$double.eps
  if (is.matrix(X))
    return(2*rowSums(X*(log(X + e) - 1)))
  else
    return(2*rowSums(apply.nonzeros(X,function (x) x*(log(x) - 1))))
}

# Compute the residuals of the first-order Karush-Kuhn-Tucker (KKT)
# conditions for Poisson non-negative matrix factorization at solution
# estimate (F,L).
poisson_nmf_kkt <- function (X, F, L, e = 1e-8) {
  if (is.matrix(X)) {
    A <- X/(tcrossprod(L,F) + e)
    return(list(F = F*(t(1 - A) %*% L),
                L = L*((1 - A) %*% F)))
  } else {
    A <- x_over_tcrossprod(X,t(L),t(F),e)
    return(list(F = F*(repmat(colSums(L),ncol(X)) - as.matrix(t(A) %*% L)),
                L = L*(repmat(colSums(F),nrow(X)) - as.matrix(A %*% F))))
  }
}

# Given a Poisson non-negative matrix factorization (F, L), compute
# the log-likelihoods for the "size factors"; that is, each vector
# element is a log-likelihood for the model t ~ Poisson(s), where
# t = sum(x) is the total sum of the counts. The size factors, s, are
# recovered from the poisson2multinom transformation.
#
#' @importFrom stats dpois
loglik_size_factors <- function (X, F, L)
  ldpois(rowSums(X),get_multinom_from_pnmf(F,L)$s)

# Compute the log-density for the Poisson distribution with rate
# lambda at x. Inputs x and lambda may be vectors, and x need not be
# integer-valued.
ldpois <- function (x, lambda)
  x*log(lambda) - lambda - lfactorial(x)
