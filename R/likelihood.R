#' @rdname likelihood
#'
#' @title Compute Topic Model Likelihoods 
#' 
#' @export
#' 
loglik_poisson_topic_model <- function (X, fit, e = 1e-15)
  loglik_topic_model_helper(X,fit,"loglik.pois",e)

#' @rdname likelihood
#' 
#' @export
#' 
loglik_multinom_topic_model <- function (X, fit, e = 1e-15)
  loglik_topic_model_helper(X,fit,"loglik.multinom",e)

#' @rdname likelihood
#' 
#' @export
#' 
deviance_poisson_topic_model <- function (X, fit, e = 1e-15) {
  # TO DO.
}

#' @rdname likelihood
#' 
#' @description Loss Function for Non-negative Matrix Factorization and
#'   Topic Modeling
#'
#' Compute negative log-likelihoods for assessing
#'   a topic model fit or quality of a non-negative matrix factorization,
#'   in which matrix X is approximated by matrix product A * B.
#'
#'   This function is mainly used internally to quickly compute
#'   log-likelihoods for model fits and objective values, and should not
#'   be used directly unless you know what you are doing. In particular,
#'   very little argument checking is performed; if you use this
#'   function, it is up to you to make sure you use it correctly.
#'
#' @param X The n x m matrix of counts or pseudocounts. It can be a
#'   sparse matrix ("dgCMatrix" class) or dense matrix ("matrix" class).
#'
#' @param A The n x k matrix of loadings (also called "activations"),
#'   where k is typically much smaller than n and m. A should be a
#'   dense matrix; that is, \code{is.matrix(A)} should be \code{TRUE}.
#'
#' @param B The k x m matrix of factors (also called "basis vectors").
#'   It should be a dense matrix; that is, \code{is.matrix(B)} should be
#'   \code{TRUE}.
#'
#' @param e A small, non-negative number added to the terms inside the
#'   logarithms to avoid computing logarithms of zero. This prevents
#'   numerical problems at the cost of introducing a very small
#'   inaccuracy in the computation.
#'
#' @param family If \code{model = "poisson"}, Poisson log-likelihoods
#'   are returned; if \code{model = "multinom"}, multinomial
#'   log-likelihoods are returned. See "Value" for details.
#'
#' @param const To avoid re-calculating terms in the log-likelihood
#'   that do not depend on either A or B, one may provide a value for
#'   this argument. Typically this argument is not provided, in which
#'   case the constant terms are calculated internally.
#' 
#' @param version If \code{version == "R"}, the computations are
#'   performed entirely in R; if \code{version == "Rcpp"}, an Rcpp
#'   implementation is used. The R version is typically faster when
#'   \code{X} is a dense matrix, whereas the Rcpp version is faster and
#'   more memory-efficient when \code{X} is a large, sparse matrices.
#'
#' @return A numeric vector with one entry per row of \code{X}. When
#'   \code{poisson = TRUE}, this vector contains the negative Poisson
#'   log-likelihoods; when \code{poisson = FALSE}, this vector contains
#'   the negative multinomial log-likelihoods.
#'
#' @seealso 
#'
#' @export
#'
#' @keywords internal
#' 
cost <- function (X, A, B, e = 1e-15, family = c("poisson","multinom"),
                  const, version) {

  # Check and process "model" and "version" input arguments.
  family  <- match.arg(family)
  poisson <- model == "poisson"
  if (missing(version)) {
    if (is.matrix(X))
      version <- "R"
    else
      version <- "Rcpp"
  }

  # When not provided, compute the terms in the log-likelihoods that
  # do not depend on A or B.
  if (missing(const)) {
    if (poisson)
      const <- loglik_poisson_const(X)
    else
      const <- loglik_multinom_const(X)
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
  return(f + const)
}

# This function provides the core implementation for the
# loglik_poisson_topic_model, loglik_multinom_topic_model and
# deviance_poisson_topic_model functions.
loglik_topic_model_helper <-
  function (X, fit,
            output.type = c("loglik.pois","loglik.mult","dev.pois"), e) {

  # Verify and process input "output.type".
  output.type <- match.arg(output.type)
      
  # Verify and process input matrix X. 
  verify.matrix(X)
  if (is.matrix(X) & is.integer(X))
    storage.mode(X) <- "double"
  
  # Input argument "fit" should be a list with elements "F" and "L".
  verify.fit(fit)
  F <- fit$F
  L <- fit$L

  # Verify and process input matrix F.
  verify.matrix(fit$F)
  if (!is.matrix(F))
    F <- as.matrix(F)
  if (is.integer(F))
    storage.mode(F) <- "double"
  if (output.type == "loglik.mult")
    F <- normalize.cols(F)
  
  # Verify and process input matrix L.
  verify.matrix(fit$L)
  if (!is.matrix(L))
    L <- as.matrix(L)
  if (is.integer(L))
    storage.mode(L) <- "double"
  if (output.type == "loglik.mult")
    L <- normalize.cols(L)

  # Verify and process input matrix L.
  if (any((rowSums(L) - 1) > eps))
    stop("Each row of input matrix \"fit$L\" should sum to 1")
  
  # Check that matrices X, F and L are compatible.
  if (!(nrow(L) == nrow(X) &
        nrow(F) == ncol(X) &
        ncol(L) == ncol(F)))
    stop("Dimensions of input arguments \"X\", \"fit$F\" and/or \"fit$L\" ",
         "do not agree")
  
  # Compute the log-likelihood or deviance for the topic model.
  if (output.type == "loglik.pois")
    return(-cost(X,L,t(F),e,TRUE))
  else if (output.type == "loglik.mult")
    return(-cost(X,L,t(F),e,FALSE))
  else if (output.type == "dev.pois") {
    # TO DO.
  }
}

# Compute the constant terms in the Poisson log-likelihoods. This is
# used by the "cost" function, but can be also used elsewhere to
# pre-compute these constants.
loglik_poisson_const <- function (X) {
  if (is.matrix(X))
    return(rowSums(lgamma(X + 1)))
  else
    return(rowSums(apply.nonzeros(X,function (x) lgamma(x + 1))))
}

# Compute the constant terms in the multinomial log-likelihoods. This
# is used by the "cost" function, but can be also used elsewhere to
# pre-compute these constants.
loglik_multinom_const <- function (X)
  loglik_poisson_const(X) - lgamma(rowSums(X) + 1)

