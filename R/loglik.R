#' @rdname altsqp
#' 
#' @export
#' 
loglik_poisson_topic_model <- function (X, fit, e = 1e-15) {

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

  # Verify and process input matrix L.
  verify.matrix(fit$L)
  if (!is.matrix(L))
    L <- as.matrix(L)
  if (is.integer(L))
    storage.mode(L) <- "double"

  # Check that matrices X, F and L are compatible.
  if (!(nrow(L) == nrow(X) &
        nrow(F) == ncol(X) &
        ncol(L) == ncol(F)))
    stop("Dimensions of input arguments \"X\", \"fit$F\" and/or \"fit$L\" ",
         "do not agree")
  
  # Compute the log-likelihood for the Poisson topic model after
  # removing terms that do not depend on L or F.
  return(-cost(X,L,t(F),e,TRUE))
}

#' @rdname altsqp
#' 
#' @export
#' 
loglik_multinom_topic_model <- function (X, fit, e = 1e-15) {

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
  eps <- 2*.Machine$double.eps
  if (any((colSums(F) - 1) > eps))
    stop("Each column of input matrix \"fit$F\" should sum to 1")
  if (!is.matrix(F))
    F <- as.matrix(F)
  if (is.integer(F))
    storage.mode(F) <- "double"

  # Verify and process input matrix L.
  verify.matrix(fit$L)
  if (any((rowSums(L) - 1) > eps))
    stop("Each row of input matrix \"fit$L\" should sum to 1")
  if (!is.matrix(L))
    L <- as.matrix(L)
  if (is.integer(L))
    storage.mode(L) <- "double"

  # Check that matrices X, F and L are compatible.
  if (!(nrow(L) == nrow(X) &
        nrow(F) == ncol(X) &
        ncol(L) == ncol(F)))
    stop(paste("Dimensions of input arguments \"X\", \"fit$F\" and/or",
               "\"fit$L\ do not agree"))
  
  # Compute the log-likelihood for the multinomial topic model after
  # removing terms that do not depend on L or F.
  return(-cost(X,L,t(F),e,FALSE))
}

#' @title Cost Function for Non-negative Matrix Factorization and
#'   Topic Modeling
#'
#' @description Compute negative log-likelihoods for assessing
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
#' @param A The n x k matrix of loadings (also called "activations").
#'   It should be a dense matrix; that is, \code{is.matrix(A)} should
#'   be \code{TRUE}.
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
#' @param model If \code{model = "poisson"}, Poisson log-likelihoods
#'   are returned; if \code{model = "multinom"}, multinomial
#'   log-likelihoods are returned. See "Value" for details.
#'
#' @param const TO DO: Add description here.
#' 
#' @param version If \code{version == "R"}, the computations are
#'   performed entirely in R; if \code{version == "Rcpp"}, an Rcpp
#'   implementation is used.  The R version is typically faster, but the
#'   Rcpp version makes much more efficient use of memory for large,
#'   sparse matrices.
#'
#' @return The return value is a vector with one entry per row of
#'   X. When \code{poisson = TRUE}, this vector contains the negative
#'   Poisson log-likelihoods; when \code{poisson = FALSE}, this vector
#'   contains the negative multinomial log-likelihoods.
#'
#' @seealso 
#'
#' @examples 
#' 
#' @export
#'
#' @keywords internal
#' 
cost <- function (X, A, B, e = 1e-15, model = c("poisson", "multinom"),
                  const = 0, version) {
  model   <- match.arg(model)
  poisson <- model == "poisson"
  if (missing(version)) {
    if (is.matrix(X))
      version <- "R"
    else
      version <- "Rcpp"
  }
  if (version == "R") {
    AB <- A %*% B
    f  <- rowSums(poisson*AB - X*log(AB + e))
  } else if (is.matrix(X))
    f <- cost_rcpp(X,A,B,e,poisson)
  else
    f <- cost_sparse_rcpp(X,A,B,e,poisson)
  return(f + const)
}
