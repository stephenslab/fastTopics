#' @title Recover Multinomial Topic Model Fit from Poisson NMF fit
#'
#' @description This function recovers parameter estimates of the
#'   multinomial topic model given parameter estimates for a Poisson
#'   non-negative matrix factorization (NMF).
#'
#' @param fit An object of class \dQuote{poisson_nmf_fit}, such as an
#'   output from \code{fit_poisson_nmf}. It does not make sense for a
#'   multinomial topic model to have less than two topics, so an error
#'   will be reported when k < 2, where k is the rank of the matrix
#'   factorization.
#'
#' @return The return value is the list \code{fit}, in which
#'   \code{fit$F} and \code{fit$L} are the parameters of the multinomial
#'   topic model; specifically, \code{fit$L[i,]} gives the topic
#'   probabilities for sample or document i, and \code{fit$F[,k]} gives
#'   the term probabilities for topic k. An additional vector
#'   \code{fit$s} of length n is returned giving the "size factors".
#' 
#' @export
#' 
poisson2multinom <- function (fit) {

  # Check input argument "fit".
  if (!inherits(fit,"poisson_nmf_fit"))
    stop("Input argument \"fit\" should be an object of class ",
         "\"poisson_nmf_fit\"")
  verify.fit(fit)
  if (ncol(fit$F) < 2 | ncol(fit$L) < 2)
    stop("Input matrices \"fit$F\" and \"fit$L\" should have 2 or more",
         "columns")
  
  # Recover F and L for the multinomial model. Here, s gives the
  # Poisson rates for generating the "document sizes".
  out   <- get_multinom_from_pnmf(fit$F,fit$L)
  fit$F <- out$F
  fit$L <- out$L
  fit$s <- out$s
  
  # Return the updated fit.
  class(fit) <- c("multinom_topic_model_fit","list")
  return(fit)
}

# Get the parameters of the multinomial topic model from the
# parameters of the Poisson NMF model.
get_multinom_from_pnmf <- function (F, L) {
  u <- colSums(F)
  F <- scale.cols(F,1/u)
  L <- scale.cols(L,u)
  s <- rowSums(L)
  L <- L / s
  return(list(F = F,L = L,s = s))
}
