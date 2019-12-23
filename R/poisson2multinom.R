#' @title Get Multinomial Topic Model from Poisson Non-Negative Matrix
#'   Factorization
#'
#' @description This function recovers parameter estimates of the
#'   multinomial topic model given parameter estimates for a Poisson
#'   non-negative matrix factorization.
#'
#' @param fit A list containing two dense, non-negative matrices,
#'   \code{fit$F} and \code{fit$L}. The former is an m x k matrix of
#'   factors, and the latter is an n x k matrix of loadings. It does not
#'   make sense for a multinomial topic model to have less than two
#'   topics, so an error will be reported when k < 2.
#'
#' @return The return value is the list \code{fit}, in which
#'   \code{fit$F} and \code{fit$L} are the parameters of the multinomial
#'   topic model; specifically, \code{fit$L[i,]} gives the topic
#'   probabilities for sample or document i, and \code{fit$F[,k]} gives
#'   the term probabilities for topic k. An additional vector
#'   \code{fit$s} of length n is returned giving the "document sizes".
#' 
#' @export
#' 
poisson2multinom <- function (fit) {

  # Verify input argument "fit".
  verify.fit(fit)
  if (inherits(fit,"multinom_topic_model"))
    stop("Input argument \"fit\" should not be a multinomial topic model")
  F <- fit$F
  L <- fit$L
  
  # Verify and process input matrix F.
  if (any(colSums(F) <= 0))
    stop("Each column of \"fit$F\" should have at least one positive entry")
  if (is.integer(F))
    storage.mode(F) <- "double"
  
  # Verify and process input matrix L. Each row of L should have at
  # least one positive entry.
  if (any(rowSums(L) <= 0))
    stop("Each row of \"fit$L\" should have at least one positive entry")
  if (is.integer(L))
    storage.mode(L) <- "double"
  
  # Check that k > 1.
  if (ncol(F) < 2)
    stop("Input matrices \"fit$F\" and \"fit$L\" should have at least ",
         "2 columns")
  
  # Recover F and L for the multinomial model. Here, s gives the
  # Poisson rates for generating the "document sizes".
  L <- scale.cols(L,colSums(F))
  s <- rowSums(L)
  L <- L / s
  F <- normalize.cols(F)

  # Update the "fit" object, and return it.
  fit$F <- F
  fit$L <- L
  fit$s <- s
  class(fit) <- c("multinom_topic_model","list")
  return(fit)
}
