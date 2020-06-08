#' @title Summarize Poisson NMF or Multinomial Topic Model Fit
#'
#' @description \code{summary} method for the \dQuote{poisson_nmf_fit}
#'   and \dQuote{multinom_topic_model_fit} classes.
#'
#' @param object An object of class \dQuote{poisson_nmf_fit} or
#'   \dQuote{multinom_topic_model_fit}. The former is usually the result
#'   of calling \code{\link{fit_poisson_nmf}}; the latter is usually the
#'   result of calling \code{link{poisson2multinom}}.
#'
#' @param \dots Additional arguments passed to the generic method.
#'
#' @method summary poisson_nmf_fit
#'
#' @return The functions \code{summary.poisson_nmf_fit} and
#'   \code{summary.multinom_topic_model_fit} compute and return a list
#'   of statistics summarizing the model fit.
#' 
#' @export
#' 
summary.poisson_nmf_fit <- function (object, ...) {
  out <- summary(poisson2multinom(object))
  class(out) <- c("summary.poisson_nmf_fit","list")
  return(out)
}

#' @rdname summary.poisson_nmf_fit
#'
#' @method summary multinom_topic_model_fit
#' 
#' @export
#' 
summary.multinom_topic_model_fit <- function (object, ...) {
  #
  # TO DO.
  #  
  class(out) <- c("summary.poisson_nmf_fit","list")
  return(out)
}

#' @rdname summary.poisson_nmf_fit
#'
#' @param x An object of class \dQuote{summary.poisson_nmf_fit},
#'   usually a result of a call to \code{summary.poisson_nmf_fit}.
#'
#' @method print summary.poisson_nmf_fit
#' 
#' @export
#' 
print.summary.poisson_nmf_fit <- function (x, ...) {
  print.summary.multinom_topic_model_fit(x,...)
  return(invisible(x))
}

#' @rdname summary.poisson_nmf_fit
#'  
#' @method print summary.multinom_topic_model_fit
#' 
#' @export
#'
print.summary.multinom_topic_model_fit <- function (x, ...) {
  #
  # TO DO.
  #  
  return(invisible(x))
}
