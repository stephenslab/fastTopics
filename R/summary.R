#' @title Summarize Poisson NMF Fit
#'
#' @description Add description of function here.
#'
#' @param object Describe input argument "object" here.
#'
#' @param \dots Additional arguments passed to the generic method.
#'
#' @method summary poisson_nmf_fit
#'
#' @return The function \code{summary.poisson_nmf_fit} computes and
#'   returns a list of statistics summarizing the model fit.
#' 
#' @export
#' 
summary.poisson_nmf_fit <- function (object, ...) {
  # TO DO.
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
  # TO DO.
  return(invisible(x))
}

