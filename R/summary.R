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
#' @method summary poisson_nmf_fit
#'
#' @return The functions \code{summary.poisson_nmf_fit} and
#' \code{summary.multinom_topic_model_fit} compute and return a list
#' of statistics summarizing the model fit. The returned list
#' includes the following elements:
#'
#' \item{n}{Describe "n" output here.}
#'
#' \item{m}{Describe "m" output here.}
#'
#' \item{k}{Describe "k" output here.}
#'
#' \item{s}{Describe "s" output here.}
#'
#' \item{numiter}{Describe "numiter" output here.}
#'
#' \item{loglik}{Describe "loglik" output here.}
#'
#' \item{dev}{Describe "dev" output here.}
#'
#' \item{res}{Describe "res" output here.}
#'
#' \item{topic.proportions}{Describe "topic.proportions" output here.}
#'
#' \item{topic.reps}{Describe "topic.reps" output here.}
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
#' @importFrom stats quantile
#' 
#' @export
#' 
summary.multinom_topic_model_fit <- function (object, ...) {
  numiter <- nrow(object$progress)
  k       <- ncol(object$F)
  out <- list(n                 = nrow(object$L),
              m                 = nrow(object$F),
              k                 = k,
              s                 = object$s,
              numiter           = numiter,
              loglik            = object$progress[numiter,"loglik"],
              dev               = object$progress[numiter,"dev"],
              res               = object$progress[numiter,"res"],
              topic.proportions = summarize_topic_proportions(object$L),
              topic.reps        = get_topic_representatives(object$L))
  class(out) <- c("summary.multinom_topic_model_fit","list")
  return(out)
}

#' @rdname summary.poisson_nmf_fit
#'
#' @param x An object of class \dQuote{summary.poisson_nmf_fit},
#'   usually a result of a call to \code{summary.poisson_nmf_fit}.
#'
#' @param \dots Additional arguments passed to the generic \code{summary}
#'   or \code{print.summary} method.
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
#' @importFrom stats quantile
#' 
#' @export
#'
print.summary.multinom_topic_model_fit <- function (x, ...) {
  k <- x$k
  cat("Model overview:\n")
  cat(sprintf("  Number of data rows, n: %d\n",x$n))
  cat(sprintf("  Number of data cols, m: %d\n",x$m))
  cat(sprintf("  Rank/Number of topics, k: %d\n",x$k))
  cat(sprintf("Evaluation of fit (%d updates of performed):\n",x$numiter))
  cat(sprintf("  Log-likelihood: %+0.12e\n",x$loglik))
  cat(sprintf("  Deviance: %+0.12e\n",x$dev))
  cat(sprintf("  Max KKT residual: %+0.6e\n",x$res))
  if (!is.null(x$s)) {
    cat(sprintf("Size factors:\n"))
    q        <- quantile(x$s)
    names(q) <- c("Min","1Q","Median","3Q","Max")
    print(q)
  }
  cat(sprintf("Topic proportions:\n"))
  print(x$topic.proportions)
  cat(sprintf("Topic representatives:\n"))
  print(round(x$topic.reps,digits = 3))
  return(invisible(x))
}

# Given a matrix of topic proportions, L, output topic proportion
# histograms as a k x n matrix, where k is the number of topics, and n
# is the number of bins in the histogram.
summarize_topic_proportions <- function (L) {
  k   <- ncol(L)
  out <- t(apply(L,2,
                 function (x) as.vector(table(cut(x,c(-1,0.1,0.5,0.9,1))))))
  if (is.null(rownames(out)))
    rownames(out) <- 1:k
  colnames(out) <- c("<0.1","0.1-0.5","0.5-0.9",">0.9")
  return(out)
}
  
# Given a matrix of topic proportions, L, return a k x k matrix, where
# k is the number of topics, in which each row is a sample (data
# matrix row). The ith row of the matrix contains the topic
# proportions for the sample "most representative" of the ith topic;
# that is, the row or sample with the largest proportion of counts
# drawn from the ith topic.
get_topic_representatives <- function (L) {
  n    <- nrow(L)
  k    <- ncol(L)
  rows <- apply(L,2,which.max)
  out  <- L[rows,]
  if (is.null(rownames(out)))
    rownames(out) <- rows
  if (is.null(colnames(out)))
    colnames(out) <- 1:k
  return(out)
}
