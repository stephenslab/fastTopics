#' @title Summarize Poisson NMF or Multinomial Topic Model Fit
#'
#' @description \code{summary} method for the \dQuote{poisson_nmf_fit}
#'   and \dQuote{multinom_topic_model_fit} classes.
#'
#' @param object An object of class \dQuote{poisson_nmf_fit} or
#'   \dQuote{multinom_topic_model_fit}. The former is usually the result
#'   of calling \code{\link{fit_poisson_nmf}}; the latter is usually the
#'   result of calling \code{\link{fit_topic_model}} or
#'   \code{\link{poisson2multinom}}.
#'
#' @method summary poisson_nmf_fit
#'
#' @return The functions \code{summary.poisson_nmf_fit} and
#' \code{summary.multinom_topic_model_fit} compute and return a list
#' of statistics summarizing the model fit. The returned list
#' includes some or all of the following elements:
#'
#' \item{n}{The number of rows in the counts matrix, typically the
#'   number of samples.}
#'
#' \item{m}{The number of columns in the counts matrix, typically the
#'   number of observed counts per sample.}
#'
#' \item{k}{The rank of the Poisson NMF or the number of topics.}
#'
#' \item{s}{A vector of length n giving the "size factor" estimates;
#'   these estimates should be equal, or close to, the total counts in
#'   each row of the counts matrix.}
#'
#' \item{numiter}{The number of loadings and/or factor updates
#'   performed.}
#'
#' \item{loglik}{The Poisson NMF log-likelihood.}
#'
#' \item{loglik.multinom}{The multinomial topic model log-likelihood.}
#'
#' \item{dev}{The Poisson NMF deviance.}
#'
#' \item{res}{The maximum residual of the Karush-Kuhn-Tucker (KKT)
#'   first-order optimality conditions. This can be used to assess
#'   convergence of the updates to a (local) solution.}
#'
#' \item{mixprops}{Matrix giving a high-level summary of the
#'   mixture proportions, in which rows correspond to topics, and
#'   columns are ranges of mixture proportionss.}
#'
#' \item{topic.reps}{A matrix in which the ith row gives the mixture
#'   proportions for the sample "most representative" of topic i; by
#'   "most representative", we mean the row (or sample) with the highest
#'   proportion of counts drawn from the topic i.}
#' 
#' @export
#' 
summary.poisson_nmf_fit <- function (object, ...) {
  out <- summary.multinom_topic_model_fit(poisson2multinom(object),...)
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
  verify.fit(object)
  numiter <- nrow(object$progress)
  out <- list(n          = nrow(object$L),
              m          = nrow(object$F),
              k          = ncol(object$F),
              s          = object$s,
              numiter    = numiter,
              loglik     = object$progress[numiter,"loglik"],
              loglik.multinom = object$progress[numiter,"loglik.multinom"],
              dev        = object$progress[numiter,"dev"],
              res        = object$progress[numiter,"res"],
              mixprops   = summarize_mixprops(object$L),
              topic.reps = get_topic_reps(object$L))
  class(out) <- c("summary.multinom_topic_model_fit","list")
  return(out)
}

#' @rdname summary.poisson_nmf_fit
#'
#' @param x An object of class \dQuote{summary.poisson_nmf_fit},
#'   usually a result of a call to \code{summary.poisson_nmf_fit}.
#'
#' @param show.mixprops If \code{TRUE}, print a summary of the mixture
#'   proportions.
#'
#' @param show.topic.reps If \code{TRUE}, print a summary of the topic
#'   representatives.
#' 
#' @param \dots Additional arguments passed to the generic \code{summary}
#'   or \code{print.summary} method.
#'
#' @method print summary.poisson_nmf_fit
#' 
#' @export
#' 
print.summary.poisson_nmf_fit <- function (x, show.mixprops = FALSE,
                                           show.topic.reps = FALSE, ...) {
  print.summary.multinom_topic_model_fit(x,FALSE,show.mixprops,
                                         show.topic.reps,...)
  return(invisible(x))
}

#' @rdname summary.poisson_nmf_fit
#'
#' @param show.size.factors If \code{TRUE}, print a summary of the
#'   size factors.
#'
#' @method print summary.multinom_topic_model_fit
#'
#' @importFrom stats quantile
#' 
#' @export
#'
print.summary.multinom_topic_model_fit <-
  function (x, show.size.factors = FALSE, show.mixprops = FALSE,
            show.topic.reps = FALSE, ...) {
  k <- x$k
  cat("Model overview:\n")
  cat(sprintf("  Number of data rows, n: %d\n",x$n))
  cat(sprintf("  Number of data cols, m: %d\n",x$m))
  cat(sprintf("  Rank/number of topics, k: %d\n",x$k))
  cat(sprintf("Evaluation of model fit (%d updates performed):\n",
              x$numiter))
  cat(sprintf("  Poisson NMF log-likelihood: %+0.12e\n",x$loglik))
  cat(sprintf("  Multinomial topic model log-likelihood: %0.12e\n",
              x$loglik.multinom))
  cat(sprintf("  Poisson NMF deviance: %+0.12e\n",x$dev))
  cat(sprintf("  Max KKT residual: %+0.6e\n",x$res))
  if (!(show.size.factors | show.mixprops | show.topic.reps))
    message("Set show.size.factors = TRUE, show.mixprops = TRUE and/or ",
            "show.topic.reps = TRUE in print(...) for more information")
  if (show.size.factors & !is.null(x$s)) {
    cat(sprintf("Size factors:\n"))
    q        <- quantile(x$s)
    names(q) <- c("Min","1Q","Median","3Q","Max")
    print(q)
  }
  if (show.mixprops) {
    cat(sprintf("Mixture proportions:\n"))
    print(x$mixprops)
  }
  if (show.topic.reps) {
    cat(sprintf("Topic representatives:\n"))
    print(round(x$topic.reps,digits = 3))
  }
  return(invisible(x))
}

#' @title Summarize and Compare Model Fits
#'
#' @description Create a table summarizing the results of fitting one
#'   or more Poisson non-negative matrix factorizations or multinomial
#'   topic models.
#'
#' @param fits An object of class \code{"poisson_nmf_fit"} or
#'   \code{"multinom_topic_model_fit"}, or a non-empty, named list in
#'   which all list elements are Poisson NMF model fits or all
#'   multinomial topic model fits.
#'
#' @return A data frame with one row per element of \code{fits}, and
#' with the following columns:
#'
#' \item{k}{The rank of the matrix factorization.}
#'
#' \item{loglik}{The log-likelihood (either Poisson NMF or multinomial topic
#'   model likelihood) achieved at the last model fitting update.}
#'
#' \item{dev}{For Poisson NMF model fits only, the deviance achieved
#'   at the last model fitting update.}
#'
#' \item{res}{The maximum residual of the Karush-Kuhn-Tucker (KKT)
#'   system achieved at the last model fitting update; small values
#'   indicate that the solution is close to a local maximum, or
#'   stationary point, of the likelihood.}
#'
#' \item{loglik.diff}{The improvement in the log-likelihood relative
#'   to the model fit with the smallest log-likelihood.}
#'
#' \item{dev.diff}{The improvement in the deviance relative to the
#'   model fit with the largest deviance.}
#' 
#' \item{nonzeros.f}{The rate of nonzeros in the factors matrix, as
#'   determined by \code{control$zero.threshold}.}
#'
#' \item{nonzeros.l}{The rate of nonzeros in the loadings matrix, as
#'   determined by \code{control$zero.threshold}.}
#'
#' \item{numiter}{The number of loadings and/or factor updates
#'   performed.}
#' 
#' \item{runtime}{The total runtime (in s) of the model fitting
#'   updates.}
#' 
#' @seealso \code{\link{fit_poisson_nmf}},
#'   \code{\link{fit_topic_model}}
#' 
#' @export
#' 
compare_fits <- function (fits) {

  # Check and process input "fits". It should either be an object of
  # class poisson_nmf_fit or multinom_topic_model_fit, or a list of
  # fit objects.
  if (inherits(fits,"poisson_nmf_fit")) {
    verify.fit(fits)
    fit.name    <- deparse(substitute(fits))          
    fits        <- list(fits)
    names(fits) <- fit.name
  } else {
    msg <- paste("Input argument \"fits\" should either be an object of",
                 "class \"poisson_nmf_fit\", or a non-empty, named list in",
                 "which all list elements are objects of class",
                 "\"poisson_nmf_fit\" or all objects of class",
                 "\"multinom_topic_model_fit\"")
    if (!(is.list(fits) & !is.null(names(fits)) & length(fits) > 0))
      stop(msg)
    if (!(all(sapply(fits,function(x)inherits(x,"poisson_nmf_fit"))) |
          all(sapply(fits,function(x)inherits(x,"multinom_topic_model_fit")))))
      stop(msg)
    if (!all(nchar(names(fits)) > 0))
      stop(msg)
    sapply(fits,verify.fit)
  }

  # Initialize the data structure.
  n   <- length(fits)
  out <- data.frame(k           = rep(0,n),
                    loglik      = rep(0,n),
                    dev         = rep(0,n),
                    loglik.diff = rep(0,n),
                    dev.diff    = rep(0,n),
                    res         = rep(0,n),
                    nonzeros.f  = rep(0,n),
                    nonzeros.l  = rep(0,n),
                    numiter     = rep(0,n),
                    runtime     = rep(0,n),
                    stringsAsFactors = FALSE)
  rownames(out) <- names(fits)
  
  # Collect the information from each fit.
  for (i in 1:n) {
    fit <- fits[[i]]
    x   <- tail(fit$progress,n = 1)
    out[i,"k"]          <- ncol(fit$F)
    out[i,"loglik"]     <- ifelse(inherits(fit,"poisson_nmf_fit"),x$loglik,
                                  x$loglik.multinom)
    out[i,"dev"]        <- ifelse(inherits(fit,"poisson_nmf_fit"),x$dev,NA)
    out[i,"res"]        <- x$res
    out[i,"nonzeros.f"] <- x$nonzeros.f
    out[i,"nonzeros.l"] <- x$nonzeros.l
    out[i,"numiter"]    <- nrow(fit$progress)
    out[i,"runtime"]    <- sum(fit$progress$timing)
  }

  # Calculate the "loglik.diff" and "dev.diff" entries.
  out$loglik.diff <- out$loglik - min(out$loglik)
  out$dev.diff    <- max(out$dev) - out$dev
  
  # Output the data frame.
  return(out)
}

# Given L, the matrix of mixture proportions, output mixture
# proportion histograms as a k x n matrix, where k is the number of
# topics, and n is the number of bins in the histogram.
summarize_mixprops <- function (L) {
  k   <- ncol(L)
  out <- t(apply(L,2,
                 function (x) as.vector(table(cut(x,c(-1,0.1,0.5,0.9,1))))))
  if (is.null(rownames(out)))
    rownames(out) <- 1:k
  colnames(out) <- c("<0.1","0.1-0.5","0.5-0.9",">0.9")
  return(out)
}
  
# Given a matrix of mixture proportions, L, return a k x k matrix,
# where k is the number of topics, in which each row is a sample (data
# matrix row). The ith row of the matrix contains the mixture
# proportions for the sample "most representative" of the ith topic;
# that is, the row or sample with the largest proportion of counts
# drawn from the ith topic.
get_topic_reps <- function (L) {
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
