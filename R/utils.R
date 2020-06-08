#' @title Extract or Re-order Multinomial Topic Model Loadings
#'
#' @description This function can be used to extract estimates for a
#'   subset of the count data, or to re-order the rows of the loadings
#'   matrix.
#'
#' @param .data Multinomial Topic Model fit; that is, an object of
#'   class \dQuote{multinom_topic_model_fit}, such as an output from
#'   \code{poisson2multinom}.
#'
#' @param loadings Loadings indices (names or numbers), corresponding
#'   to rows of the counts matrix. If not specified, all loadings are
#'   returned.
#' 
#' @param ... Other arguments passed to the generic select function.
#' 
#' @return A multinomial topic model fit containing the selected
#'   loadings only.
#'
#' @seealso \code{\link{fit_poisson_nmf}},
#'   \code{\link{poisson2multinom}}
#'
#' @importFrom dplyr select
#'
#' @method select multinom_topic_model_fit
#'  
#' @export
#' 
select.multinom_topic_model_fit <- function (.data, loadings, ...) {
  n <- nrow(.data$L)
  if (missing(loadings))
    loadings <- 1:n
  tryCatch({  
    .data$L  <- .data$L[loadings,,drop = FALSE]
    .data$Ln <- .data$Ln[loadings,,drop = FALSE]
    .data$Ly <- .data$Ly[loadings,,drop = FALSE]
  },error = function (e) stop("Invalid selection of loadings"))
  return(.data)
}

#' @title Summarize and Compare Poisson NMF Model Fits
#'
#' @description Create a table summarizing the results of fitting one
#'    or more Poisson NMF models.
#'
#' @param fits An object of class \code{"poisson_nmf_fit"}, or a
#'   non-empty, named list in which each list element is an object of
#'   class \code{"poisson_nmf_fit"}.
#'
#' @return A data frame with one row per element of \code{fits}, and
#' with the following columns:
#'
#' \item{k}{The rank of the matrix factorization.}
#'
#' \item{loglik}{The log-likelihood achieved at the last model fitting
#'   update.}
#'
#' \item{dev}{The deviance achieved at the last model fitting update.}
#'
#' \item{res}{The maximum residual of the Karush-Kuhn-Tucker (KKT)
#'   system achieved at the last model fitting update; small values
#'   indicate that the solution is close to a local maximum, or
#'   stationary point, of the likelihood.}
#'
#' \item{nonzeros.f}{The rate of nonzeros in the factors matrix, as
#'   determined by \code{control$zero.threshold}.}
#'
#' \item{nonzeros.l}{The rate of nonzeros in the loadings matrix, as
#'   determined by \code{control$zero.threshold}.}
#' 
#' \item{runtime}{The total runtime (in s) of the model fitting
#'   updates.}
#' 
#' @seealso \code{\link{fit_poisson_nmf}}
#' 
#' @export
#' 
compare_poisson_nmf_fits <- function (fits) {

  # CHECK & PROCESS INPUTS
  # ----------------------
  # Check and process input "fits". It should either be an object of
  # class poisson_nmf_fit, or a list of poisson_nmf_fit objects.
  if (inherits(fits,"poisson_nmf_fit")) {
    fit.name    <- deparse(substitute(fits))          
    fits        <- list(fits)
    names(fits) <- fit.name
  } else {
    msg <- paste("Input argument \"fits\" should either be an object of",
                 "class \"poisson_nmf_fit\", or a non-empty, named list in",
                 "which each list element is an object of class",
                 "\"poisson_nmf_fit\"")
    if (!(is.list(fits) & !is.null(names(fits)) & length(fits) > 0))
      stop(msg)
    if (!all(sapply(fits,function (x) inherits(x,"poisson_nmf_fit"))))
      stop(msg)
    if (!all(nchar(names(fits)) > 0))
      stop(msg)
  }

  # CREATE SUMMARY TABLE
  # --------------------
  # Initialize the data structure.
  n   <- length(fits)
  out <- data.frame(k       = rep(0,n),
                    loglik  = rep(0,n),
                    dev     = rep(0,n),
                    res     = rep(0,n),
                    runtime = rep(0,n),
                    stringsAsFactors = FALSE)
  rownames(out) <- names(fits)
  
  # Collect the information from each fit.
  for (i in 1:n) {
    fit <- fits[[i]]
    x   <- tail(fit$progress,n = 1)
    out[i,"k"]          <- ncol(fit$F)
    out[i,"loglik"]     <- x$loglik
    out[i,"dev"]        <- x$dev
    out[i,"res"]        <- x$res
    out[i,"nonzeros.f"] <- x$nonzeros.f
    out[i,"nonzeros.l"] <- x$nonzeros.l
    out[i,"runtime"]    <- sum(fit$progress$timing)
  }

  # Output the data frame.
  return(out)
}
