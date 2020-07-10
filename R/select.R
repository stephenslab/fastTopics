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
#' @param \dots Other arguments passed to the generic select function.
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

