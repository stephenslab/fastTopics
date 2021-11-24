#' @rdname select_loadings
#'
#' @title Extract or Re-order Data Rows in Poisson NMF or Multinomial Topic Model Fit
#'
#' @description This function can be used to extract estimates for a
#'   subset of the count data, or to re-order the rows of the loadings
#'   matrix.
#'
#' @param .data Poisson NMF or Multinomial Topic Model fit; that is,
#'   an object of class \dQuote{poisson_nmf_fit} or
#'   \dQuote{multinom_topic_model_fit}, such as an output from
#'   \code{\link{fit_poisson_nmf}} or \code{\link{fit_topic_model}}.
#'
#' @param loadings Indices (names or numbers) giving data rows to
#'   keep. If not specified, all rows are kept.
#' 
#' @param \dots Other arguments passed to the generic select function.
#' 
#' @return A Poisson NMF or multinomial topic model fit containing the
#'   selected data rows only.
#'
#' @importFrom dplyr select
#'
#' @aliases select
#' 
#' @method select poisson_nmf_fit
#'
#' @export
#' 
select.poisson_nmf_fit <- function (.data, loadings, ...)
  select_loadings(.data,loadings,...)

#' @rdname select_loadings
#' 
#' @method select multinom_topic_model_fit
#'
#' @export
#' 
select.multinom_topic_model_fit <- function (.data, loadings, ...)
  select_loadings(.data,loadings,...)

#' @rdname select_loadings
#'
#' @export
#'
select_loadings <- function (.data, loadings, ...) {
  if (!(inherits(.data,"poisson_nmf_fit") |
        inherits(.data,"multinom_topic_model_fit")))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\" or ",
         "\"multinom_topic_model_fit\"")
  verify.fit(.data)
  n <- nrow(.data$L)
  if (missing(loadings))
    loadings <- 1:n
  tryCatch({  
    .data$L  <- .data$L[loadings,,drop = FALSE]
    .data$Ln <- .data$Ln[loadings,,drop = FALSE]
    .data$Ly <- .data$Ly[loadings,,drop = FALSE]
    .data$s  <- .data$s[loadings]
  },error = function (e) stop("Invalid selection of loadings"))
  return(.data)
}
