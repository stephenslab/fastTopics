#' @title Recover Binomial Topic Model Fit from Poisson NMF fit
#'
#' @examples
#' # See the vignette for an example.
#' 
#' @export
#'
poisson2binom <- function (fit) {

  # Check input argument "fit".
  if (inherits(fit,"binom_topic_model_fit"))
    return(fit)
  if (!inherits(fit,"poisson_nmf_fit"))
    stop("Input argument \"fit\" should be an object of class ",
         "\"poisson_nmf_fit\"")
  verify.fit(fit)
  if (ncol(fit$F) < 2 | ncol(fit$L) < 2)
    stop("Input matrices \"fit$F\" and \"fit$L\" should have 2 or more",
         "columns")

  # TO DO.
}
