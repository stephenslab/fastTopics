#' @title Title Goes Here
#'
#' @description Describe function here.
#'
#' @param fit Describe input argument "fit" here.
#'
#' @param dims Describe input argument "dims" here.
#'
#' @param pca Describe input argument "pca" here.
#'
#' @param normalize Describe input argument "normalize" here.
#'
#' @param verbose Describe input argument "verbose" here.
#' @param ... Describe other input arguments here.
#'
#' @return Describe the return value here.
#' 
#' @importFrom Rtsne Rtsne
#' 
#' @export
#' 
tsne_from_topics <- function (fit, dims = 2, pca = FALSE, normalize = FALSE,
                              verbose = FALSE, ...) {
  if (!inherits(fit,"poisson_nmf_fit"))
    stop("Input \"fit\" should be an object of class \"poisson_nmf_fit\"")
  out <- Rtsne(fit$L,dims = dims,pca = pca,normalize = normalize,
               verbose = verbose,...)
  Y   <- out$Y
  rownames(Y) <- rownames(fit$L)
  colnames(Y) <- paste("d",1:dimes)
  return(Y)
}
