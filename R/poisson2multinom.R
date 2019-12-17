#' @title Add title here.
#'
#' @description Add description here.
#'
#' @param fit Add description of "fit" here.
#' 
#' @export
#' 
poisson2multinom <- function (fit) {

  # Input argument "fit" should be a list with elements "F" and "L".
  verify.fit(fit)
  F <- fit$F
  L <- fit$L

  # Verify and process input matrix F. Each column of F should have at
  # least one positive entry.
  verify.matrix(fit$F)
  if (any(colSums(F) <= 0))
    stop("Each column of \"fit$F\" should have at least one positive entry")
  if (!is.matrix(F))
    F <- as.matrix(F)
  if (is.integer(F))
    storage.mode(F) <- "double"
  
  # Verify and process input matrix L. Each row of L should have at
  # least one positive entry.
  verify.matrix(fit$L)
  if (any(rowSums(L) <= 0))
    stop("Each row of \"fit$L\" should have at least one positive entry")
  if (!is.matrix(L))
    L <- as.matrix(L)
  if (is.integer(L))
    storage.mode(L) <- "double"
  
  # Check that matrices F and L are compatible.
  if (ncol(L) != ncol(F))
    stop("Dimensions of \"fit$F\" and \"fit$L\" do not agree")
  
  # Convert the parameters---the factors F and the loadings L---for
  # the Poisson model to the factors and loadings for the multinomial
  # model. The return value "s" gives the Poisson rates for generating
  # the "document" sizes.
  L <- t(t(L) * colSums(F))
  s <- rowSums(L)
  L <- L / s
  F <- scale.cols(F)

  # Update the "fit" object, and return it.
  fit$F <- F
  fit$L <- L
  fit$s <- s
  return(fit)
}
