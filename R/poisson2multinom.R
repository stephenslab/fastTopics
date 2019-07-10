#' @rdname altsqp
#' 
#' @export
#' 
poisson2multinom <- function (fit) {

  # Input argument "fit" should be a list with elements "F" and "L".
  verify.fit(fit)
  F <- fit$F
  L <- fit$L

  # Verify and process input F.
  verify.matrix(fit$F)
  if (!is.matrix(F))
    F <- as.matrix(F)
  if (is.integer(F))
    storage.mode(F) <- "double"
  
  # Verify and process input L.
  verify.matrix(fit$L)
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
  return(list(F = F,L = L,s = s))
}
