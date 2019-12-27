# Compute a maximum-likelihood estimate (MLE) of the mixture weights
# in a Poisson mixture model by iterating the multinomial mixture
# model EM updates for a fixed number of iterations. This is mainly
# used for testing the C++ implementation. See the comments attached
# to the "poismixem" C++ function for an explanation of the inputs.
poismixem <- function (L, w, x0, numiter) {
  x <- x0

  # Recover the mixture weights of the multinomial mixture model from
  # the mixture weights of the Poisson mixture model. Here, s is the
  # "scale factor".
  s <- sum(L %*% x)
  u <- colSums(L)
  L <- normalize.cols(L)
  x <- x*u/s

  # Perform one or more EM updates for the multinomial mixture model.
  x <- mixem(L,w,x,numiter)

  # Recover the mixture weights of the Poisson mixture model from the
  # mixture weights of the multinomial mixture model.
  s <- sum(w)
  return(s*x/u)
}

