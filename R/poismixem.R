# Compute a maximum-likelihood estimate (MLE) of the mixture weights
# in a Poisson mixture model by iterating the multinomial mixture
# model EM updates for a fixed number of iterations. This is mainly
# used for testing the C++ implementation. See the comments attached
# to the "poismixem" C++ function for an explanation of the inputs.
poismixem <- function (L, w, x0, numiter) {
  x   <- x0
  s   <- sum(w)
  out <- p2m(L,x)
  x   <- mixem(out$L,w,out$x,numiter)
  return(m2p(L,x,s))
}

# Recover the mixture weights of the multinomial mixture model from
# the mixture weights of the Poisson mixture model. The "scale factor"
# (s) is also returned. This is mainly used for testing the C++
# implementation.
p2m <- function (L, x) {
  s <- sum(L %*% x)
  u <- colSums(L)
  return(list(L = normalize.cols(L),x = x*u/s,s = s))
}

# Recover the mixture weights of the Poisson mixture model from the
# mixture weights of the multinomial mixture model, plus the "scale
# factor" (s). This is mainly used for testing the C++
# implementation.
m2p <- function (L, x, s)
  s*x/colSums(L)
