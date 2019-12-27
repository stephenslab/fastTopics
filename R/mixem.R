# Compute a maximum-likelihood estimate (MLE) of the mixture
# proportions in the multinomial mixture model by iterating the EM
# updates for a fixed number of iterations. This is mainly used for
# testing the C++ implementation. See the comments attached to the
# "mixem" C++ function for an explanation of the inputs.
mixem <- function (L, w, x0, numiter) {
  L <- normalize.cols(L)
  x <- x0/sum(x0)
  for (i in 1:numiter)
    x <- mixem.update(L,w,x)
  return(x)
}

# Perform a single EM update for the multinomial mixture model. This
# is mainly used for testing the C++ implementation.
mixem.update <- function (L, w, x) {
  e <- 1e-15
  
  # Compute the posterior mixture assignment probabilities. A small
  # number is added to the posterior probabilities to prevent any
  # divisions by zero. This is the "E step".
  P <- scale.cols(L,x)
  P <- normalize.rows.by.max(P) + e
  P <- normalize.rows(P)

  # Update the mixture weights. This is the "M step".
  w <- w/sum(w)
  return(drop(w %*% P))
}

# Find the MLE for the special case when only one of the counts is
# nonzero.
mixture.one.nonzero <- function (L, w) {
  x <- rep(0,ncol(L))
  x[which.max(w %*% L)] <- 1
  return(x)
}
