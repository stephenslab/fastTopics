# Compute a maximum-likelihood estimate (MLE) of the mixture
# proportions in the multinomial mixture model by iterating the EM
# updates for a fixed number of iterations. This is mainly used for
# testing the C++ implementation. See the comments attached to the
# "mixem" C++ function for an explanation of the inputs.
mixem <- function (L, w, x0, numiter) {
  L1 <- normalize.cols(L)
  x  <- x0
  for (i in 1:numiter)
    x <- mixem.update(L1,w,x)
  return(x)
}

# Perform a single EM update for the multinomial mixture model. This
# is mainly used for testing the C++ implementation.
mixem.update <- function (L1, w, x) {
  e <- 1e-15
  x <- x/sum(x)
  w <- w/sum(w)
  
  # Compute the posterior mixture assignment probabilities. A small
  # number is added to the posterior probabilities to prevent any
  # divisions by zero. This is the "E step".
  P <- scale.cols(L1,x)
  P <- normalize.rows.by.max(P) + e
  P <- normalize.rows(P)

  # Update the mixture weights. This is the "M step".
  return(drop(w %*% P))
}

# Find the MLE for the special case when only one of the counts is
# positive.
mixture.one.nonzero <- function (L, w) {
  L1   <- normalize.cols(L)
  m    <- ncol(L)
  x    <- rep(0,m)
  j    <- which.max(w %*% L1)
  x[j] <- 1
  return(x)
}
