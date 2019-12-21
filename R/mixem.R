# Compute maximum-likelihood estimates of the mixture proportions in a
# mixture model by iterating the EM updates for a fixed number of
# iterations. This is mainly used for testing the C++ implementation.
mixem <- function (L, w, x0, numiter, e = 1e-15) {

  # Get the initial estimate of the solution.
  x <- x0

  # Iterate the E and M steps.
  for (i in 1:numiter)
    x <- mixem.update(L,w,x,e)

  # Return the updated mixture proportions.
  return(x)
}

# Perform a single EM update. This is mainly used for testing the C++
# implementation.
mixem.update <- function (L, w, x, e = 1e-15) {

  # Compute the posterior mixture assignment probabilities. This is
  # the "E step". 
  P <- scale.cols(L,x)
  P <- normalize.rows(P) + e
  P <- P / rowSums(P)

  # Update the mixture weights. This is the "M step".
  w <- w/sum(w)
  return(drop(w %*% P))
}
