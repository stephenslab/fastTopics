# Compute maximum-likelihood estimates of the mixture proportions in a
# mixture model by iterating the EM updates for a fixed number of
# iterations.
mixem <- function (L, w, x0, numiter, e = 1e-15) {

  # Get the initial estimate of the solution.
  x <- x0

  # Iterate the E and M steps.
  n <- nrow(L)
  e <- rep(e,n)
  for (i in 1:numiter)
    x <- mixem.update(L,w,x,e)

  # Return (1) the estimate of the solution and (2) the value of the
  # objective at this estimate.
  return(list(x = x,value = mixobjective(L,w,x,e)))
}

# Perform a single EM update.
mixem.update <- function (L, w, x, e) {

  # Compute the n x m matrix of posterior mixture assignment
  # probabilities (L is an n x m matrix). This is the "E step".
  P <- scale.cols(L,x) + e
  P <- P / rowSums(P)

  # Update the mixture weights. This is the "M step".
  return(drop(w %*% P))
}
