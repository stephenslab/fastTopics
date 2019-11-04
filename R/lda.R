# TO DO: Explain here briefly what this function does, and how to use
# it.
lda <- function (X, L, F, alpha = rep(1,ncol(F)), numiter = 1000) {

  # Get the number of rows (n) and columns (m) of X, and the number of
  # topics.
  n <- nrow(X)
  m <- ncol(X)
  k <- ncol(F)
    
  # This variable is used to keep track of the algorithm's progress;
  # it stores the value of the objective (the variational lower bound,
  # or "ELBO") at each iteration.
  value <- rep(0,numiter)
  
  # Iterate the E and M steps.
  for (iter in 1:numiter) {

    # E STEP
    # ------
    # Update the expected topic counts (N) and expected word counts (M).
    N <- matrix(0,n,k)
    M <- matrix(0,m,k)
    for (i in 1:n)
      for (j in 1:m) {
        pij   <- F[j,] * exp(digamma(L[i,]))
        pij   <- pij / sum(pij)
        N[i,] <- N[i,] + X[i,j]*pij
        M[j,] <- M[j,] + X[i,j]*pij
      }

    # M STEP
    # ------
    # Update the topic proportions (loadings).
    L <- alpha + N

    # Update the word probabilities (factors).
    F <- scale.cols(M)
    
    # Compute the variational lower bound at the current solution.
    # TO DO.
  }

  # Return the estimates of the topic proportions (L) and word
  # probabilities (F), and the value of the objective at each
  # iteration ("value").
  return(list(F = F,L = L,value = value))
}

# TO DO: Explain here what this function does, and how to use it.
elbo.lda <- function (X, F, L, alpha) {
  # TO DO.
}
