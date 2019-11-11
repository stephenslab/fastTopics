# TO DO: Explain here briefly what this function does, and how to use
# it.
#
#' @keywords internal
#' 
#' @export
#'
lda <- function (X, F, L, alpha = rep(1,ncol(F)), numiter = 1000) {

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
  cat("iter --objective(ELBO)-- max.diff\n")
  for (iter in 1:numiter) {
    L0 <- L
    F0 <- F

    # E STEP
    # ------
    # Update the expected topic counts (N) and expected word counts (M).
    N <- matrix(0,n,k)
    M <- matrix(0,m,k)
    for (i in 1:n) {
      P     <- scale.cols(F,exp(digamma(L[i,])))
      P     <- P / rowSums(P)
      N[i,] <- X[i,] %*% P
      M     <- M + X[i,] * P
    }

    # M STEP
    # ------
    # Update the topic proportions (loadings).
    L <- alpha + N

    # Update the word probabilities (factors).
    F <- scale.cols(M + 1e-6)
    
    # Compute the variational lower bound at the current solution.
    value[iter] <- elbo.lda(X,F,L,alpha)
    cat(sprintf("%4d %+0.12e %0.2e\n",iter,value[iter],
                max(max(abs(L - L0)),max(abs(F - F0)))))
  }

  # Return the estimates of the topic proportions (L) and word
  # probabilities (F), and the value of the objective at each
  # iteration ("value").
  return(list(F = F,L = L,value = value))
}

# TO DO: Explain here what this function does, and how to use it.
elbo.lda <- function (X, F, L, alpha) {
  n <- nrow(X)
  f <- rep(0,n)
  for (i in 1:n) {
    L[i,] <- L[i,] * (sum(alpha) + sum(X[i,]))
    P     <- scale.cols(F,exp(digamma(L[i,])))
    P     <- P / rowSums(P)
    u     <- digamma(L[i,]) - digamma(sum(L[i,]))
    f[i]  <- (lgamma(sum(alpha)) - lgamma(sum(L[i,]))
              + sum(lgamma(L[i,])) - sum(lgamma(alpha))
              + sum((alpha - L[i,]) * u)
              + sum(X[i,] %*% (scale.cols(P,u) + P*log(F) - P*log(P))))
  }
  return(f)
}
