# This function decomposes the input matrix X = A*B by nonnegative
# matrix factorization (NMF) based on the beta-divergence criterion
# (negative Poisson log-likelihood) and multiplicative update
# rules. All entries of initial estimates A and B should be
# positive. This is adapted from the MATLAB code by D. Kitamura
# (http://d-kitamura.net).
betanmf <- function (X, A, B, numiter = 1000, e = 1e-8, verbose = TRUE) {
  n <- nrow(X)
  p <- ncol(X)
  E <- matrix(1,n,p)
  progress <- data.frame(iter = 1:numiter,objective = 0,max.diff  = 0,
                         timing = 0)
      
  # Repeat until we reach the number of requested iterations.
  if (verbose)
    cat("iter         objective max.diff\n")
  for (i in 1:numiter) {

    # Save the current estimates of the factors and loadings.
    A0 <- A
    B0 <- B

    timing <- system.time({
    
      # Update the loadings ("activations").
      A <- A * (((X / (A %*% B)) %*% t(B)) / (E %*% t(B)))
      A <- pmax(A,e)
    
      # Update the factors ("basis vectors").
      B <- B * ((t(A) %*% (X / (A %*% B))) / (t(A) %*% E))
      B <- pmax(B,e)
    })

    # Compute the value of the objective (cost) function at the
    # current estimates of the factors and loadings.
    f <- cost(X,A %*% B,e)
    d <- max(max(abs(A - A0)),max(abs(B - B0)))
    progress[i,"objective"] <- f
    progress[i,"max.diff"]  <- d
    progress[i,"timing"]    <- timing["elapsed"]
    if (verbose)
      cat(sprintf("%4d %0.10e %0.2e\n",i,f,d))
  }

  return(list(A = A,B = B,value = f,progress = progress))
}
