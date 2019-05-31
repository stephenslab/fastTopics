# This function decomposes the input matrix X = A*B by nonnegative
# matrix factorization (NMF) based on the beta-divergence criterion
# (negative Poisson log-likelihood) and multiplicative update
# rules. All entries of initial estimates A and B should be
# positive. This is adapted from the MATLAB code by D. Kitamura
# (http://d-kitamura.net).
betanmf <- function (X, A, B, numiter = 1000, e = 1e-15, verbose = TRUE) {
  n <- nrow(X)
  m <- ncol(X)
  E <- matrix(1,n,m)
  value <- data.frame(iter      = 1:numiter,
                         objective = 0,
                         max.diff  = 0,
                         timing    = 0)
      
  # Iteratively apply the multiplicative updates.
  if (verbose)
    cat("        objective max.diff time(s)\n")
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
    value[i] <- cost(X,A %*% B,e)
    if (verbose) {
      d <- max(max(abs(A/rowMeans(A) - A0/rowMeans(A0))),
               max(abs(scale.cols(B,1/colMeans(B)) -
                       scale.cols(B0,1/rowMeans(B0)))))
      cat(sprintf("%0.10e %0.2e %07.2f\n",value[i],d,timing["elapsed"]))
    }
  }

  return(list(A = A,B = B,value = value))
}
