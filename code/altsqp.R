# TO DO: Explain here what this function does, and how to use it.
altsqp <- function (X, F, L, numiter = 100, nem = 1, e = 1e-15,
                    verbose = TRUE) {
  n <- nrow(X)
  m <- ncol(X)
  progress <- data.frame(iter      = 1:numiter,
                         objective = 0,
                         max.diff  = 0,
                         timing    = 0)
  
  # Iteratively apply the EM And SQP updates.
  if (verbose)
    cat("iter         objective max.diff\n")
  for (iter in 1:numiter) {

    # Save the current estimates of the factors and loadings.
    F0 <- F
    L0 <- L

    timing <- system.time({

      # Update the loadings ("activations").
      for (i in 1:n) {
        fi    <- cost.poismix(F,X[i,],L[i,],e)
        out   <- fitpoismix.update(F,X[i,],L[i,],fi,e = e,
                                   qp.solver = "activeset")
        L[i,] <- out$x
      }
      
      # Update the factors ("basis vectors").
      for (j in 1:m) {
        fj    <- cost.poismix(L,X[,j],F[j,],e)
        out   <- fitpoismix.update(L,X[,j],F[j,],fj,e = e,
                                   qp.solver = "activeset")
        F[j,] <- out$x
      }
    })

    # Compute the value of the objective (cost) function at the
    # current estimates of the factors and loadings.
    f <- cost(X,tcrossprod(L,F),e)
    d <- max(max(abs(F - F0)),max(abs(L - L0)))
    progress[iter,"objective"] <- f
    progress[iter,"max.diff"]  <- d
    progress[iter,"timing"]    <- timing["elapsed"]
    if (verbose)
      cat(sprintf("%4d %+0.10e %0.2e\n",iter,f,d))
  }
  
  return(list(F = F,L = L,value = f,progress = progress))
}
