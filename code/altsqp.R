# TO DO: Explain here what this function does, and how to use it.
altsqp <- function (X, F, L, numiter = 1000, e = 1e-8, verbose = TRUE) {
  n <- nrow(X)
  m <- ncol(X)
  progress <- data.frame(iter = 1:numiter,objective = 0,
                         max.diff = 0,timing = 0)
  
  # Repeat until we reach the number of requested iterations.
  if (verbose)
    cat("iter         objective max.diff\n")
  for (iter in 1:numiter) {

    # Save the current estimates of the factors and loadings.
    F0 <- F
    L0 <- L

    timing <- system.time({

      # Update the loadings ("activations").
      L <- altsqp.update.loadings(X,F,L,e)
      # rows <- splitIndices(n,nc)
      # L <- mclapply(rows,function(i) altsqp.update.loadings(X[i,],F,L[i,],e))
      # L <- do.call(rbind,L)
      # L[unlist(rows),] <- L

      # Update the factors ("basis vectors").
      F <- altsqp.update.factors(X,F,L,e)
      # cols <- splitIndices(m,nc)
      # F <- mclapply(cols,function(j) altsqp.update.factors(X[,j],F[j,],L,e))
      # F <- do.call(rbind,F)
      # F[unlist(cols),] <- F
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

# This implements the loadings update step in the alternating SQP method.
altsqp.update.loadings <- function (X, F, L, e) {
  n <- nrow(X)
  for (i in 1:n) {
    fi    <- cost.poismix(F,X[i,],L[i,],e)
    out   <- fitpoismix.update(F,X[i,],L[i,],fi,e = e,qp.solver = "activeset")
    L[i,] <- out$x
  }
  return(L)
}

# This implements the factors update step in the alternating SQP method.
altsqp.update.factors <- function (X, F, L, e) {
  m <- ncol(X)
  for (j in 1:m) {
    fj    <- cost.poismix(L,X[,j],F[j,],e)
    out   <- fitpoismix.update(L,X[,j],F[j,],fj,e = e,qp.solver = "activeset")
    F[j,] <- out$x
  }
  return(F)
}

