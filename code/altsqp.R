# Compute maximum-likelihood estimates for the Poisson topic model;
# equivalently, find a non-negative matrix factorization X = L*F' that
# optimizes the beta divergence objective.
altsqp <- function (X, F, L, numiter = 100, nem = 1, nsqp = 4, tol = 1e-10,
                    zero.threshold = 0, zero.searchdir = 1e-15,
                    suffdecr = 0.01, stepsizereduce = 0.75,
                    minstepsize = 1e-10, e = 1e-15, verbose = TRUE) {
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
        if (nem > 0)
          L[i,] <- altsqp.update.em(F,X[i,],L[i,],nem,e)
        if (nsqp > 0)
          L[i,] <- altsqp.update.sqp(F,X[i,],L[i,],nsqp,e,tol,zero.threshold,
                                     zero.searchdir,suffdecr,stepsizereduce,
                                     minstepsize)
      }

      # Update the factors ("basis vectors").
      for (j in 1:m) {
        if (nem > 0)
          F[j,] <- altsqp.update.em(L,X[,j],F[j,],nem,e)
        if (nsqp > 0)
          F[j,] <- altsqp.update.sqp(L,X[,j],F[j,],nsqp,e,tol,zero.threshold,
                                     zero.searchdir,suffdecr,stepsizereduce,
                                     minstepsize)
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

# Run a fixed number of EM updates for the alternating SQP method.
altsqp.update.em <- function (B, w, y, numiter, e) {
  ws  <- sum(w)
  bs  <- colSums(B)  
  out <- mixem(scale.cols(B,ws/bs),w/ws,y*bs/ws,numiter,e)
  return(out$x*ws/bs)
}

# Run a fixed number of SQP updates for the alternating SQP method.
altsqp.update.sqp <- function (B, w, y, numiter, e, tol, zero.threshold,
                               zero.searchdir, suffdecr, stepsizereduce,
                               minstepsize) {
  ws  <- sum(w)
  bs  <- colSums(B)  
  out <- mixsqp(scale.cols(B,ws/bs),w/ws,y*bs/ws,numiter,e,tol,
                zero.threshold,zero.searchdir,suffdecr,stepsizereduce,
                minstepsize)
  return(out$x*ws/bs)
}
