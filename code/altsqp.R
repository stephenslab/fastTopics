# Compute maximum-likelihood estimates for the Poisson topic model;
# equivalently, find a non-negative matrix factorization X = L*F' that
# optimizes the beta divergence objective.
altsqp <- function (X, F, L, numiter = 100, nem = 1, nsqp = 4, nc = 1,
                    tol = 1e-10, zero.threshold = 0, zero.searchdir = 1e-15,
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
      if (nc == 1)
        L <- altsqp.update.loadings(X,F,L,nem,nsqp,e,tol,zero.threshold,
                                  zero.searchdir,suffdecr,stepsizereduce,
                                  minstepsize)
      else {
        rows <- splitIndices(n,nc)
        L    <- mclapply(rows,function (i)
                  altsqp.update.loadings(X[i,],F,L[i,],nem,nsqp,e,tol,
                                         zero.threshold,zero.searchdir,
                                         suffdecr,stepsizereduce,minstepsize))
        L    <- do.call(rbind,L)
        L[unlist(rows),] <- L
      }
          
      # Update the factors ("basis vectors").
      if (nc == 1)
        F <- altsqp.update.factors(X,F,L,nem,nsqp,e,tol,zero.threshold,
                                   zero.searchdir,suffdecr,stepsizereduce,
                                   minstepsize)
      else {
        cols <- splitIndices(m,nc)
        F    <- mclapply(cols,function (j)
                  altsqp.update.factors(X[,j],F[j,],L,nem,nsqp,e,tol,
                                        zero.threshold,zero.searchdir,
                                        suffdecr,stepsizereduce,minstepsize))
        F    <- do.call(rbind,F)
        F[unlist(cols),] <- F
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

# Update all the loadings with the factors remaining fixed.
altsqp.update.loadings <- function (X, F, L, nem, nsqp, e, tol, zero.threshold,
                                    zero.searchdir, suffdecr, stepsizereduce,
                                    minstepsize) {
  n <- nrow(X)
  for (i in 1:n) {
    if (nem > 0)
      L[i,] <- altsqp.update.em(F,X[i,],L[i,],nem,e)
    if (nsqp > 0)
      L[i,] <- altsqp.update.sqp(F,X[i,],L[i,],nsqp,e,tol,zero.threshold,
                                 zero.searchdir,suffdecr,stepsizereduce,
                                 minstepsize)
  }
  return(L)  
}

# Update all the factors with the loadings remaining fixed.
altsqp.update.factors <- function (X, F, L, nem, nsqp, e, tol, zero.threshold,
                                   zero.searchdir, suffdecr, stepsizereduce,
                                   minstepsize) {
  m <- ncol(X)
  for (j in 1:m) {
    if (nem > 0)
      F[j,] <- altsqp.update.em(L,X[,j],F[j,],nem,e)
    if (nsqp > 0)
      F[j,] <- altsqp.update.sqp(L,X[,j],F[j,],nsqp,e,tol,zero.threshold,
                                 zero.searchdir,suffdecr,stepsizereduce,
                                 minstepsize)
  }
  return(F)
}

# Run a fixed number of EM updates for the alternating SQP method.
altsqp.update.em <- function (B, w, y, numiter, e) {

  # Remove any counts that are exactly zero.
  ws <- sum(w)
  bs <- colSums(B)
  i  <- which(w > 0)
  w  <- w[i]
  B  <- B[i,]

  # Run the EM updates for the modified problem.
  out <- mixem(scale.cols(B,ws/bs),w/ws,y*bs/ws,numiter,e)

  # Recover the solution to the unmodified problem.
  return(out$x*ws/bs)
}

# Run a fixed number of SQP updates for the alternating SQP method.
altsqp.update.sqp <- function (B, w, y, numiter, e, tol, zero.threshold,
                               zero.searchdir, suffdecr, stepsizereduce,
                               minstepsize) {
    
  # Remove any counts that are exactly zero.
  ws <- sum(w)
  bs <- colSums(B)  
  i  <- which(w > 0)
  w  <- w[i]
  B  <- B[i,]

  # Run the SQP updates for the modified problem.
  out <- mixsqp(scale.cols(B,ws/bs),w/ws,y*bs/ws,numiter,e,tol,
                zero.threshold,zero.searchdir,suffdecr,stepsizereduce,
                minstepsize)

  # Recover the solution to the unmodified problem.
  return(out$x*ws/bs)
}
